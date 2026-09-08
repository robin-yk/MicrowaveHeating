"""Three-dimensional steady laminar Navier–Stokes in a square channel.

Taylor–Hood P2 velocity/P1 pressure tetrahedra; the convective term is retained.
The solved reference gas is incompressible with constant rho and mu. This is
not a variable-density low-Mach solver. Thermal expansion is not represented.
"""
from dataclasses import dataclass
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from skfem import MeshTet, Basis, FacetBasis, ElementVector, ElementTetP2, ElementTetP1, BilinearForm, Functional, asm
from skfem.helpers import grad, ddot, div, dot


@BilinearForm
def viscous(u,v,w):
    return ddot(grad(u),grad(v))


@BilinearForm
def convection(u,v,w):
    return dot(np.einsum('j...,ij...->i...',w.wind,grad(u)),v)


@BilinearForm
def continuity(u,q,w):
    return -div(u)*q


def square_profile(x,y):
    """Fully developed square-duct profile normalized to mean velocity one."""
    value=np.zeros_like(x)
    mean=0.
    for n in range(1,100,2):
        # Stable cosh(n*pi*(x-.5))/cosh(n*pi/2).
        ratio=(np.exp(-n*np.pi*x)+np.exp(-n*np.pi*(1-x)))/(1+np.exp(-n*np.pi))
        value += (1-ratio)*np.sin(n*np.pi*y)/n**3
        mean += 2/(n**4*np.pi)*(1-2*np.tanh(n*np.pi/2)/(n*np.pi))
    return value/mean


@dataclass
class ChannelSolution:
    velocity_basis: object
    pressure_basis: object
    velocity: np.ndarray
    pressure: np.ndarray
    width: float
    mean_speed: float
    pressure_scale: float
    diagnostics: dict

    def sample(self,points):
        # points are local physical coordinates, z increases from inlet.
        p=np.asarray(points).T/self.width
        u=np.asarray(self.velocity_basis.probes(p)@self.velocity).reshape(3,-1).T*self.mean_speed
        pressure=np.asarray(self.pressure_basis.probes(p)@self.pressure)*self.pressure_scale
        return u,pressure


def solve_channel(width,length,mean_speed,rho,mu,cross=3,axial=12,
                  tolerance=1e-8,max_iterations=60,fully_developed=False):
    if not np.all(np.isfinite([width,length,mean_speed,rho,mu,tolerance])) or min(width,length,rho,mu,tolerance) <= 0 or mean_speed < 0 or cross < 2 or axial < 2:
        raise ValueError("Invalid channel flow dimensions, properties, or mesh")
    re=rho*mean_speed*width/mu
    if re > 200:
        raise ValueError("Current unstabilized laminar channel solver is restricted to Re <= 200")
    aspect=length/width
    mesh=MeshTet.init_tensor(np.linspace(0,1,cross+1),np.linspace(0,1,cross+1),np.linspace(0,aspect,axial+1))
    vb=Basis(mesh,ElementVector(ElementTetP2()),intorder=4)
    pb=Basis(mesh,ElementTetP1(),intorder=4)
    wall=mesh.facets_satisfying(lambda x: np.isclose(x[0],0)|np.isclose(x[0],1)|np.isclose(x[1],0)|np.isclose(x[1],1),boundaries_only=True)
    inlet=mesh.facets_satisfying(lambda x: np.isclose(x[2],0),boundaries_only=True)
    outlet=mesh.facets_satisfying(lambda x: np.isclose(x[2],aspect),boundaries_only=True)
    fixed=np.unique(np.r_[vb.get_dofs(wall).all(),vb.get_dofs(inlet).all(),vb.get_dofs(outlet).all() if fully_developed else []]).astype(int)
    values=np.zeros(vb.N+pb.N)
    for facets in ([inlet,outlet] if fully_developed else [inlet]):
        dofs=vb.get_dofs(facets).all(['u^3'])
        x,y=vb.doflocs[:2,dofs]
        values[dofs]=square_profile(x,y) if fully_developed else 36*x*(1-x)*y*(1-y)
    values[vb.get_dofs(wall).all()]=0
    @Functional
    def throughput(w):
        return w.velocity[2]
    fb=FacetBasis(mesh,ElementVector(ElementTetP2()),facets=inlet,intorder=5)
    inlet_integral=asm(throughput,fb,velocity=fb.interpolate(values[:vb.N]))
    values[:vb.N]/=inlet_integral
    if fully_developed:
        gauge=pb.get_dofs(outlet).all()[0]
        fixed=np.r_[fixed,vb.N+gauge]
    free=np.setdiff1d(np.arange(len(values)),fixed)
    K=asm(viscous,vb); B=asm(continuity,vb,pb)
    def operator(u):
        A=K+re*asm(convection,vb,wind=vb.interpolate(u))
        return sparse.bmat([[A,B.T],[B,None]],format='csc')
    state=values.copy(); history=[]
    for it in range(max_iterations):
        A=operator(state[:vb.N])
        candidate=values.copy()
        candidate[free]=spsolve(A[free][:,free],-(A@values)[free])
        if not np.all(np.isfinite(candidate)):
            raise RuntimeError("Navier–Stokes linear system failed")
        nonlinear=operator(candidate[:vb.N])@candidate
        residual=float(np.linalg.norm(nonlinear[free])/max(np.linalg.norm((K@candidate[:vb.N])),1e-15))
        change=float(np.max(abs(candidate[:vb.N]-state[:vb.N])))
        history.append({'iteration':it+1,'relative_residual':residual,'velocity_step':change})
        state=candidate if re<20 else .7*candidate+.3*state
        if residual<tolerance and change<tolerance:
            state=candidate
            break
    else:
        raise RuntimeError(f"Navier–Stokes did not converge: residual={residual:g}")
    divergence=B@state[:vb.N]
    @Functional
    def pressure_integral(w):
        return w.pressure
    end_pressure=[]
    for facets in (inlet,outlet):
        fp=FacetBasis(mesh,ElementTetP1(),facets=facets,intorder=4)
        end_pressure.append(float(asm(pressure_integral,fp,pressure=fp.interpolate(state[vb.N:]))))
    diagnostics={'model':'3D incompressible Navier-Stokes (P2/P1)',
        'reynolds_hydraulic':float(re),'iterations':len(history),'relative_nonlinear_residual':residual,
        'weak_continuity_l2':float(np.linalg.norm(divergence)),
        'velocity_unknowns':int(vb.N),'pressure_unknowns':int(pb.N),'tetrahedra':int(mesh.nelements),
        'rho_kg_m3':float(rho),'mu_pa_s':float(mu),'mean_speed_m_s':float(mean_speed),
        'pressure_drop_pa':(end_pressure[0]-end_pressure[1])*mu*mean_speed/width,
        'outlet_condition':'velocity profile + pressure gauge' if fully_developed else 'zero natural traction (grad-u formulation)',
        'history':history}
    return ChannelSolution(vb,pb,state[:vb.N],state[vb.N:],width,mean_speed,mu*mean_speed/width,diagnostics)
