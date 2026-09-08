"""Conservative conduction, flow enthalpy transport, and enclosure radiation."""
from dataclasses import dataclass
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve, spilu, gmres, LinearOperator
from .materials import SIGMA


def conduction(grid, conductivity, ambient, contact_fraction=None, contact_h=None, channel_fraction=None):
    n = grid.n; links = grid.links
    i,j = links[:,0].astype(int),links[:,1].astype(int)
    k = np.asarray(conductivity)
    resistance = links[:,4]/k[i]+links[:,5]/k[j]
    if contact_fraction is not None and contact_h is not None:
        # Contact is attached only to faces separating dominant bed and quartz.
        bed,quartz = contact_fraction
        contact = ((bed[i]>.5)&(quartz[j]>.5))|((quartz[i]>.5)&(bed[j]>.5))
        resistance = resistance + contact/contact_h
    G = links[:,3]/resistance
    if channel_fraction is not None:
        # Channel end planes are flow boundaries, not stagnant downstream cells.
        end=(links[:,2]==2)&((channel_fraction[i]>.5) != (channel_fraction[j]>.5))
        G[end]=0
    diag = np.bincount(np.r_[i,j],weights=np.r_[G,G],minlength=n)
    b = grid.boundary; cells = b[:,0].astype(int)
    gb = b[:,3]*k[cells]/b[:,4]
    np.add.at(diag,cells,gb)
    rhs = np.zeros(n); np.add.at(rhs,cells,gb*ambient)
    diag[~grid.active] = 1.; rhs[~grid.active] = ambient
    A = sparse.diags(diag)+sparse.coo_matrix((np.r_[-G,-G],(np.r_[i,j],np.r_[j,i])),shape=(n,n)).tocsr()
    return A,rhs,gb


def advection(grid, flow, inlet_temperature):
    """Upwind enthalpy flux on the exact same faces as the Darcy solve."""
    n = grid.n
    left,right = flow.links[:,0].astype(int),flow.links[:,1].astype(int)
    flux = flow.mass_flux*flow.cp
    upstream = np.where(flux>=0,left,right)
    downstream = np.where(flux>=0,right,left)
    g = abs(flux)
    A = sparse.coo_matrix((np.r_[g,-g,flow.outlet_mass*flow.cp],
         (np.r_[upstream,downstream,flow.outlet_cells],np.r_[upstream,upstream,flow.outlet_cells])),shape=(n,n)).tocsr()
    rhs = np.zeros(n); np.add.at(rhs,flow.inlet_cells,flow.inlet_mass*flow.cp*inlet_temperature)
    return A,rhs


def radiation(temperature, ambient, area, cavity_area, tube_emissivity, wall_emissivity):
    """Diffuse-gray convex tube inside an isothermal gray enclosure.

    Tube patches see only the enclosure. The wall radiosity is uniform. This
    closure includes wall reflection, but excludes quartz spectral transmission.
    """
    tk = np.asarray(temperature)+273.15; ta = ambient+273.15
    total = area.sum()
    if total == 0 or tube_emissivity == 0:
        return np.zeros_like(tk),np.zeros_like(tk),np.zeros_like(tk)
    f = total/cavity_area
    if f > 1:
        raise ValueError("Tube area exceeds the enclosing cavity area")
    mean_black = np.dot(area,SIGMA*tk**4)/total
    irradiation = (wall_emissivity*SIGMA*ta**4+(1-wall_emissivity)*f*tube_emissivity*mean_black)/(wall_emissivity+(1-wall_emissivity)*f*tube_emissivity)
    heat = tube_emissivity*area*(SIGMA*tk**4-irradiation)
    secant = tube_emissivity*area*SIGMA*(tk+ta)*(tk*tk+ta*ta)
    rhs = secant*ambient+tube_emissivity*area*(irradiation-SIGMA*ta**4)
    return heat,secant,rhs


@dataclass
class ThermalResult:
    temperature: np.ndarray
    boundary_w: float
    radiation_w: float
    gas_w: float
    residual: float
    max_cell_residual_w: float


class Thermal:
    def __init__(self, grid, materials, cfg):
        self.grid,self.materials,self.cfg = grid,materials,cfg
        self._preconditioner = None

    def operators(self, temperature, properties, flow):
        cfg,m,g = self.cfg,self.materials,self.grid
        ambient = cfg["thermal"]["ambient_c"]
        A,rhs,gb = conduction(g,properties["k"],ambient,
            (m.bed,m.quartz),cfg["thermal"].get("contact_h_w_m2k"),m.fluid_fraction if cfg.get("channels") and flow is not None else None)
        if cfg.get("channels") and flow is not None:
            cells=flow.inlet_cells
            gi=2*properties["k"][cells]*g.volume[cells]/g.cell_lengths[cells,2]**2
            diag=np.zeros(g.n); np.add.at(diag,cells,gi)
            A+=sparse.diags(diag); rhs+=diag*cfg["gas"]["inlet_c"]
        if flow is not None:
            adv,adv_rhs = advection(g,flow,cfg["gas"]["inlet_c"])
            A += adv; rhs += adv_rhs
        qr,hr,rr = radiation(temperature,ambient,m.tube_area,m.cavity_area,
            cfg["quartz"]["thermal_emissivity"],cfg["thermal"]["wall_emissivity"])
        return A+sparse.diags(hr),rhs+rr,gb,qr

    def step(self, temperature, properties, flow, source_w_m3):
        A,rhs,_,_ = self.operators(temperature,properties,flow)
        rhs = rhs+source_w_m3*self.grid.volume
        if self.cfg["solver"].get("thermal_linear_solver","direct") == "krylov":
            scale=1/np.sqrt(A.diagonal())
            D=sparse.diags(scale); balanced=(D@A@D).tocsc(); b=scale*rhs
            for attempt in range(2):
                if self._preconditioner is None or attempt:
                    ilu=spilu(balanced,drop_tol=1e-4,fill_factor=12)
                    self._preconditioner=LinearOperator(balanced.shape,ilu.solve)
                y,info=gmres(balanced,b,x0=temperature/scale,M=self._preconditioner,
                    rtol=1e-12,atol=1e-14,restart=60,maxiter=30)
                t=scale*y
                if info == 0 and np.max(abs(A@t-rhs)) < 1e-10*max(np.max(abs(rhs)),1.):
                    break
            else:
                raise RuntimeError("Preconditioned thermal solve failed its explicit residual check")
        else:
            t = spsolve(A.tocsc(),rhs)
        if not np.all(np.isfinite(t)):
            raise RuntimeError("Thermal linear solve returned nonfinite temperatures")
        return t

    def audit(self, temperature, properties, flow, source_w_m3):
        g=self.grid
        A,rhs,gb,qr = self.operators(temperature,properties,flow)
        cell_residual = A@temperature-rhs-source_w_m3*g.volume
        qb = float(np.dot(gb,temperature[g.boundary[:,0].astype(int)]-self.cfg["thermal"]["ambient_c"]))
        if self.cfg.get("channels") and flow is not None:
            cells=flow.inlet_cells
            gi=2*properties["k"][cells]*g.volume[cells]/g.cell_lengths[cells,2]**2
            qb+=float(gi@(temperature[cells]-self.cfg["gas"]["inlet_c"]))
        qg = float(flow.cp*(flow.outlet_mass@temperature[flow.outlet_cells]-flow.inlet_mass.sum()*self.cfg["gas"]["inlet_c"])) if flow is not None else 0.
        deposited = float(source_w_m3@g.volume)
        residual = (deposited-qb-qr.sum()-qg)/max(deposited,1e-9)
        return ThermalResult(temperature,qb,float(qr.sum()),qg,float(residual),float(np.max(abs(cell_residual[g.active]))))
