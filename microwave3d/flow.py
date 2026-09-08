"""Conservative 3D Darcy mass flux in a homogenized packed bed.

The solved mass flux enters thermal advection. This is porous CFD, not a
Navier-Stokes calculation of monolith channels or external natural convection.
"""
from dataclasses import dataclass
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from .materials import gas_properties


@dataclass
class FlowResult:
    pressure_pa: np.ndarray
    links: np.ndarray
    mass_flux: np.ndarray
    inlet_cells: np.ndarray
    inlet_mass: np.ndarray
    outlet_cells: np.ndarray
    outlet_mass: np.ndarray
    cp: float
    pressure_drop_pa: float
    mass_residual: float
    reynolds_particle_max: float


def darcy(grid, fraction, temperature, cfg):
    gas, mat = cfg["gas"],cfg["material"]
    rho,mu,cp,mdot = gas_properties(temperature,gas)
    ids = np.flatnonzero((fraction>0)&grid.active)
    if not len(ids):
        raise ValueError("No sample cells for porous flow")
    numbering = np.full(grid.n,-1); numbering[ids] = np.arange(len(ids))
    eps, dp = mat["porosity"],mat["particle_diameter_m"]
    permeability = eps**3*dp**2/(150*(1-eps)**2)
    mobility = rho*permeability/mu
    axes = np.asarray(mat.get("permeability_ratios",[1.,1.,1.]))
    links = grid.links
    left,right = links[:,0].astype(int),links[:,1].astype(int)
    take = (fraction[left]>0)&(fraction[right]>0)
    links = links[take]; left,right = left[take],right[take]
    area = links[:,3]*np.minimum(fraction[left],fraction[right])
    conductance = area/(links[:,4]/mobility[left]+links[:,5]/mobility[right])*axes[links[:,2].astype(int)]
    # Top and bottom bed faces are isobars. Closed side faces carry no mass.
    zz = grid.xyz[ids,2]; top,bottom = zz.max(),zz.min()
    inlet = ids[np.isclose(zz,top,rtol=0,atol=1e-12)]
    outlet = ids[np.isclose(zz,bottom,rtol=0,atol=1e-12)]
    end_g = lambda cells: (mobility[cells]*axes[2]*fraction[cells]*grid.volume[cells]
                          /grid.cell_lengths[cells,2]/(grid.cell_lengths[cells,2]/2))
    gin,gout = end_g(inlet),end_g(outlet)
    il,ir = numbering[left],numbering[right]
    diagonal = np.bincount(np.r_[il,ir],weights=np.r_[conductance,conductance],minlength=len(ids))
    np.add.at(diagonal,numbering[inlet],gin); np.add.at(diagonal,numbering[outlet],gout)
    A = sparse.diags(diagonal)+sparse.coo_matrix((np.r_[-conductance,-conductance],
               (np.r_[il,ir],np.r_[ir,il])),shape=(len(ids),len(ids))).tocsr()
    rhs = np.zeros(len(ids)); np.add.at(rhs,numbering[inlet],gin)
    phi = spsolve(A.tocsc(),rhs)
    unit_rate = np.sum(gin*(1-phi[numbering[inlet]]))
    if unit_rate <= 0:
        raise RuntimeError("Darcy unit-flow conductance is not positive")
    drop = mdot/unit_rate
    pressure = np.zeros(grid.n); pressure[ids] = drop*phi
    flux = conductance*(pressure[left]-pressure[right])
    fin = gin*(drop-pressure[inlet]); fout = gout*pressure[outlet]
    residual = np.zeros(grid.n)
    np.add.at(residual,left,flux); np.add.at(residual,right,-flux)
    np.add.at(residual,inlet,-fin); np.add.at(residual,outlet,fout)
    err = float(np.max(abs(residual))/max(mdot,1e-30))
    if err > 1e-7 or drop > .01*gas["pressure_pa"]:
        raise RuntimeError("Darcy flow failed conservation or ΔP/P exceeds the 1% constant-pressure density limit")
    rep = float(np.max(abs(flux)/np.maximum(area,1e-30)*dp/((mu[left]+mu[right])/2))) if len(flux) else 0.
    if rep > 1:
        raise ValueError("Particle Reynolds number exceeds the creeping porous-flow scope")
    return FlowResult(pressure,links,flux,inlet,fin,outlet,fout,cp,float(drop),err,rep)


def cell_velocity(grid, fraction, temperature, cfg, flow):
    """Area-weighted reconstruction of superficial velocity for display/export."""
    rho,_,_,_=gas_properties(temperature,cfg["gas"])
    velocity=np.zeros((grid.n,3)); weights=np.zeros((grid.n,3))
    for link,flux in zip(flow.links,flow.mass_flux):
        i,j,axis=map(int,link[:3]); area=link[3]*min(fraction[i],fraction[j])
        velocity[i,axis] += flux/rho[i]; velocity[j,axis] += flux/rho[j]
        weights[i,axis] += area; weights[j,axis] += area
    for cells,mass in ((flow.inlet_cells,flow.inlet_mass),(flow.outlet_cells,flow.outlet_mass)):
        area=fraction[cells]*grid.volume[cells]/grid.cell_lengths[cells,2]
        velocity[cells,2] -= mass/rho[cells]; weights[cells,2] += area
    return np.divide(velocity,weights,out=np.zeros_like(velocity),where=weights>0)
