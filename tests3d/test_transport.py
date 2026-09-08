import numpy as np
from scipy.sparse.linalg import spsolve
from microwave3d.grid import Grid
from microwave3d.thermal import conduction,advection,radiation
from microwave3d.flow import darcy
from microwave3d.materials import gas_properties,SIGMA


def test_manufactured_3d_conduction_second_order():
    errors=[]
    for n in (6,12,24):
        grid=Grid(*(np.linspace(0,.02,n+1) for _ in range(3)))
        exact=20+30*np.prod(np.sin(np.pi*grid.xyz/.02),axis=1)
        q=3*(np.pi/.02)**2*5*(exact-20)
        A,rhs,gb=conduction(grid,np.full(grid.n,5.),20.)
        t=spsolve(A.tocsc(),rhs+q*grid.volume)
        errors.append(np.sqrt(np.mean((t-exact)**2)))
        qb=gb@(t[grid.boundary[:,0].astype(int)]-20)
        assert abs(qb-q@grid.volume)/(q@grid.volume)<1e-11
    assert errors[0]/errors[1]>3.8 and errors[1]/errors[2]>3.9


def test_darcy_analytic_pressure_drop_and_cell_mass_balance(cfg):
    grid=Grid(*(np.linspace(-.005,.005,6) for _ in range(3)))
    temperature=np.full(grid.n,200.)
    flow=darcy(grid,np.ones(grid.n),temperature,cfg)
    rho,mu,cp,mdot=gas_properties(temperature,cfg["gas"])
    por=cfg["material"]["porosity"]; dp=cfg["material"]["particle_diameter_m"]
    K=por**3*dp**2/(150*(1-por)**2)
    expected=mdot*mu[0]*.01/(rho[0]*K*.01**2)
    assert abs(flow.pressure_drop_pa/expected-1)<1e-12
    assert flow.mass_residual<1e-12
    A,rhs=advection(grid,flow,200.)
    assert np.max(abs(A@temperature-rhs))<1e-14


def test_advection_source_matches_outlet_enthalpy(cfg):
    grid=Grid(*(np.linspace(-.005,.005,5) for _ in range(3)))
    flow=darcy(grid,np.ones(grid.n),np.full(grid.n,20.),cfg)
    A,rhs=advection(grid,flow,20.)
    source=np.full(grid.n,.1/grid.n)
    t=spsolve(A.tocsc(),rhs+source)
    gas=flow.cp*(flow.outlet_mass@t[flow.outlet_cells]-flow.inlet_mass.sum()*20)
    assert abs(gas-.1)<1e-12
    # No artificial radial gas equilibration: independent streamlines keep their
    # temperatures when their heat sources differ.
    source[grid.xyz[:,0]<0]*=2
    t=spsolve(A.tocsc(),rhs+source)
    out=flow.outlet_cells
    assert t[out[grid.xyz[out,0]<0]].mean()>t[out[grid.xyz[out,0]>0]].mean()


def test_radiation_gray_enclosure_matches_two_surface_formula():
    area=np.array([.01,.02,.03]); cavity=.8; et=.85; ew=.2
    q,h,rhs=radiation(np.full(3,500.),20.,area,cavity,et,ew)
    exact=SIGMA*((500+273.15)**4-(20+273.15)**4)/(1/(et*area.sum())+(1-ew)/(ew*cavity))
    assert abs(q.sum()/exact-1)<1e-12
    np.testing.assert_allclose(h*500-rhs,q,rtol=1e-12)
    cold,_,_=radiation(np.full(3,20.),20.,area,cavity,et,ew)
    assert np.max(abs(cold))<1e-12
