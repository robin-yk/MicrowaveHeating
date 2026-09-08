import json
from pathlib import Path
import numpy as np
import pytest
from microwave3d.navier_stokes import solve_channel
from microwave3d.channels import NestedMesh,geometry
from microwave3d.grid import reactor_grid
from microwave3d.materials import Materials
from microwave3d.materials import EPS0
from microwave3d.solver import Coupled


def monolith():
    cfg=json.loads((Path(__file__).parents[1]/'examples3d/itaca-monolith-ns.json').read_text())
    cfg['mesh'].update(core_xy=4,outer_xy=2,core_z=4,outer_z=2)
    cfg['channels'].update(centres_m=[[0,0]],thermal_cells_across=2,thermal_cells_axial=6,flow_cells_across=2,flow_cells_axial=6)
    cfg['gas']['flow_sccm']=5
    return cfg


def test_square_duct_pressure_converges_to_analytic():
    a=.0013; length=.015; speed=.05; mu=1.76e-5
    expected=28.45415377*mu*length*speed/a**2
    errors=[]
    for n in (2,3,4):
        s=solve_channel(a,length,speed,1.16,mu,cross=n,axial=6,fully_developed=True)
        _,p=s.sample([[a/2,a/2,0],[a/2,a/2,length]])
        errors.append(abs((p[0]-p[1])/expected-1))
        assert s.diagnostics['relative_nonlinear_residual']<1e-8
        assert s.diagnostics['weak_continuity_l2']<1e-8
        u,_=s.sample([[0,a/2,length/2],[a,a/2,length/2]])
        np.testing.assert_allclose(u,0,atol=1e-10)
    assert errors[2]<errors[1]<errors[0]
    assert errors[-1]<.006


def test_nested_transfer_conserves_phase_watts():
    cfg=monolith(); em=reactor_grid(cfg); transfer=NestedMesh(em,cfg)
    m=Materials(transfer.grid,cfg)
    phase=m.bed
    coarse=transfer.average(phase)*123.4
    fine=transfer.deposit(coarse,phase)
    assert fine@transfer.grid.volume == pytest.approx(coarse@em.volume,rel=1e-12)
    assert np.all(fine[phase==0]==0)
    channel_volume=m.fluid_fraction@transfer.grid.volume
    assert channel_volume == pytest.approx(cfg['channels']['width_m']**2*cfg['geometry']['sample_length_m'],rel=1e-10)
    np.testing.assert_allclose(transfer.average(np.ones(transfer.grid.n)),1,rtol=1e-9)


def test_channel_coupling_energy_mass_and_zero_power():
    model=Coupled(monolith())
    zero,_=model.solve(0)
    assert zero['converged']
    summary,fields=model.solve(3)
    assert summary['converged']
    assert abs(summary['residuals']['thermal_relative'])<1e-4
    assert summary['residuals']['flow_mass_relative']<1e-10
    assert summary['flow']['diagnostics']['model'].startswith('3D incompressible Navier-Stokes')
    assert summary['power_w']['gas_enthalpy']>0
    assert summary['temperature_c']['sample_max']>summary['temperature_c']['channel_gas_outlet_mass_average']>20
    assert fields['q_sample_w_m3']@model.grid.volume == pytest.approx(summary['power_w']['sample'],rel=1e-10)
    assert fields['q_quartz_w_m3']@model.grid.volume == pytest.approx(summary['power_w']['quartz'],rel=1e-10)
    assert np.all(fields['q_sample_w_m3'][fields['channel_fraction']>.999]==0)


def test_invalid_channel_geometry_rejected():
    cfg=monolith(); cfg['channels']['centres_m']=[[0,0],[0,0]]
    with pytest.raises(ValueError,match='overlap'):geometry(cfg)
    cfg['channels']['centres_m']=[[.005,0]]
    with pytest.raises(ValueError,match='corners'):geometry(cfg)


def test_krylov_matches_direct_thermal_solution():
    cfg=monolith(); cfg['solver']['thermal_linear_solver']='direct'
    model=Coupled(cfg)
    t=np.full(model.grid.n,20.)
    p,e,w,q,flow,a=model.evaluate(t,3)
    direct=model.thermal.step(t,p,flow,q)
    model.cfg['solver']['thermal_linear_solver']='krylov'
    iterative=model.thermal.step(t,p,flow,q)
    np.testing.assert_allclose(iterative,direct,atol=1e-7,rtol=1e-9)


def test_subcell_heat_tracks_local_temperature_dependent_loss():
    cfg=monolith(); cfg['gas']['enabled']=False
    cfg['material']['table']=[[20,9.72,1.,50],[1200,10.,4.,50]]
    model=Coupled(cfg)
    t=100+1000*model.grid.xyz[:,2]
    p,e,w,q,flow,a=model.evaluate(t,3)
    sample=model.thermal_source(e,w,'sample',p)
    expected=.5*2*np.pi*e[0].frequency_hz*EPS0*p['bed_epp']*model.display_em(e[0].e2)
    np.testing.assert_allclose(sample,expected,rtol=1e-10,atol=1e-9)
