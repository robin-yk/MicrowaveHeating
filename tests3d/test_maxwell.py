import numpy as np
from microwave3d.grid import Grid
from microwave3d.electromagnetics import Maxwell
from microwave3d.materials import EPS0,MU0


def box(n=4):
    return Grid(*(np.linspace(-length/2,length/2,n+1) for length in [.06,.045,.065]))


def props(grid,loss=.2):
    return {"ep":np.ones(grid.n),"epp":np.full(grid.n,loss),"bed_epp":np.full(grid.n,loss),"quartz_epp":np.zeros(grid.n)}


def test_rectangular_cavity_vector_eigenmode_and_second_order(cfg):
    errors=[]
    for n in (4,8,16):
        grid=box(n); cfg["em"]["wall_model"]="pec"
        model=Maxwell(grid,cfg)
        voltage=np.zeros(grid.ne)
        # TE101: Ey is normal to the y walls, zero tangentially at x,z walls.
        for i,j,k in np.ndindex(grid.edges[1].shape):
            voltage[grid.edges[1][i,j,k]]=grid.dx[1][j]*np.sin(np.pi*i/n)*np.sin(np.pi*k/n)
        dx,dz=grid.dx[0][0],grid.dx[2][0]
        wd2=((2*np.sin(np.pi/(2*n))/dx)**2+(2*np.sin(np.pi/(2*n))/dz)**2)/(EPS0*MU0)
        me=EPS0*grid.mass(np.ones(grid.n))*voltage
        free=grid.free_pec
        assert np.linalg.norm((model.K@voltage-wd2*me)[free])/np.linalg.norm((wd2*me)[free]) < 1e-12
        exact=np.sqrt((np.pi/.06)**2+(np.pi/.065)**2)/np.sqrt(EPS0*MU0)
        errors.append(abs(np.sqrt(wd2)/exact-1))
    assert errors[0]/errors[1]>3.8 and errors[1]/errors[2]>3.9


def test_passive_port_balance_and_power_scaling(cfg):
    grid=box(); model=Maxwell(grid,cfg)
    a=model.solve(props(grid),2.404e9,2)
    b=model.solve(props(grid),2.404e9,8)
    assert abs(a.power_residual)<1e-9
    assert a.sample_w>0 and a.wall_w>0
    assert 0 <= a.reflected_w <= a.incident_w
    np.testing.assert_allclose(b.q_sample,4*a.q_sample,rtol=1e-9,atol=1e-15)
    np.testing.assert_allclose(b.voltage,2*a.voltage,rtol=1e-9,atol=1e-15)
    assert abs(a.s11-b.s11)<1e-9


def test_lossless_pec_cavity_reflects_all_power(cfg):
    cfg["em"]["wall_model"]="pec"
    grid=box(); model=Maxwell(grid,cfg)
    result=model.solve(props(grid,0),2.404e9,3)
    assert result.sample_w==0 and result.wall_w==0
    assert abs(result.reflected_w-3)<1e-8


def test_heat_quadrature_equals_dielectric_loss(cfg):
    grid=box(); p=props(grid)
    p["epp"]=np.linspace(.01,.5,grid.n); p["bed_epp"]=p["epp"]
    e=Maxwell(grid,cfg).solve(p,2.4e9,2)
    direct=.5*2*np.pi*2.4e9*EPS0*np.dot(grid.mass(p["epp"]),abs(e.voltage)**2)
    assert abs(direct-e.sample_w)/e.sample_w < 1e-12


def test_zero_power_and_incoherent_dwell(cfg):
    grid=box(); model=Maxwell(grid,cfg)
    zero=model.solve(props(grid),2.4e9,0)
    assert np.max(zero.e2)==0 and zero.sample_w==0
    driven=model.solve(props(grid),2.4e9,1)
    assert abs(zero.s11-driven.s11)<1e-12
    runs,weights=model.spectrum(props(grid),[2.35e9,2.45e9],4,[1,3])
    assert sum(weights)==1
    accepted=sum(w*(e.sample_w+e.quartz_w+e.wall_w+e.reflected_w) for w,e in zip(weights,runs))
    assert abs(accepted-4)<1e-8


def test_aperture_and_ohmic_loss_are_in_power_balance(cfg):
    cfg["em"]["tube_apertures"]=True
    grid=box(); p=props(grid,.05); p["bed_sigma"]=np.full(grid.n,.1)
    result=Maxwell(grid,cfg).solve(p,2.4e9,4)
    assert result.aperture_w>0
    total=result.sample_w+result.quartz_w+result.wall_w+result.aperture_w+result.reflected_w
    assert abs(total-4)<1e-8
    ohmic=.5*np.dot(p["bed_sigma"]*result.e2,grid.volume)
    assert 0<ohmic<result.sample_w
