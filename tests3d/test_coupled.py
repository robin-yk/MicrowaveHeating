import numpy as np
import pytest
from microwave3d.solver import Coupled,validate


def small(cfg):
    cfg["mesh"].update(core_xy=4,outer_xy=2,core_z=4,outer_z=2)
    return cfg


def test_zero_power_is_ambient(cfg):
    model=Coupled(small(cfg))
    summary,fields=model.solve(0)
    assert summary["converged"]
    np.testing.assert_allclose(fields["temperature_c"],20,atol=1e-8)


def test_coupled_solution_closes_with_fresh_maxwell(cfg):
    model=Coupled(small(cfg))
    summary,fields=model.solve(3)
    assert summary["converged"],summary["history"][-1]
    assert abs(summary["residuals"]["global_relative"])<1e-4
    assert summary["temperature_c"]["sample_max"]>20
    assert summary["power_w"]["sample"]<3
    p,em,w,q,flow,audit=model.evaluate(fields["temperature_c"],3)
    assert abs(audit.residual)<cfg["solver"]["power_tolerance"]
    assert np.max(abs(fields["q_sample_w_m3"]-em[0].q_sample))<1e-10
    assert flow.mass_residual<1e-9


def test_nonconvergence_is_not_success(cfg):
    cfg=small(cfg);cfg["solver"]["max_iterations"]=1
    summary,_=Coupled(cfg).solve(3)
    assert not summary["converged"]


def test_bad_inputs_and_temperature_extrapolation_fail(cfg):
    cfg=small(cfg); cfg["material"]["porosity"]=1
    with pytest.raises(ValueError):validate(cfg)
    cfg["material"]["porosity"]=.6
    model=Coupled(cfg)
    with pytest.raises(ValueError,match="outside material"):
        model.materials.at(np.full(model.grid.n,2000.))
