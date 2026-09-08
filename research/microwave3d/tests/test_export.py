import json
import numpy as np
from microwave3d.solver import Coupled
from microwave3d.cli import save,main


def test_exports_include_real_fields_and_offline_viewer(cfg,tmp_path):
    cfg["mesh"].update(core_xy=4,outer_xy=2,core_z=4,outer_z=2)
    model=Coupled(cfg);summary,fields=model.solve(0)
    save(tmp_path,summary,fields)
    assert json.loads((tmp_path/"result.json").read_text())["converged"]
    arrays=np.load(tmp_path/"fields.npz")
    np.testing.assert_array_equal(arrays["temperature_c"],fields["temperature_c"])
    assert "VECTORS superficial_velocity_m_s" in (tmp_path/"fields.vtk").read_text()
    html=(tmp_path/"viewer.html").read_text()
    assert "window.RESULT=" in html and "/*__RESULT_DATA__*/" not in html
    assert "experimental validation pending" in html


def test_cli_returns_failure_for_unconverged_run(cfg,tmp_path):
    cfg["mesh"].update(core_xy=4,outer_xy=2,core_z=4,outer_z=2)
    cfg["solver"]["max_iterations"]=1
    config=tmp_path/"case.json";config.write_text(json.dumps(cfg))
    assert main(["--config",str(config),"--power","5","--output",str(tmp_path/"out")])==2
    assert not json.loads((tmp_path/"out"/"result.json").read_text())["converged"]
