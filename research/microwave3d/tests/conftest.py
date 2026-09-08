import json
from pathlib import Path
import pytest


@pytest.fixture
def cfg():
    cfg=json.loads((Path(__file__).parents[1]/"examples/itaca-cylinder.json").read_text())
    cfg["em"]["tube_apertures"]=False
    return cfg
