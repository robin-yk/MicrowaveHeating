# Microwave3D | paused research project

Development is paused. This directory preserves the vector Maxwell solver,
Darcy and square-channel Navier–Stokes options, conjugate heat transfer,
input provenance, numerical tests, and offline result viewer.
Flow assumes constant reference density/viscosity and equal channel flow.
Actual coax geometry, thermal expansion, and measured material curves remain
outside the implemented example. See the detailed methodology below.

## Independent setup

Run these commands from the repository root:

```sh
cd research/microwave3d
python3 -m venv .venv
.venv/bin/python -m pip install -e '.[test]'
.venv/bin/python -m pytest -q
.venv/bin/python -m microwave3d.cli --config examples/itaca-monolith-ns.json --output results3d/monolith
```

Open the generated `results3d/monolith/viewer.html`. The HTML inside the
`microwave3d/` package is a template, not a computed result or live solver.
No Python dependencies are needed to run the separate 2D browser application.

## Directory map

| Path | Purpose |
| --- | --- |
| `microwave3d/` | Python solver and viewer template |
| `examples/` | Input cases and their provenance |
| `tests/` | Python numerical regression tests |
| `tools/verification/` | Empty-cavity geometric benchmark |
| `docs/` | Equations, implementation scope, verification records |
| `pyproject.toml` | Independent Python package configuration |
| `requirements-tested.txt` | Versions used in the recorded checks |
| `results3d/` | New local outputs, ignored by Git |

- [Setup and equations](docs/MICROWAVE3D.md)
- [Channel flow and heat transfer](docs/MICROWAVE3D-CHANNELS.md)
- [Verification record](docs/MICROWAVE3D-VERIFICATION.md)
- [Input provenance](examples/README.md)

## Migration from the earlier root layout

The Python import name `microwave3d` is unchanged. Source files moved under this
project; `examples3d/` became `examples/`, and `tests3d/` became `tests/`.
Commands in this project's documents are relative to this directory.
Historical calculation values and input records have not been rewritten.

Existing exports at the repository root `results3d/` were not moved or deleted,
so previously opened result URLs continue to work. Old root-level source-template
URLs no longer identify the viewer template; use a generated result instead.
Original pre-move source paths remain available on `feat/microwave3d-emthermal`.

If reusing the earlier root virtual environment, update its editable install:

```sh
# From the repository root:
.venv/bin/python -m pip install -e './research/microwave3d[test]'
# Then run from the 3D project directory:
cd research/microwave3d
../../.venv/bin/python -m pytest -q
```
