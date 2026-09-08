# Microwave Heating 2D Model

## New: local 3D electromagnetic and thermal research solver

The `microwave3d/` Python package adds a full vector Maxwell calculation,
temperature-dependent heat deposition, a 3D thermal field, and optional coupled
packed-bed Darcy flow or resolved square-channel Navier–Stokes flow. It includes incident-power/frequency sweeps and an offline
3D field viewer. The ITACA cylindrical cavity envelope is based on the supplied
Malhotra paper; the probe and material example remain assumptions. This is a
numerically tested research prototype, not an experimentally validated reactor model.

For the 21-channel example, run `python -m microwave3d.cli --config examples3d/itaca-monolith-ns.json --output results3d/monolith-ns`.
The channel solver retains 3D inertia and viscosity, uses no-slip walls, and
couples conservative gas enthalpy transport to a separate solid/gas thermal mesh.
It assumes constant flow density/viscosity and equal channel flow; manifolds and
thermal expansion are not resolved. See the [channel methodology](docs/MICROWAVE3D-CHANNELS.md).

See [setup, equations, and limitations](docs/MICROWAVE3D.md) and
[input provenance](examples3d/README.md). The existing 2D application below is
retained independently.

## Existing 2D application

Steady 2D temperature field of a microwave-heated powder bed, with a dielectric
response that follows the local temperature and a Helmholtz field solve, run
entirely in the browser. By Yeonsu Kwak (Vlachos Lab, University of Delaware).

Open `apps/microwave/index.html` in a browser. There is no build step and
nothing to install.

## Solver API

The numeric core is a dependency-free, DOM-free ES module, importable in Node
or in another page:

```js
import { solve2D, transportNumbers, materialProfiles } from "./apps/microwave/solver.js";
```

`tests/microwave-solver.test.js` doubles as a worked example of its inputs and
outputs.

## Testing and verification

```bash
npm test                  # Node regression suite
npm run verify:microwave  # grid convergence and the analytic benchmarks
npm run verify:field      # the Helmholtz field solve
npm run verify:calibrate  # what the calibration is actually fitting
```

`docs/VERIFICATION.md` records what each study found, including the grid
convergence the gas-exchange term breaks and the finding that the published
calibration is partly fitting mesh error. `docs/LIMITS.md` states what the
model does not do.

## History

This repository was split out of
[Electrification Suite](https://github.com/robin-yk/Electrification-Suite),
which kept the Joule heating and pulsed heating tools. The commit history of
the microwave solver came across with it.

## License

MIT (c) Yeonsu Kwak
