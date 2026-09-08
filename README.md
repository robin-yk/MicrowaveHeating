# Microwave Heating Models

The main browser application is **Microwave2D**. **Microwave3D is paused** and
preserved as a separate research project; merging its history into main does
not replace the 2D application or resume 3D development.

| Project | Location | Status / entry point |
| --- | --- | --- |
| Microwave2D | [apps/microwave/](apps/microwave/) | Existing browser application; open index.html |
| Microwave3D | [research/microwave3d/](research/microwave3d/) | Paused Python research solver with its own examples, tests, and docs |
| Single-screen UI experiment | [claude/microwave-simple-ui-ujh2kh](https://github.com/robin-yk/microwave-2D/tree/claude/microwave-simple-ui-ujh2kh) | Separate branch; not merged into main |

Root `package.json`, `tests/`, `tools/`, `docs/`, and `assets/` belong to 2D.
Python configuration and all tracked 3D files live under `research/microwave3d/`.
Existing ignored `results3d/` exports and the local `.venv/` are preserved at the
repository root. They are local artifacts, not another source project.

See the [3D run instructions and path migration](research/microwave3d/README.md).
The root browser entry still opens the existing 2D application.

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
