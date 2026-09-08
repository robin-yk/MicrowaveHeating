# Microwave2D

Steady 2D temperature field of a microwave-heated powder bed, with a dielectric
response that follows the local temperature and a Helmholtz field solve, run
entirely in the browser. By Yeonsu Kwak (Vlachos Lab, University of Delaware).

Open `index.html` in a browser. There is no build step and
nothing to install.

## Solver API

The numeric core is a dependency-free, DOM-free ES module, importable in Node
or in another page:

```js
import { solve2D, transportNumbers, materialProfiles } from "./solver.js";
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


## Project layout

All paths in this README and the technical documents are relative to `apps/microwave/`.

| Path | Contents |
| --- | --- |
| `index.html` | Browser interface |
| `solver.js` | Numerical model; independent of the interface |
| `tests/` | Node regression tests |
| `tools/verification/` | Grid, field, and calibration studies |
| `tools/si/` | Supplementary-note reproduction scripts |
| `docs/` | Limits, verification records, and supplementary note |
| `assets/` | Project image assets |
| `citation.ris` | Citation download |
| `package.json`, `vite.config.js` | App commands and optional build setup |

From this directory, `npm test` runs the tests. From the repository root,
`npm test` forwards to this same project. Optional browser development uses
`npm install` at the repository root, followed by `npm run dev`.
`npm run build` writes this project's `dist/` directory.

The page and solver URLs remain `apps/microwave/index.html` and
`apps/microwave/solver.js`. Earlier root-level `tests/`, `tools/`, `docs/`, and
`assets/` now live here. Scientific calculations and recorded results are unchanged.

[Repository overview](../../README.md) | [MIT license](../../LICENSE)
