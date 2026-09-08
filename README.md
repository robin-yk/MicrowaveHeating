# Microwave Heating Models

## Start here

| What you need | Where to go |
| --- | --- |
| Use the 2D browser model | [apps/microwave/index.html](apps/microwave/index.html) |
| Understand or change 2D | [2D project guide](apps/microwave/README.md) |
| Inspect the paused 3D work | [3D research guide](research/microwave3d/README.md) |
| Read 2D model limitations | [2D limits](apps/microwave/docs/LIMITS.md) |

## Repository structure

Each project owns its code, tests, tools, and documentation.

```text
apps/
  microwave/          2D browser application
    index.html        Interface
    solver.js         Numerical model
    tests/            Regression tests
    tools/            Verification and note reproduction
    docs/             Equations, limits, and verification records
    assets/           Project images
research/
  microwave3d/        Paused 3D research project
    microwave3d/      Python package and viewer template
    examples/         Input cases
    tests/            Numerical tests
    tools/            Verification scripts
    docs/             Methodology and verification records
index.html            Existing redirect to the 2D application
package.json          Root shortcuts to the 2D workspace
```

## Run the 2D application

The existing root and `apps/microwave/` page addresses are unchanged.
For local development, from the repository root:

```sh
npm install
npm run dev
npm test
npm run build
```

Root commands forward to the 2D project. Its build output is
`apps/microwave/dist/`. Verification commands remain `npm run verify:microwave`,
`npm run verify:field`, `npm run verify:calibrate`, and `npm run si:note`.
The 3D project has its own Python setup; see its guide rather than these npm commands.

## Development status and preserved work

- **2D:** existing browser application; model limitations are documented in its project.
- **3D:** development paused. Source and numerical records are preserved on main,
  separately from 2D. Merging the code does not resume development.
- **Single-screen UI experiment:** remains on the separate
  [claude/microwave-simple-ui-ujh2kh branch](https://github.com/robin-yk/microwave-2D/tree/claude/microwave-simple-ui-ujh2kh), not merged.
- **Earlier 3D layout:** retained on `feat/microwave3d-emthermal`.

Existing local `results3d/` exports and `.venv/` were not deleted or relocated.
They are ignored artifacts, not source projects. Generated HTML results keep
working at their existing paths. New project outputs stay inside their project
and are excluded from Git.

## History and license

Split from [Electrification Suite](https://github.com/robin-yk/Electrification-Suite),
which retained the Joule and pulsed-heating tools. The microwave solver's history
is preserved. By Yeonsu Kwak (Vlachos Lab, University of Delaware).
[MIT license](LICENSE).
