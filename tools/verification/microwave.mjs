// Numerical verification of the microwave 2D solver
// (apps/microwave/solver.js: solve2D). Three studies, printed as markdown:
//
//   1. Radial-parabola benchmark: public knobs alone can reduce solve2D to
//      pure conduction with a uniform bed source (manual constant bed k,
//      radiation off, gas exchange off, near-uniform field and penetration
//      depth), so the mid-plane bed profile must match T(r)-T(0) =
//      -q r^2/(4 k_bed); the residual is finite-bed axial leakage, which
//      must shrink as the bed is made longer.
//   2. Physical-case grid convergence: the app's default calibrated case on
//      a doubling grid sequence; Richardson extrapolation for the observed
//      order and the default 30×60 grid's distance from the extrapolated
//      answer, plus energy-balance closure per grid.
//   3. Darcy flow: discrete mass conservation of the flow field (the solver
//      reports its own worst per-cell imbalance) per grid.
//
// A full-domain manufactured solution is deliberately out of scope here: the
// tube-gas and outside-air conductivities are temperature-dependent by design
// and are not overridable through public parameters, so a uniform-k domain
// cannot be produced without modifying shipped code. The shared FV pattern
// (2-point flux, harmonic interface resistance) is MMS-verified in the Joule
// solver; this file verifies the microwave implementation against analytic
// conduction and grid-refinement behavior.
//
// Run: node tools/verification/microwave.mjs
"use strict";
import { solve2D, kBed } from "../../apps/microwave/solver.js";
import { richardson, markdownTable, sci, fix } from "./common.mjs";

// The page's default inputs: the reduced-rutile TiO₂ calibration profile with
// every field at its index.html default, with the same derived quantities
// (volume, bulk density, void fractions) parameters() computes.
export function defaultParams(overrides = {}) {
  const clamp = (v, lo, hi) => Math.min(hi, Math.max(lo, v));
  const rhoSolid = 4230, refMass = 1.150e-3, refVolume = 1.18e-6;
  const p = {
    rhoSolid, cpSolid: 690, P: 26, frequency: 2.404e9,
    D: 0.010, H: 0.015, mass: 1.150, gas: "He", flow: 50, pressure: 1e5,
    dp: 50e-6, Ta: 20,
    Nr: 30, Nz: 60, domainWidth: 0.03, domainHeight: 0.03,
    k200: 1.70, k500: 0.62, k800: 0.58, kzRatio: 1.20, hContact: 1000,
    kq: 1.40, tq: 1e-3, airFactor: 3, boundaryMode: "automatic", hBoundary: 12,
    epsTube: 0.85, radArea: 1.80, gasTransferMode: "automatic", gasEff: 0.80,
    dielectricMode: "looyenga", bedKMode: "automatic",
    fieldWr: 1.20, fieldWz: 1.20, fbgR: 0, fbgZ: 0,
    diel: [[20, 3.49, 0.3204], [100, 5.8, 3.0], [200, 6.2, 4.0], [400, 6.5, 4.8], [600, 6.7, 5.4], [725, 7.2, 6.2]],
    maxIter: 6000, tol: 3e-4, omega: 1.05,
    ...overrides,
  };
  if (!Object.prototype.hasOwnProperty.call(p, "volume") || overrides.volume === undefined) p.volume = Math.PI * p.D * p.D * p.H / 4 * 1e6;
  const rhoBulk = (p.mass / 1000) / (p.volume * 1e-6);
  p.rhoBulk = rhoBulk;
  p.voidFraction = clamp(1 - rhoBulk / rhoSolid, 0.01, 0.99);
  p.referenceVoidFraction = clamp(1 - (refMass / refVolume) / rhoSolid, 0.01, 0.99);
  return p;
}

// ---------------------------------------------------------------- study 1
// Pure-conduction reduction. Everything that is not conduction is switched
// off through public parameters; the microwave source is flattened by huge
// field widths and a huge penetration depth (tiny loss factor).
export function radialParabola(HoverD, { Nr = 30, Nz = 60 } = {}) {
  const kConst = 0.5;
  const D = 0.011, H = HoverD * D;
  const p = defaultParams({
    D, H, Nr, Nz, domainWidth: 0.03, domainHeight: H + 0.011,
    k200: kConst, k500: kConst, k800: kConst, bedKMode: "manual", kzRatio: 1,
    gasTransferMode: "manual", gasEff: 1e-9,        // no gas heat exchange
    epsTube: 0,                                     // no wall radiation
    boundaryMode: "manual", hBoundary: 10,
    fieldWr: 1e4, fieldWz: 1e4,                     // flat field shape
    dielectricMode: "manual",
    diel: [[20, 5, 1e-6], [2000, 5, 1e-6]],         // huge penetration depth
    P: 3, tol: 1e-6, maxIter: 60000, omega: 1.5, flowMode: "off",
  });
  const sol = solve2D(p);
  if (!sol.converged) throw new Error(`solve2D did not converge (maxDelta ${sol.maxDelta})`);
  const { T, R, dr, dz } = sol;
  const qVol = p.P / sol.V;
  const jmid = Math.floor(p.Nz / 2);
  let worst = 0, bedCells = 0;
  for (let i = 1; (i + 0.5) * dr < R; i++) {
    const r = (i + 0.5) * dr, r0 = 0.5 * dr;
    const exact = -qVol * (r * r - r0 * r0) / (4 * kConst);
    const got = T[jmid][i] - T[jmid][0];
    worst = Math.max(worst, Math.abs(got - exact) / Math.abs(exact));
    bedCells++;
  }
  const dTexact = qVol * R * R / (4 * kConst);
  return { worstRelative: worst, centerToSurfaceK: dTexact, qVol, bedCells };
}

// ---------------------------------------------------------------- study 2
// Default calibrated case on a doubling grid sequence. domainWidth/Height are
// physical lengths, so the domain is identical on every grid by construction.
export function physicalConvergence() {
  const rows = [];
  for (const level of [0, 1, 2]) {
    const p = defaultParams({
      Nr: 30 << level, Nz: 60 << level,
      tol: 3e-4 / (1 << level), maxIter: 60000, flowMode: "on",
    });
    const sol = solve2D(p);
    if (!sol.converged) throw new Error(`grid ${p.Nr}x${p.Nz} did not converge`);
    // Refinement has a floor. Once the mesh cell drops below the packing unit
    // cell the grid is resolving structure the continuum model does not carry,
    // and an extrapolation through such a level converges an equation that has
    // stopped describing the bed. This profile's 50 um powder leaves room at
    // 120x240 (h/dp = 2.5); the coarser 194 um SiC would not, and a study on it
    // has to stop at 60x120.
    if (sol.homogenization.resolvedBelowUnitCell) {
      throw new Error(`grid ${p.Nr}x${p.Nz} refines below the ${(p.dp * 1e6).toFixed(0)} um unit cell `
        + `(h/dp = ${sol.homogenization.cellPerParticle.toFixed(2)}); the continuum model does not reach there`);
    }
    rows.push({
      grid: `${p.Nr}×${p.Nz}`,
      center: sol.center, wall: sol.wall, Tavg: sol.Tavg,
      balance: Math.abs(sol.balance) / p.P,
      massImbalance: sol.darcy ? sol.darcy.maxMassImbalance / sol.darcy.massFlow : null,
      it: sol.it, cellPerParticle: sol.homogenization.cellPerParticle,
      linearResidual: sol.linearResidual,
    });
  }
  const rich = {
    center: richardson(rows[0].center, rows[1].center, rows[2].center),
    Tavg: richardson(rows[0].Tavg, rows[1].Tavg, rows[2].Tavg),
  };
  return { rows, rich };
}

// ---------------------------------------------------------------- report
function main() {
  console.log("## Microwave 2D solver: numerical verification\n");

  console.log("### 1. Radial parabola (pure-conduction reduction)\n");
  const parabolaRows = [];
  for (const hd of [4, 8, 16]) {
    const p = radialParabola(hd);
    parabolaRows.push([`H/D = ${hd}`, fix(p.centerToSurfaceK, 2) + " K", sci(p.worstRelative, 2)]);
  }
  console.log(markdownTable(["case", "analytic center→surface rise", "worst relative mismatch"], parabolaRows) + "\n");

  console.log("### 2. Default calibrated case: grid sensitivity\n");
  const { rows, rich } = physicalConvergence();
  console.log(markdownTable(
    ["grid", "center T (°C)", "wall T (°C)", "avg bed T (°C)", "energy closure", "Darcy mass imbalance", "sweeps"],
    rows.map((r) => [r.grid, fix(r.center, 2), fix(r.wall, 2), fix(r.Tavg, 2), sci(r.balance, 2), r.massImbalance === null ? "—" : sci(r.massImbalance, 2), r.it]),
  ) + "\n");
  const last = rows[rows.length - 1];
  const asymptotic = rich.Tavg.p > 0.5 && rich.Tavg.p < 3;
  if (asymptotic) {
    console.log(`Richardson (avg bed T): order ${fix(rich.Tavg.p, 2)}, extrapolated ${fix(rich.Tavg.qExtrap, 2)} °C, finest-grid error ${sci(rich.Tavg.fineError, 2)}`);
  } else {
    console.log(`Richardson order not meaningful here (grid-to-grid differences are not yet asymptotic: the near-wall exponential source deposition and T^4 wall radiation are resolved progressively with the grid). Reporting sensitivity against the finest grid instead:`);
  }
  console.log(`Default 30×60 grid vs finest ${last.grid}: center T differs by ${fix(Math.abs(rows[0].center - last.center), 2)} K (${fix(100 * Math.abs(rows[0].center - last.center) / (last.center - 20), 2)}% of the rise), avg bed T by ${fix(Math.abs(rows[0].Tavg - last.Tavg), 2)} K, wall T by ${fix(Math.abs(rows[0].wall - last.wall), 2)} K.`);
}

if (import.meta.url === `file://${process.argv[1]}`) main();
