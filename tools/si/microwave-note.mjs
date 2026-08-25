// Reproduces every table in docs/si/microwave-thermal-note.md from the shipped
// microwave solver. Nothing here is fitted at run time: the thermal parameters
// are the page's calibrated defaults, and the only per-sample inputs are the
// measured tap density and the measured dielectric function.
//
// With no flag the default tables are printed. Each flag selects one study
// instead; pass --all to add the default tables back alongside.
//
// Timings below are measured on a 4-core container. Everything except --grid
// runs on the default 30x60 grid at roughly six solves a minute, so cost tracks
// the solve count; --grid is the outlier because four of its twelve solves are
// at 120x240.
//
// Run: node tools/si/microwave-note.mjs             (Tables S1, S3, S5, S6; ~5 min)
//      node tools/si/microwave-note.mjs --band      (Table S4;              ~7 min)
//      node tools/si/microwave-note.mjs --invert    (Table S7;             ~20 min)
//      node tools/si/microwave-note.mjs --grid      (thin-skin check;    ~1 h 45)
"use strict";
import { solve2D } from "../../apps/microwave/solver.js";
import { defaultParams } from "../verification/microwave.mjs";
import { markdownTable, fix } from "../verification/common.mjs";

const clamp = (v, lo, hi) => Math.min(hi, Math.max(lo, v));

// The calibration series: P_sample (W) -> [T_wall, T_center] in °C.
const SERIES = {
  5: [167.3, 185.5], 10: [287.1, 354.6], 14: [360, 490], 17: [402.7, 580],
  20: [425.5, 618], 23: [468, 710], 26: [499.6, 799.2],
};

// The reduced-rutile ε(T) table the calibration is built on, and its
// room-temperature anchor. Per-sample tables are this shape rescaled to the
// sample's own measured room-temperature ε′, ε″.
const BASE_DIEL = defaultParams().diel;
const [REF_EP, REF_EPP] = [3.49, 0.3204];

// name, tap density (g/mL), solid density (g/mL), ε′(RT), ε″(RT),
// reported P_sample (W), measured T_wall (°C)
const SAMPLES = [
  ["a-TiO₂-R600",  0.405, 3.89,  3.26,  0.0040,  0.79,  94.0],
  ["r-TiO₂-R600",  0.975, 4.23,  3.49,  0.3204, 14.88, 402.7],
  ["r-TiO₂-R1000", 0.757, 4.23,  5.21,  1.2813, 17.48, 436.0],
  ["r-TiO₂-R1100", 1.185, 4.23, 12.93,  6.5240, 17.55, 386.0],
  ["Ti₂O₃",        2.340, 4.49, 36.38, 52.7100,  7.99, 221.0],
];

// The constant residual shared by the three low-loss samples (Table S6); read
// as a calibration difference between the power-sweep cavity and the dual-mode
// cavity, and removed before inverting.
const CAVITY_OFFSET = 35;
const GAMMA = 3e-3; // cavity filling factor, Supplementary Note 2b

// A sample at its own packing: the ε table is already the bed's, so reference
// and actual void fraction coincide and Looyenga mixing is a no-op.
function sampleParams(tapDensity, rhoSolid, diel, overrides = {}) {
  const p = defaultParams({
    rhoSolid: rhoSolid * 1000, mass: tapDensity * 1.18,
    dielectricMode: "manual", diel, ...overrides,
  });
  p.rhoBulk = (p.mass / 1000) / (p.volume * 1e-6);
  p.voidFraction = clamp(1 - p.rhoBulk / p.rhoSolid, 0.01, 0.99);
  p.referenceVoidFraction = p.voidFraction;
  return p;
}

const rescale = (ep, epp) => BASE_DIEL.map(([T, e1, e2]) => [T, e1 * ep / REF_EP, e2 * epp / REF_EPP]);

// Every study here costs tens of minutes, so a flag selects one rather than
// appending it to a full run: asking for the grid check should not first
// recompute four tables you already have.
const flag = (f) => process.argv.includes(f);
const selective = ["--band", "--grid", "--invert"].some(flag);
const cheapTables = !selective || flag("--all");

// ------------------------------------------------------------- Table S1
if (cheapTables) {
console.log("### Table S_x.1 — calibration against the two-thermometer sweep\n");
let sw = 0, sc = 0, n = 0;
const t1 = [];
for (const P of Object.keys(SERIES).map(Number)) {
  const s = solve2D(defaultParams({ P }));
  const [ew, ec] = SERIES[P];
  sw += (s.wall - ew) ** 2; sc += (s.center - ec) ** 2; n++;
  t1.push([P, fix(s.wall, 1), fix(ew, 1), fix(s.wall - ew, 1),
           fix(s.center, 1), fix(ec, 1), fix(s.center - ec, 1),
           fix(s.Tavg, 1), fix(s.center - s.Tavg, 1)]);
}
console.log(markdownTable(
  ["P_sample (W)", "T_wall model", "T_wall meas.", "Δ", "T_center model", "T_center meas.", "Δ", "⟨T⟩_V", "T_center − ⟨T⟩_V"], t1));
console.log(`\nRMSE: T_wall ${fix(Math.sqrt(sw / n), 1)} K, T_center ${fix(Math.sqrt(sc / n), 1)} K\n`);
}

// ------------------------------------------------------------- Table S3
if (cheapTables) {
console.log("### Table S_x.3 — power partitioning\n");
const t3 = [];
for (const P of [5, 10, 17, 26]) {
  const s = solve2D(defaultParams({ P }));
  const pct = (x) => fix(100 * x / P, 1) + "%";
  t3.push([P, pct(s.qrad), pct(s.qBoundary), pct(s.qgas), pct(s.qRadialBed), pct(s.qAxialBed)]);
}
console.log(markdownTable(
  ["P_sample (W)", "quartz radiation", "outer convection", "gas enthalpy", "→ radial (internal)", "→ axial (internal)"], t3));

// Sweep-flow sensitivity quoted in the same section: the helium is not a
// meaningful heat-removal path, so flow rate is not a temperature handle.
console.log("\nSweep-flow sensitivity at 26 W:");
for (const flow of [25, 50, 100]) {
  const s = solve2D(defaultParams({ P: 26, flow }));
  console.log(`  ${String(flow).padStart(3)} sccm  T_wall ${fix(s.wall, 1)} °C  T_center ${fix(s.center, 1)} °C  gas share ${fix(100 * s.qgas / 26, 1)}%`);
}
}

// ------------------------------------------------------------- Table S4
// Field-width sensitivity. Only w_r is varied; the remaining calibrated
// parameters stay at their jointly fitted values, so this measures how much
// each thermometer *responds* to the deposition profile, not a re-fit.
if (flag("--band")) {
  console.log("\n### Table S_x.4 — sensitivity to the fitted field width\n");
  const t4b = [];
  for (const wr of [0.8, 1.0, 1.2, 1.5, 2.0, 3.0, 6.0, 20.0]) {
    let bw = 0, bc = 0, bn = 0;
    for (const P of Object.keys(SERIES).map(Number)) {
      const s = solve2D(defaultParams({ P, fieldWr: wr }));
      bw += (s.wall - SERIES[P][0]) ** 2; bc += (s.center - SERIES[P][1]) ** 2; bn++;
    }
    const s26 = solve2D(defaultParams({ P: 26, fieldWr: wr }));
    t4b.push([fix(wr, 1), fix(Math.sqrt(bw / bn), 1), fix(Math.sqrt(bc / bn), 1),
              fix(s26.center, 0), fix(s26.Tavg, 0), fix(s26.center - s26.Tavg, 0)]);
  }
  console.log(markdownTable(
    ["w_r", "RMSE T_wall (K)", "RMSE T_center (K)", "T_center (°C)", "⟨T⟩_V (°C)", "gap (K)"], t4b));
}

// ------------------------------------------------------------- Table S5
if (cheapTables) {
console.log("\n### Table S_x.5 — response to the deposition profile at fixed power\n");
const t4 = [];
for (const x of [1, 2, 4, 8, 16, 32, 64]) {
  const p = defaultParams({ P: 26, dielectricMode: "manual", diel: BASE_DIEL.map(([T, e1, e2]) => [T, e1, e2 * x]) });
  p.referenceVoidFraction = p.voidFraction;
  const s = solve2D(p);
  // hottest bed cell anywhere in the bed; radius from the cell size, not the
  // full domain width (the domain is wider than the tube).
  let best = -Infinity, bi = 0;
  for (let j = 0; j < s.T.length; j++) for (let i = 0; i < s.T[j].length; i++) {
    if (s.material[j][i] !== 2) continue;
    if (s.T[j][i] > best) { best = s.T[j][i]; bi = i; }
  }
  // fraction of the deposited power landing in the outer 20 % of the radius
  let pOut = 0, pTot = 0;
  for (let j = 0; j < s.T.length; j++) for (let i = 0; i < s.T[j].length; i++) {
    if (s.material[j][i] !== 2) continue;
    const vol = Math.PI * ((((i + 1) * s.dr) ** 2) - ((i * s.dr) ** 2)) * s.dz;
    const q = s.heat[j][i] * vol;
    pTot += q; if ((i + 0.5) * s.dr > 0.8 * s.R) pOut += q;
  }
  t4.push([`×${x}`, fix(s.dpCenter * 1000, 2), fix(s.dpCenter / s.R, 2), fix(s.wall, 0),
           fix(s.center, 0), fix(best, 0), fix((bi + 0.5) * s.dr * 1000, 1),
           fix(100 * pOut / pTot, 0) + "%"]);
}
console.log(markdownTable(
  ["ε″ scaling", "δ_p (mm)", "δ_p/R", "T_wall (°C)", "T_center (°C)", "T_max (°C)", "r(T_max) (mm)", "power in outer 20%"], t4));
}

// ------------------------------------------------------------- Table S6
if (cheapTables) {
console.log("\n### Table S_x.6 — per-sample reconstruction\n");
const t5 = [];
for (const [name, tap, rhoS, ep, epp, P, Tmeas] of SAMPLES) {
  const p = sampleParams(tap, rhoS, rescale(ep, epp), { P });
  const s = solve2D(p);
  t5.push([name, fix(tap, 3), fix(p.voidFraction, 3), fix(P, 2), fix(s.dpCenter / s.R, 2),
           fix(s.wall, 0), fix(Tmeas, 0), fix(s.wall - Tmeas, 0)]);
}
console.log(markdownTable(
  ["sample", "tap ρ (g/mL)", "void", "P_sample (W)", "δ_p/R", "T_wall model", "T_wall meas.", "Δ (K)"], t5));
}

// -------------------------------------------------- thin-skin grid check
// The two most reduced samples have a penetration depth of order one radial
// cell on the default grid. Section S_x.7 rests on T_wall being insensitive to
// that, so check it directly on the doubling sequence of Table S_x.2.
if (flag("--grid")) {
  console.log("\n### Thin-skin grid check — T_wall at the reported power\n");
  const tg = [];
  for (const [name, tap, rhoS, ep, epp, P] of SAMPLES.slice(1)) {
    const diel = rescale(ep, epp);
    const row = [name];
    let dpR = 0;
    for (const lv of [0, 1, 2]) {
      const s = solve2D(sampleParams(tap, rhoS, diel, {
        P, Nr: 30 << lv, Nz: 60 << lv, tol: 3e-4 / (1 << lv), maxIter: 60000,
      }));
      row.push(fix(s.wall, 1));
      dpR = s.dpCenter / s.R;
    }
    tg.push([row[0], fix(dpR, 2), row[1], row[2], row[3],
             fix(Math.abs(Number(row[1]) - Number(row[3])), 1)]);
  }
  console.log(markdownTable(
    ["sample", "δ_p/R", "10 cells/R", "20 cells/R", "40 cells/R", "drift (K)"], tg));
}

// ------------------------------------------------------------- Table S7
if (flag("--invert")) {
  console.log("\n### Table S_x.7 — power implied by the measured wall temperature\n");
  const t6 = [];
  for (const [name, tap, rhoS, ep, epp, P, Tmeas] of SAMPLES) {
    const diel = rescale(ep, epp);
    let lo = 0.05, hi = 25, mid = 0;
    for (let k = 0; k < 26; k++) {
      mid = (lo + hi) / 2;
      const s = solve2D(sampleParams(tap, rhoS, diel, { P: mid }));
      if (s.wall + CAVITY_OFFSET < Tmeas) lo = mid; else hi = mid;
    }
    const s0 = solve2D(sampleParams(tap, rhoS, diel, { P }));
    t6.push([name, fix(GAMMA * epp, 4), fix(s0.dpCenter / s0.R, 2),
             fix(P, 2), fix(mid, 2), fix(mid / P, 2)]);
  }
  console.log(markdownTable(
    ["sample", "γε″", "δ_p/R", "P_sample from CPT (W)", "P_sample implied by T_wall (W)", "ratio"], t6));
}
