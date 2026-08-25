// Refit the calibrated bed parameters against the solved field.
//
// The page fits four numbers -- the bed conductivity anchors k200, k500, k800
// and the radiation-area multiplier radArea -- by minimising the combined RMSE
// against a measured power sweep of [P_abs, T_wall, T_centre]. Those fits were
// made with the heating source shaped by a fitted Gaussian whose axis-to-edge
// power ratio is 1.809 where the solved field gives 1.264, so whatever error a
// 43%-too-peaked source produced had to be absorbed somewhere, and the bed
// conductivity is the parameter positioned to absorb it.
//
// This runs the same bounded coordinate search the page runs, from the same
// starting point and against the same data, once per source model, and prints
// both fits side by side. The comparison is the point: if the fitted
// conductivities move, they were carrying the source error; if the RMSE
// improves, the solved field describes the measurements better.
//
// The search is the page's own -- one parameter at a time, both directions,
// step fractions of each bound -- extended with two finer passes so that neither
// model is judged on a stopping point rather than an optimum. Both models get
// the identical schedule.
//
// Run: node tools/verification/microwave-calibrate.mjs [--material sic|rutile] [--quick]
"use strict";
import { materialProfiles, parseRows, solve2D, clamp } from "../../apps/microwave/solver.js";
import { markdownTable, fix } from "./common.mjs";

const FRACTIONS = [0.20, 0.08, 0.03, 0.012];

export function profileFor(name) {
  const entries = Object.entries(materialProfiles);
  const hit = entries.find(([key, profile]) => key === name
    || new RegExp(name, "i").test(profile.label) || new RegExp(name, "i").test(profile.formula || ""));
  if (!hit) throw new Error(`no material profile matching "${name}" (have: ${entries.map(([k]) => k).join(", ")})`);
  return { key: hit[0], profile: hit[1] };
}

export function baseParameters(profile, overrides = {}) {
  const d = profile.defaults, D = d.diameter / 1000, H = d.length / 1000;
  const volume = Math.PI * D * D * H / 4 * 1e6, rhoBulk = (d.mass / 1000) / (volume * 1e-6);
  const p = {
    rhoSolid: profile.rhoSolid, P: d.pabs, frequency: d.frequency * 1e9,
    D, H, volume, mass: d.mass, tq: d["tube-thickness"] / 1000,
    gas: d.gas, flow: d.flow, pressure: d["gas-pressure"] * 1e5,
    dp: d["particle-diameter"] * 1e-6, Ta: d.ambient,
    Nr: 30, Nz: 60, domainWidth: 0.03, domainHeight: 0.03,
    k200: d.k200, k500: d.k500, k800: d.k800, kzRatio: d["kz-ratio"],
    hContact: d["h-contact"], kq: d["k-quartz"], airFactor: d["air-factor"],
    boundaryMode: d["boundary-mode"], hBoundary: d["h-boundary"],
    epsTube: d.emissivity, radArea: d["rad-area"],
    gasTransferMode: d["gas-transfer-mode"], gasEff: d["gas-eff"],
    dielectricMode: d["dielectric-mode"], bedKMode: d["bed-k-mode"],
    fieldWr: d["field-wr"], fieldWz: d["field-wz"], fieldMode: profile.fieldMode,
    fbgR: d["fbg-r"] / 1000, fbgZ: d["fbg-z"] / 1000,
    diel: parseRows(profile.dielectric, 3),
    rhoBulk, voidFraction: clamp(1 - rhoBulk / profile.rhoSolid, 0.05, 0.95),
    maxIter: 6000, tol: 3e-4, omega: 1.05, flowMode: "off",
    ...overrides,
  };
  p.referenceVoidFraction = p.voidFraction;
  return p;
}

const rootMeanSquare = (values) => Math.sqrt(values.reduce((sum, v) => sum + v * v, 0) / Math.max(1, values.length));

// The page searched on a 10x30 mesh. With the linear system relaxed by
// Gauss-Seidel that was the only affordable choice, and it put the objective 24 K
// from the converged centre temperature while the residual being minimised was
// about 10 K -- so the search moved parameters to cancel discretisation error.
// The Krylov solve makes 30x60 cost about a second, where the same error is
// 3.5 K and sits under the residual, so the search mesh is now the report mesh.
export function evaluate(base, rows, coarse = true) {
  const results = [];
  for (const row of rows) {
    const p = { ...base, P: Math.max(0, row[0]),
      Nr: coarse ? (base.searchNr ?? 30) : base.Nr, Nz: coarse ? (base.searchNz ?? 60) : base.Nz,
      maxIter: base.maxIter, tol: base.tol, omega: base.omega };
    const s = solve2D(p);
    results.push({ P: row[0], wallObs: row[1], centreObs: row[2], centre: s.fbg, wall: s.wall, converged: s.converged });
  }
  const powered = results.filter((r) => r.P > 0);
  const centreErrors = powered.map((r) => r.centre - r.centreObs), wallErrors = powered.map((r) => r.wall - r.wallObs);
  return { results, centreRMSE: rootMeanSquare(centreErrors), wallRMSE: rootMeanSquare(wallErrors),
    totalRMSE: rootMeanSquare(centreErrors.concat(wallErrors)) };
}

export function calibrate(profile, rows, { fieldMode, fractions = FRACTIONS } = {}) {
  const bounds = profile.fitBounds;
  const base = baseParameters(profile, fieldMode ? { fieldMode } : {});
  let score = evaluate(base, rows, true).totalRMSE;
  const start = score, path = [];
  for (const fraction of fractions) {
    for (const key of Object.keys(bounds)) {
      let chosen = base[key], chosenScore = score;
      for (const direction of [-1, 1]) {
        const trial = { ...base }, span = bounds[key][1] - bounds[key][0];
        trial[key] = clamp(base[key] + direction * fraction * span, bounds[key][0], bounds[key][1]);
        const trialScore = evaluate(trial, rows, true).totalRMSE;
        if (trialScore < chosenScore) { chosen = trial[key]; chosenScore = trialScore; }
      }
      base[key] = chosen; score = chosenScore;
    }
    path.push({ fraction, score });
  }
  return { fitted: base, coarseRMSE: score, startRMSE: start, path, final: evaluate(base, rows, false) };
}

function main() {
  const argv = process.argv.slice(2);
  const name = argv.includes("--material") ? argv[argv.indexOf("--material") + 1] : "SiC";
  const fractions = argv.includes("--quick") ? [0.20, 0.08] : FRACTIONS;
  const { profile } = profileFor(name);
  const rows = parseRows(profile.experiments, 3);

  console.log(`## Recalibrating ${profile.label}\n`);
  console.log(`dataset: ${profile.calibrationLabel} · ${rows.length} rows, ${rows.filter((r) => r[0] > 0).length} powered\n`);

  const runs = [
    { label: "fitted Gaussian source (shipped)", fieldMode: undefined },
    { label: "solved Helmholtz field", fieldMode: "helmholtz" },
  ].map((run) => ({ ...run, ...calibrate(profile, rows, { fieldMode: run.fieldMode, fractions }) }));

  console.log(markdownTable(
    ["source model", "k200", "k500", "k800", "radArea", "coarse RMSE", "final centre RMSE", "final wall RMSE", "final combined"],
    runs.map((r) => [r.label,
      fix(r.fitted.k200, 4), fix(r.fitted.k500, 4), fix(r.fitted.k800, 4), fix(r.fitted.radArea, 4),
      fix(r.coarseRMSE, 2), fix(r.final.centreRMSE, 2), fix(r.final.wallRMSE, 2), fix(r.final.totalRMSE, 2)]),
  ) + "\n");

  const shipped = baseParameters(profile);
  console.log(markdownTable(
    ["parameter", "shipped default", "refit on fitted source", "refit on solved field", "solved vs shipped"],
    Object.keys(profile.fitBounds).map((k) => [k,
      fix(shipped[k], 4), fix(runs[0].fitted[k], 4), fix(runs[1].fitted[k], 4),
      `${fix(100 * (runs[1].fitted[k] / Math.max(shipped[k], 1e-30) - 1), 1)}%`]),
  ) + "\n");

  console.log(`Combined RMSE from the shipped defaults before any search: ${fix(runs[0].startRMSE, 2)} °C on the fitted`
    + ` source, ${fix(runs[1].startRMSE, 2)} °C on the solved field.`);
  console.log(`Both models ran the identical schedule [${fractions.join(", ")}].`);
}

if (import.meta.url === `file://${process.argv[1]}`) main();
