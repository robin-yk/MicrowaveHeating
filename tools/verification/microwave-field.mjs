// Control and consequence for the microwave frequency-domain field solve.
//
// The analytic anchor is the classical one: an infinite lossy dielectric
// cylinder of radius R sitting in a uniform axial E field carries
//
//     E_z(r) / E_z(R)  =  J0(k r) / J0(k R),     k = k0 sqrt(eps)
//
// with J0 of a complex argument. Refraction concentrates the field on the axis
// whenever |kR| is an appreciable fraction of the first zero at 2.405, and that
// is the behaviour the shipped Beer-Lambert skin cannot produce at all: it can
// only attenuate inward from the surface. Comparing the two against the exact
// answer is the point of this file.
//
// The shape is normalised at r = R before comparison, so the test does not
// depend on how much the finite domain perturbs the field outside the load.
//
// Run: node tools/verification/microwave-field.mjs
"use strict";
import { bedMesh, solveField2D, penetrationDepth, dielectric, c0 } from "../../apps/microwave/solver.js";
import { markdownTable, fix, sci } from "./common.mjs";

// J0(z) for complex z by its Maclaurin series; |z| stays below 1 here, where the
// series is exact to machine precision within a few terms.
function besselJ0(zr, zi) {
  let sumR = 1, sumI = 0, termR = 1, termI = 0;
  for (let m = 1; m < 60; m++) {
    const z2r = zr * zr - zi * zi, z2i = 2 * zr * zi, scale = -1 / (4 * m * m);
    const nextR = (termR * z2r - termI * z2i) * scale, nextI = (termR * z2i + termI * z2r) * scale;
    termR = nextR; termI = nextI; sumR += termR; sumI += termI;
    if (Math.hypot(termR, termI) < 1e-18) break;
  }
  return { re: sumR, im: sumI };
}

const EPS_REAL = 7.959225, EPS_LOSS = 0.398952;   // SiC at 20 C, 2.404 GHz

function input({ D = 0.010, H = 0.19, Nr = 30, Nz = 60, domainWidth = 0.03, domainHeight = 0.20 } = {}) {
  const flat = [[0, EPS_REAL, EPS_LOSS], [2000, EPS_REAL, EPS_LOSS]];
  return {
    D, H, tq: 0.001, Nr, Nz, domainWidth, domainHeight,
    frequency: 2.404e9, diel: flat, dielectricMode: "manual",
    voidFraction: 0.4, referenceVoidFraction: 0.4,
    fieldWr: 1.2, fieldWz: 1.2,
  };
}

export function control(overrides = {}) {
  const p = input(overrides), mesh = bedMesh(p);
  const T = Array.from({ length: p.Nz }, () => Array(p.Nr).fill(20));
  const field = solveField2D({ p, T, mesh });
  const k0 = 2 * Math.PI * p.frequency / c0;
  const magnitude = Math.hypot(EPS_REAL, EPS_LOSS), angle = Math.atan2(-EPS_LOSS, EPS_REAL) / 2;
  const root = Math.sqrt(magnitude), nRe = root * Math.cos(angle), nIm = root * Math.sin(angle);
  const kRe = k0 * nRe, kIm = k0 * nIm, R = p.D / 2;
  const jMid = Math.floor(p.Nz / 2), rows = [];
  // Last bed cell on the mid-plane row is the normalisation point.
  let iEdge = 0;
  for (let i = 0; i < p.Nr; i++) if (mesh.material[jMid][i] === 2) iEdge = i;
  const at = (i) => Math.sqrt(field.magnitude[jMid * p.Nr + i]);
  const edge = at(iEdge), rEdge = (iEdge + 0.5) * mesh.dr;
  const jEdge = besselJ0(kRe * rEdge, kIm * rEdge), edgeMag = Math.hypot(jEdge.re, jEdge.im);
  for (let i = 0; i <= iEdge; i++) {
    const r = (i + 0.5) * mesh.dr, j = besselJ0(kRe * r, kIm * r);
    const exact = Math.hypot(j.re, j.im) / edgeMag, solved = at(i) / edge;
    // What the shipped model claims at the same point, on the same scale. Both
    // factors, not just the skin: heat = eps'' * gaussian * attenuation, so the
    // product is the power weighting and must be compared against the exact
    // solution squared. Quoting the skin alone -- as an earlier version of this
    // file did -- makes the shipped source look edge-peaked when the fitted
    // Gaussian in front of it dominates and peaks it hard on the axis.
    const delta = penetrationDepth(20, p);
    const shipped = (rr) => Math.exp(-Math.pow(rr / (Math.max(0.05, p.fieldWr) * R), 2)) * Math.exp(-(R - rr) / delta);
    const fitted = shipped(r) / shipped(rEdge);
    const skin = Math.exp(-(R - r) / delta) / Math.exp(-(R - rEdge) / delta);
    rows.push({ r, solved, exact, fitted, skin, exactPower: exact * exact, error: solved / exact - 1 });
  }
  return { rows, kR: Math.hypot(kRe * R, kIm * R), field, p, mesh,
    lambdaMaterial: (c0 / p.frequency) / nRe, delta: penetrationDepth(20, p) };
}

function main() {
  console.log("## Microwave field solve against the exact lossy-cylinder solution\n");
  const base = control();
  console.log(`|kR| = ${fix(base.kR, 4)}   lambda in material ${fix(base.lambdaMaterial * 1000, 1)} mm`
    + `   penetration depth ${fix(base.delta * 1000, 0)} mm   bed radius ${fix(base.p.D / 2 * 1000, 1)} mm`);
  console.log(`field solve: ${base.field.iterations} iterations, residual ${sci(base.field.relativeResidual, 2)}\n`);

  console.log(markdownTable(
    ["r (mm)", "solved |E|/|E(R)|", "exact J0 ratio", "solved − exact",
      "exact POWER", "shipped fitted POWER", "skin factor alone"],
    base.rows.filter((_, i) => i % 2 === 0).map((r) => [
      fix(r.r * 1000, 2), fix(r.solved, 5), fix(r.exact, 5), sci(r.error, 2),
      fix(r.exactPower, 4), fix(r.fitted, 4), fix(r.skin, 4),
    ]),
  ) + "\n");

  const worst = base.rows.reduce((a, b) => (Math.abs(b.error) > Math.abs(a.error) ? b : a));
  console.log(`Worst shape error against the exact solution: ${sci(worst.error, 3)} at r = ${fix(worst.r * 1000, 2)} mm.`);
  const axis = base.rows[0];
  console.log(`On the axis, power density relative to the bed edge:`);
  console.log(`  exact (refraction)      ${fix(axis.exactPower, 4)}x`);
  console.log(`  shipped fitted source   ${fix(axis.fitted, 4)}x   -> ${fix(100 * (axis.fitted / axis.exactPower - 1), 0)}% too peaked`);
  console.log(`  its skin factor alone   ${fix(axis.skin, 4)}x   (edge-peaked, and only ${fix(100 * Math.abs(axis.skin - 1), 1)}% either way)`);
  console.log(`\nThe centre peaking is real physics -- refraction into a subwavelength load -- but`);
  console.log(`the shipped model produces it through a fitted Gaussian width, not from the field,`);
  console.log(`and at fieldWr = ${fix(base.p.fieldWr, 2)} that fit overshoots the exact answer by half again.\n`);

  console.log("### Grid convergence of the field shape\n");
  const levels = [{ Nr: 15, Nz: 30 }, { Nr: 30, Nz: 60 }, { Nr: 60, Nz: 120 }];
  const axisValues = levels.map((level) => {
    const run = control(level);
    const first = run.rows[0];
    return [`${level.Nr}x${level.Nz}`, fix(first.solved, 6), fix(first.exact, 6), sci(first.error, 3)];
  });
  console.log(markdownTable(["grid", "solved on axis", "exact", "error"], axisValues) + "\n");
}

if (import.meta.url === `file://${process.argv[1]}`) main();
