// Regression tests for the microwave-heating numeric core
// (microwave-solver.js). Pure Node, no browser: run with `node --test`.
import test from "node:test";
import assert from "node:assert/strict";
import {
  materialProfiles, clamp, dielectric, kBed, penetrationDepth, solve2D, transportNumbers,
  axisymmetricTensor, homogenizationValidity, porousContinuumClosures, darcyVelocity
} from "../microwave-solver.js";

// Mirrors microwave.html's parameters() for the two embedded material
// profiles, using each profile's own defaults instead of guessing values,
// so solves land in the same regime the page itself calibrates against.
function makeParams(profileId, overrides = {}) {
  const profile = materialProfiles[profileId];
  assert.ok(profile, `unknown test profile: ${profileId}`);
  const d = profile.defaults;
  const D = d.diameter / 1000, H = d.length / 1000;
  const diel = d.diel ?? profile.dielectric.trim().split(/\n+/).map((line) => line.split(",").map(Number));
  const p = {
    materialId: profileId, materialLabel: profile.label, materialFormula: profile.formula,
    rhoSolid: profile.rhoSolid, cpSolid: profile.cpSolid,
    P: d.pabs, frequency: d.frequency * 1e9, volume: Math.PI * D * D * H / 4 * 1e6, D, H,
    mass: d.mass, gas: d.gas, flow: d.flow, pressure: d["gas-pressure"] * 1e5, dp: d["particle-diameter"] * 1e-6,
    Ta: d.ambient, Nr: 30, Nz: 60, domainWidth: 0.03, domainHeight: 0.03,
    k200: d.k200, k500: d.k500, k800: d.k800, kzRatio: d["kz-ratio"], hContact: d["h-contact"],
    kq: d["k-quartz"], tq: d["tube-thickness"] / 1000, airFactor: d["air-factor"],
    boundaryMode: d["boundary-mode"], hBoundary: d["h-boundary"], epsTube: d.emissivity, radArea: d["rad-area"],
    gasTransferMode: d["gas-transfer-mode"], gasEff: d["gas-eff"], dielectricMode: d["dielectric-mode"],
    bedKMode: d["bed-k-mode"], fieldWr: d["field-wr"], fieldWz: d["field-wz"],
    fbgR: d["fbg-r"] / 1000, fbgZ: d["fbg-z"] / 1000, diel,
    maxIter: 6000, tol: 3e-4, omega: 1.05,
    ...overrides
  };
  const rhoBulk = (p.mass / 1000) / (p.volume * 1e-6);
  const refMass = d.mass / 1000, refVolume = d.volume * 1e-6;
  p.rhoBulk = rhoBulk;
  p.voidFraction = clamp(1 - rhoBulk / p.rhoSolid, 0.01, 0.99);
  p.referenceVoidFraction = clamp(1 - (refMass / refVolume) / p.rhoSolid, 0.01, 0.99);
  return p;
}

test("dielectric() Looyenga mixing reduces to the reference table when void fraction matches the reference", () => {
  const p = makeParams("rutile-reduced-600c-30m");
  const e = dielectric(200, p);
  // reference void fraction equals the actual void fraction for the default mass/volume,
  // so the mixed dielectric should reproduce the table entry at 200 C closely.
  assert.ok(Math.abs(e.ep - 6.2) < 0.05, `ep=${e.ep}`);
  assert.ok(Math.abs(e.epp - 4.0) < 0.05, `epp=${e.epp}`);
});

test("kBed() and penetrationDepth() stay finite and positive across the calibrated range", () => {
  for (const profileId of Object.keys(materialProfiles)) {
    const p = makeParams(profileId);
    for (const T of [20, 200, 500, 800]) {
      const k = kBed(T, p);
      const dp = penetrationDepth(T, p);
      assert.ok(Number.isFinite(k) && k > 0, `${profileId} kBed(${T}) = ${k}`);
      assert.ok(Number.isFinite(dp) && dp > 0, `${profileId} penetrationDepth(${T}) = ${dp}`);
    }
  }
});

test("solve2D() converges for the reduced-rutile TiO2 default profile", () => {
  const p = makeParams("rutile-reduced-600c-30m");
  const sol = solve2D(p);
  assert.ok(sol.converged, `did not converge in ${sol.it + 1} iterations, maxDelta=${sol.maxDelta}`);
  assert.ok(Number.isFinite(sol.center) && sol.center > p.Ta);
  assert.ok(Number.isFinite(sol.wall) && sol.wall > p.Ta);
  const relBalance = Math.abs(sol.balance) / Math.max(sol.p.P, 1e-30);
  assert.ok(relBalance < 0.01, `energy balance residual too large: ${relBalance}`);
});

test("solve2D() converges for the SiC default profile", () => {
  const p = makeParams("sic-60-100-mesh");
  const sol = solve2D(p);
  assert.ok(sol.converged, `did not converge in ${sol.it + 1} iterations, maxDelta=${sol.maxDelta}`);
  const relBalance = Math.abs(sol.balance) / Math.max(sol.p.P, 1e-30);
  assert.ok(relBalance < 0.01, `energy balance residual too large: ${relBalance}`);
});

test("solve2D() is deterministic for identical inputs", () => {
  const p = makeParams("sic-60-100-mesh");
  const a = solve2D(p);
  const b = solve2D(p);
  assert.equal(a.center, b.center);
  assert.equal(a.wall, b.wall);
  assert.equal(a.it, b.it);
});

test("transportNumbers() returns finite dimensionless groups for both profiles", () => {
  for (const profileId of Object.keys(materialProfiles)) {
    const p = makeParams(profileId);
    const sol = solve2D(p);
    const n = transportNumbers(sol);
    for (const [key, value] of Object.entries(n)) {
      assert.ok(Number.isFinite(value), `${profileId}: transportNumbers().${key} = ${value}`);
    }
  }
});

test("solve2D() power balance holds across a small absorbed-power sweep", () => {
  const failures = [];
  for (const profileId of Object.keys(materialProfiles)) {
    for (const pabs of [5, 15, 25]) {
      const p = makeParams(profileId, { P: pabs });
      const sol = solve2D(p);
      const relBalance = Math.abs(sol.balance) / Math.max(sol.p.P, 1e-30);
      if (!sol.converged || relBalance > 0.01 || !Number.isFinite(sol.center)) {
        failures.push(`${profileId} P=${pabs}W: converged=${sol.converged}, relBalance=${relBalance}`);
      }
    }
  }
  assert.deepEqual(failures, []);
});

test("solve2D() models tube gas entering at the top and exiting at the bottom, matching the 'Gas flow' arrow drawn into the bed from above", () => {
  // Crank flow rate up so convective cooling dominates the (axially symmetric) microwave
  // source: gas is coldest right where it enters, so that end of the bed should run
  // measurably cooler than the end it exits from - the opposite pattern would mean the
  // gas-row loop is walking z backwards again.
  const p = makeParams("sic-60-100-mesh", { flow: materialProfiles["sic-60-100-mesh"].defaults.flow * 20 });
  const sol = solve2D(p);
  assert.ok(sol.converged);
  let topJ = -1, botJ = -1;
  for (let j = 0; j < p.Nz; j++) if (sol.material[j][0] === 2) { if (topJ === -1) botJ = j; topJ = j; }
  assert.ok(topJ > botJ, "expected at least one bed row found");
  const Ttop = sol.T[topJ][0], Tbot = sol.T[botJ][0];
  assert.ok(Ttop < Tbot, `top-of-bed (z high, gas inlet) should run cooler than bottom-of-bed (z low, gas outlet) under strong forced convection: Ttop=${Ttop}, Tbot=${Tbot}`);
});

test("solve2D() clamps temperature and stays finite under a wildly oversized absorbed power", () => {
  const p = makeParams("sic-60-100-mesh", { P: 5000 });
  const sol = solve2D(p);
  assert.ok(Number.isFinite(sol.center) && Number.isFinite(sol.wall) && Number.isFinite(sol.Tmax));
  // solve2D's own relaxation update clamps every cell to [Ta, 2500 C]
  assert.ok(sol.Tmax <= 2500 + 1e-6, `Tmax=${sol.Tmax} exceeded the solver's clamp ceiling`);
});

test("solve2D() with zero absorbed power settles at ambient with no NaN", () => {
  const p = makeParams("rutile-reduced-600c-30m", { P: 0 });
  const sol = solve2D(p);
  assert.ok(sol.converged, `did not converge in ${sol.it + 1} iterations, maxDelta=${sol.maxDelta}`);
  assert.ok(Math.abs(sol.center - p.Ta) < 1, `center=${sol.center} should settle near ambient=${p.Ta}`);
  assert.ok(Math.abs(sol.wall - p.Ta) < 1, `wall=${sol.wall} should settle near ambient=${p.Ta}`);
});

test("Supplementary Note 3 closures expose axisymmetric effective tensors", () => {
  const p = makeParams("rutile-reduced-600c-30m", { kzRatio: 1.5, permeabilityLongitudinalRatio: 2 });
  const closures = porousContinuumClosures(400, p);
  assert.deepEqual(closures.thermalDispersion, axisymmetricTensor(kBed(400, p) * 1.5, kBed(400, p)));
  assert.equal(closures.permeability[0][0], 2 * closures.permeability[1][1]);
  assert.ok(closures.rhoCpEffective > 0);
  assert.deepEqual(darcyVelocity([-1000, 0, 0], closures.permeability, 2e-5), [closures.permeability[0][0] * 1000 / 2e-5, 0, 0]);
});

test("homogenization validity checks both macro and microwave scale separation", () => {
  assert.deepEqual(homogenizationValidity({ unitCellLength: 1e-4, macroLength: 1e-2, wavelength: 1e-1 }), { macroRatio: 0.01, waveRatio: 0.001, valid: true });
  assert.equal(homogenizationValidity({ unitCellLength: 2e-3, macroLength: 1e-2, wavelength: 1e-1 }).valid, false);
});
