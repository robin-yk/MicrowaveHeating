// Regression tests for the microwave-heating numeric core
// (microwave-solver.js). Pure Node, no browser: run with `node --test`.
import test from "node:test";
import assert from "node:assert/strict";
import {
  materialProfiles, clamp, dielectric, kBed, penetrationDepth, solve2D, transportNumbers,
  axisymmetricTensor, homogenizationValidity, porousContinuumClosures, darcyVelocity,
  bedMesh, darcyField, darcyPermeability, packedBedTransport, gasState, ERGUN_VISCOUS_CONSTANT
} from "../microwave-solver.js";

// Bed extent on a given mesh: the row/column span of material 2, plus the axial
// length and cross-section the analytic Darcy checks below need.
function bedExtent(p, mesh) {
  let jMin = Infinity, jMax = -1, iMax = -1;
  for (let j = 0; j < p.Nz; j++) for (let i = 0; i < p.Nr; i++) if (mesh.material[j][i] === 2) {
    jMin = Math.min(jMin, j); jMax = Math.max(jMax, j); iMax = Math.max(iMax, i);
  }
  let area = 0; for (let i = 0; i <= iMax; i++) area += mesh.areasZ[i];
  return { jMin, jMax, iMax, rows: jMax - jMin + 1, length: (jMax - jMin + 1) * mesh.dz, area };
}
function uniformField(p, value) {
  return Array.from({ length: p.Nz }, () => Array(p.Nr).fill(value));
}

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
  const velocity = darcyVelocity([-1000, 0, 0], closures.permeability, 2e-5);
  assert.deepEqual(velocity, [closures.permeability[0][0] * 1000 / 2e-5, 0, 0]);
  // Pin the sign of the vanishing components: -0 breaks strict comparisons and
  // would print as "-0.00 m/s".
  assert.ok(Object.is(velocity[1], 0) && Object.is(velocity[2], 0), `transverse components must be +0, got ${velocity[1]} and ${velocity[2]}`);
});

test("closure anisotropy ratios are independent and default to isotropic", () => {
  const p = makeParams("rutile-reduced-600c-30m", { epsLongitudinalRatio: 3, kzRatio: undefined });
  const e = dielectric(400, p), closures = porousContinuumClosures(400, p);
  // Stretching eps' must leave eps'' alone: the loss ratio does not inherit the
  // real-part ratio, so the loss tangent is not silently made isotropic.
  assert.equal(closures.permittivityReal[0][0], e.ep * 3);
  assert.equal(closures.permittivityLoss[0][0], closures.permittivityLoss[1][1]);
  // An omitted kzRatio falls back to isotropic instead of producing NaN.
  assert.equal(closures.thermalDispersion[0][0], closures.thermalDispersion[1][1]);

  const lossy = porousContinuumClosures(400, makeParams("rutile-reduced-600c-30m", { epsLossLongitudinalRatio: 2 }));
  assert.equal(lossy.permittivityLoss[0][0], 2 * lossy.permittivityLoss[1][1]);
  assert.equal(lossy.permittivityReal[0][0], lossy.permittivityReal[1][1]);
  assert.throws(() => porousContinuumClosures(400, makeParams("rutile-reduced-600c-30m", { epsLossLongitudinalRatio: 0 })), RangeError);
});

test("homogenization validity checks both macro and microwave scale separation", () => {
  assert.deepEqual(homogenizationValidity({ unitCellLength: 1e-4, macroLength: 1e-2, wavelength: 1e-1 }), { macroRatio: 0.01, waveRatio: 0.001, valid: true });
  assert.equal(homogenizationValidity({ unitCellLength: 2e-3, macroLength: 1e-2, wavelength: 1e-1 }).valid, false);
});

test("darcyPermeability() uses the same constant as the Ergun viscous term", () => {
  const p = makeParams("rutile-reduced-600c-30m");
  const eps = p.voidFraction, solid = 1 - eps;
  assert.equal(ERGUN_VISCOUS_CONSTANT, 150);
  assert.equal(darcyPermeability(p), Math.pow(eps, 3) * p.dp * p.dp / (150 * solid * solid));
  // The point of pinning the constant: mu*u/K must reproduce the Ergun viscous
  // term exactly, so Darcy's law and the reported dP describe one bed.
  const state = gasState(400, p), K = darcyPermeability(p);
  const ergunViscous = 150 * state.mu * solid * solid * state.u / (Math.pow(eps, 3) * p.dp * p.dp);
  assert.ok(Math.abs(state.mu * state.u / K - ergunViscous) <= 1e-12 * ergunViscous,
    `mu*u/K = ${state.mu * state.u / K} should equal the Ergun viscous term ${ergunViscous}`);
  assert.equal(darcyPermeability({ ...p, ergunViscousConstant: 180 }), darcyPermeability(p) * 150 / 180);
  assert.throws(() => darcyPermeability({ ...p, ergunViscousConstant: 0 }), RangeError);
  // packedBedTransport must report the same permeability it drops into dP.
  assert.equal(packedBedTransport(400, p).permeability, darcyPermeability(p));
});

test("solve2D() legacy temperatures are untouched by the flow-field diagnostic", () => {
  // Stage 1 adds a Darcy solve but does not couple it to the energy equation, so
  // these must stay bit-for-bit what they were before the flow solver existed.
  // Values captured on the autoFit mesh from the commit that added them.
  const p = makeParams("rutile-reduced-600c-30m", { Nr: 15, Nz: 30, maxIter: 3500, tol: 1.5e-3, omega: 1.08 });
  const sol = solve2D(p);
  const expected = {
    center: 826.1530412131374, wall: 502.0457021673843, Tavg: 639.8001257911388,
    Tout: 552.3499012252433, qgas: 0.41140303109916077, qBoundary: 7.101922904136165,
    qrad: 18.492709840234983
  };
  for (const [key, want] of Object.entries(expected)) {
    assert.ok(Math.abs(sol[key] - want) <= 1e-9 * Math.abs(want),
      `${key} drifted: got ${sol[key]}, expected ${want}`);
  }
  assert.equal(sol.it, 107);
});

test("darcyField() reproduces the analytic Darcy pressure drop on an isothermal bed", () => {
  const p = makeParams("rutile-reduced-600c-30m");
  const mesh = bedMesh(p), extent = bedExtent(p, mesh), T = 300;
  const field = darcyField({ p, T: uniformField(p, T), material: mesh.material, dr: mesh.dr, dz: mesh.dz, areasZ: mesh.areasZ });
  assert.ok(field.converged);
  const state = gasState(T, p), K = darcyPermeability(p);
  // Uniform mobility makes the discrete solution exact, so dP = mu*u*L/K holds to
  // round-off rather than to discretization error.
  const u = field.massFlow / (state.rho * extent.area);
  const analytic = state.mu * u * extent.length / K;
  assert.ok(Math.abs(field.dP - analytic) <= 1e-12 * analytic,
    `dP=${field.dP} should match mu*u*L/K=${analytic}`);
  // An isothermal bed cannot redistribute flow: every column carries the same
  // superficial velocity and nothing crosses a radial face.
  assert.ok(Math.abs(field.maldistribution - 1) < 1e-12, `maldistribution=${field.maldistribution}`);
  let maxRadial = 0;
  for (let j = 0; j < p.Nz; j++) for (let i = 0; i <= p.Nr; i++) maxRadial = Math.max(maxRadial, Math.abs(field.fluxR[j][i]));
  assert.ok(maxRadial <= 1e-12 * field.massFlow, `radial flux ${maxRadial} should vanish on an isothermal bed`);
});

test("darcyField() conserves mass discretely and delivers the prescribed flow", () => {
  const p = makeParams("rutile-reduced-600c-30m");
  const mesh = bedMesh(p), extent = bedExtent(p, mesh);
  // Temperature varying in both r and z, so mobility varies axially and the
  // radial coupling in the pressure equation is actually exercised.
  const T = Array.from({ length: p.Nz }, (_, j) => Array.from({ length: p.Nr }, (_, i) => {
    if (mesh.material[j][i] !== 2) return 300;
    const rn = (i + .5) * mesh.dr / (p.D / 2), zn = (j - extent.jMin) / Math.max(extent.rows - 1, 1);
    return 200 + 700 * (1 - rn) * (1 - Math.abs(2 * zn - 1));
  }));
  const field = darcyField({ p, T, material: mesh.material, dr: mesh.dr, dz: mesh.dz, areasZ: mesh.areasZ });
  assert.ok(field.converged && field.it > 0, `expected the pressure solve to iterate, it=${field.it}`);
  assert.ok(field.maxMassImbalance <= 1e-10 * field.massFlow,
    `per-cell mass imbalance ${field.maxMassImbalance} must vanish for the energy equation to be conservative`);
  assert.ok(Math.abs(field.outletMassFlow - field.massFlow) <= 1e-10 * field.massFlow,
    `outlet flow ${field.outletMassFlow} should equal the prescribed ${field.massFlow}`);
  let maxRadial = 0;
  for (let j = 0; j < p.Nz; j++) for (let i = 0; i <= p.Nr; i++) maxRadial = Math.max(maxRadial, Math.abs(field.fluxR[j][i]));
  assert.ok(maxRadial > 1e-4 * field.massFlow, `axially varying mobility should drive radial flow, got ${maxRadial}`);
  // Warm starting from the converged potential must land on the same answer.
  const warm = darcyField({ p, T, material: mesh.material, dr: mesh.dr, dz: mesh.dz, areasZ: mesh.areasZ, phi0: field.phi });
  assert.ok(Math.abs(warm.dP - field.dP) <= 1e-9 * field.dP, `warm start dP=${warm.dP} vs cold ${field.dP}`);
});

test("darcyField() pushes flow away from a hot core and scales with axial permeability", () => {
  const p = makeParams("rutile-reduced-600c-30m");
  const mesh = bedMesh(p), extent = bedExtent(p, mesh);
  // Hot on the axis, cool at the wall: helium viscosity rises with temperature,
  // so the central columns resist more and carry less flow. A single plug-flow
  // stream cannot represent this, which is the reason for solving the field.
  const T = Array.from({ length: p.Nz }, (_, j) => Array.from({ length: p.Nr }, (_, i) =>
    mesh.material[j][i] === 2 ? 800 - 600 * ((i + .5) * mesh.dr / (p.D / 2)) : 300));
  const field = darcyField({ p, T, material: mesh.material, dr: mesh.dr, dz: mesh.dz, areasZ: mesh.areasZ });
  assert.ok(field.converged);
  assert.ok(field.columnVelocity[0] < field.columnVelocity[extent.iMax],
    `hot axis should run slower than the cool wall: ${field.columnVelocity[0]} vs ${field.columnVelocity[extent.iMax]}`);
  assert.ok(field.maldistribution > 1.1, `expected real maldistribution, got ${field.maldistribution}`);
  // Doubling the longitudinal permeability halves the pressure drop.
  const stretched = darcyField({ p: { ...p, permeabilityLongitudinalRatio: 2 }, T, material: mesh.material, dr: mesh.dr, dz: mesh.dz, areasZ: mesh.areasZ });
  assert.ok(Math.abs(stretched.dP / field.dP - 0.5) < 1e-9, `dP ratio=${stretched.dP / field.dP}, expected 0.5`);
  assert.throws(() => darcyField({ p: { ...p, permeabilityLongitudinalRatio: -1 }, T, material: mesh.material, dr: mesh.dr, dz: mesh.dz, areasZ: mesh.areasZ }), RangeError);
});

test("darcyField() degrades gracefully at zero flow and can be switched off", () => {
  const p = makeParams("rutile-reduced-600c-30m", { flow: 0 });
  const mesh = bedMesh(p);
  const field = darcyField({ p, T: uniformField(p, 300), material: mesh.material, dr: mesh.dr, dz: mesh.dz, areasZ: mesh.areasZ });
  for (const key of ["dP", "massFlow", "outletMassFlow", "uMax", "uMin", "maxMassImbalance"]) {
    assert.equal(field[key], 0, `${key} should be exactly zero with no flow`);
  }
  assert.equal(field.maldistribution, 1);
  const off = solve2D(makeParams("rutile-reduced-600c-30m", { P: 0, flowMode: "off" }));
  assert.equal(off.darcy, null);
  const on = solve2D(makeParams("rutile-reduced-600c-30m", { P: 0 }));
  assert.ok(on.darcy && Number.isFinite(on.darcy.dP), "the flow field should be reported by default");
});
