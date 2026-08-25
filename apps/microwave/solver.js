// Microwave heating numeric core: dielectric mixing, bed conductivity, gas
// transport, and the 2D axisymmetric FVM steady-state thermal solver behind
// microwave.html. Pure functions of plain input objects, with no DOM access,
// so this module can be imported directly by microwave.html (as an ES
// module) or by a Node test runner.
"use strict";

export const sigma = 5.670374419e-8;
export const c0 = 299792458;

export const clamp = (x, a, b) => Math.max(a, Math.min(b, x));

export const materialProfiles = {
  "rutile-reduced-600c-30m": {
    label: "Reduced rutile TiO₂ · H₂ 600 °C, 30 min", formula: "TiO₂", rhoSolid: 4230, calibrationLabel: "Embedded reduced-rutile TiO₂ power sweep", dielectricNote: "Reduced-rutile ε′(T), ε″(T) values are an approximate digitization of the supplied curve; replace them with the numerical source before publication use.", fitBounds: { k200: [.2, 3], k500: [.15, 2], k800: [.1, 1.5], radArea: [.3, 3] },
    dielectric: `20,3.49,0.3204
100,5.8,3.0
200,6.2,4.0
400,6.5,4.8
600,6.7,5.4
725,7.2,6.2`,
    experiments: `0,20,20
5,167.3,185.5
10,287.1,354.6
14,360,490
17,402.7,580
20,425.5,618
23,468,710
26,499.6,799.2`,
    cpSolid: 690,
    defaults: { pabs: 26, frequency: 2.404, volume: 1.18, diameter: 10, length: 15, mass: 1.150, gas: "He", flow: 50, ambient: 20, "tube-thickness": 1, "gas-pressure": 1, "particle-diameter": 50, "fbg-r": 0, "fbg-z": 0, k200: 1.70, k500: .62, k800: .58, "kz-ratio": 1.20, "h-contact": 1000, "k-quartz": 1.40, "air-factor": 3, "boundary-mode": "automatic", "h-boundary": 12, emissivity: .85, "rad-area": 1.80, "gas-transfer-mode": "automatic", "gas-eff": .80, "dielectric-mode": "looyenga", "bed-k-mode": "automatic", "field-wr": 1.20, "field-wz": 1.20 }
  },
  "sic-60-100-mesh": {
    label: "SiC 60–100 mesh · as received", formula: "SiC", rhoSolid: 3210, calibrationLabel: "Embedded SiC 60–100 mesh power sweep · T_wall at quartz outer surface", dielectricNote: "SiC ε′(T), ε″(T) values are 20–800 °C half-degree-bin averages extracted from Sheet1 at approximately 2.404 GHz.", fitBounds: { k200: [.5, 30], k500: [.5, 30], k800: [.5, 30], radArea: [1, 12] },
    dielectric: `20,7.959225,0.398952
25,7.969336,0.385909
50,8.005840,0.344083
100,8.056020,0.269573
150,8.124843,0.219999
200,8.175771,0.181880
250,8.217561,0.156341
300,8.258681,0.138804
350,8.298917,0.126258
400,8.340075,0.118962
450,8.382080,0.113153
500,8.433152,0.111549
550,8.482923,0.110614
600,8.542403,0.112061
650,8.588825,0.115461
700,8.655626,0.117419
750,8.731564,0.124630
800,8.834589,0.132493`,
    experiments: `0,20,20
10,185.5,225
20,275.3,340.5
30,331,411
40,395,480`,
    cpSolid: 750,
    defaults: { pabs: 40, frequency: 2.404, volume: 1.18, diameter: 10, length: 15, mass: 1.628, gas: "He", flow: 65.3, ambient: 20, "tube-thickness": 1, "gas-pressure": 1, "particle-diameter": 194, "fbg-r": 0, "fbg-z": 0, k200: 4.00, k500: 18.00, k800: 18.00, "kz-ratio": 1.00, "h-contact": 1500, "k-quartz": 1.40, "air-factor": 3, "boundary-mode": "automatic", "h-boundary": 12, emissivity: .85, "rad-area": 6.00, "gas-transfer-mode": "automatic", "gas-eff": .80, "dielectric-mode": "looyenga", "bed-k-mode": "automatic", "field-wr": 1.20, "field-wz": 1.20 }
  }
};

export function parseRows(text, n) {
  return text.trim().split(/\n+/).map(line => line.trim().split(/[\t, ]+/).map(Number)).filter(row => row.length >= n && row.slice(0, n).every(Number.isFinite)).sort((a, b) => a[0] - b[0]);
}
export function interpolate(x, rows, col) {
  if (!rows.length) return NaN;
  if (x <= rows[0][0]) return rows[0][col];
  if (x >= rows[rows.length - 1][0]) return rows[rows.length - 1][col];
  for (let i = 1; i < rows.length; i++) if (x <= rows[i][0]) return rows[i - 1][col] + (rows[i][col] - rows[i - 1][col]) * (x - rows[i - 1][0]) / (rows[i][0] - rows[i - 1][0]);
  return rows[rows.length - 1][col];
}
export function complexCubeRoot(z) { const mag = Math.hypot(z.re, z.im), angle = Math.atan2(z.im, z.re) / 3, root = Math.cbrt(mag); return { re: root * Math.cos(angle), im: root * Math.sin(angle) }; }
export function complexCube(z) { return { re: z.re * z.re * z.re - 3 * z.re * z.im * z.im, im: 3 * z.re * z.re * z.im - z.im * z.im * z.im }; }
export function dielectricReference(T, p) { return { ep: interpolate(T, p.diel, 1), epp: Math.max(1e-9, interpolate(T, p.diel, 2)) }; }
export function dielectric(T, p) {
  const ref = dielectricReference(T, p); if (p.dielectricMode === "manual") return ref;
  const eRef = complexCubeRoot({ re: ref.ep, im: -ref.epp }), gasRoot = { re: 1, im: 0 }, solidRoot = { re: (eRef.re - p.referenceVoidFraction * gasRoot.re) / (1 - p.referenceVoidFraction), im: (eRef.im - p.referenceVoidFraction * gasRoot.im) / (1 - p.referenceVoidFraction) };
  const mixed = complexCube({ re: (1 - p.voidFraction) * solidRoot.re + p.voidFraction * gasRoot.re, im: (1 - p.voidFraction) * solidRoot.im + p.voidFraction * gasRoot.im });
  return { ep: Math.max(1.000001, mixed.re), epp: Math.max(1e-9, -mixed.im) };
}
export function kBedReference(T, p) { return Math.exp(interpolate(T, [[200, Math.log(p.k200)], [500, Math.log(p.k500)], [800, Math.log(p.k800)]], 1)); }
export function gasConductivity(temp, gas) { const k0 = gas === "He" ? .1513 : .0258; return k0 * Math.pow((temp + 273.15) / 293.15, .72); }
export function maxwellEucken(kSkeleton, kGas, voidFraction) { return kSkeleton * (2 * kSkeleton + kGas - 2 * voidFraction * (kSkeleton - kGas)) / (2 * kSkeleton + kGas + voidFraction * (kSkeleton - kGas)); }
export function inferredSkeletonK(kReference, kGas, referenceVoid) {
  if (kReference <= kGas) return Math.max(kReference, 1e-9);
  let lo = Math.max(kGas, 1e-9), hi = Math.max(kReference * 4, kGas * 8);
  while (maxwellEucken(hi, kGas, referenceVoid) < kReference && hi < 1e6) hi *= 2;
  for (let i = 0; i < 70; i++) { const mid = (lo + hi) / 2; if (maxwellEucken(mid, kGas, referenceVoid) < kReference) lo = mid; else hi = mid; }
  return (lo + hi) / 2;
}
export function kBed(T, p) {
  const kReference = kBedReference(T, p); if (p.bedKMode === "manual") return kReference;
  const kGas = gasConductivity(T, p.gas), kSkeleton = inferredSkeletonK(kReference, kGas, p.referenceVoidFraction);
  return Math.max(1e-6, maxwellEucken(kSkeleton, kGas, p.voidFraction));
}
export function penetrationDepth(T, p) {
  const e = dielectric(T, p), k0 = 2 * Math.PI * p.frequency / c0;
  const alpha = k0 * Math.sqrt(Math.max(0, (Math.hypot(e.ep, e.epp) - e.ep) / 2));
  return alpha > 1e-15 ? 1 / (2 * alpha) : 1e9;
}
export function gasCapacityRate(p) {
  const cp = p.gas === "He" ? 20.786 : 29.124;
  return p.flow * 1e-6 / 60 / 0.022414 * cp;
}
export function gasState(tempC, p, pressure = p.pressure) {
  const Rg = 8.314462618, TK = tempC + 273.15, M = p.gas === "He" ? .004002602 : .0280134, Cpm = p.gas === "He" ? 20.786 : 29.124, mu0 = p.gas === "He" ? 1.96e-5 : 1.76e-5, S = p.gas === "He" ? 79.4 : 111;
  const mu = mu0 * Math.pow(TK / 293.15, 1.5) * (293.15 + S) / (TK + S), rho = pressure * M / (Rg * TK), cp = Cpm / M, kg = gasConductivity(tempC, p.gas), Pr = mu * cp / kg, ndot = p.flow * 1e-6 / 60 / .022414, Qactual = ndot * Rg * TK / pressure, u = Qactual / (Math.PI * p.D * p.D / 4);
  return { TK, M, Cpm, mu, rho, cp, kg, Pr, ndot, Qactual, u, pressure };
}
// Blake-Kozeny permeability. Its constant has to be the same one the Ergun
// viscous term uses, or Darcy's law and the reported pressure drop describe
// two different beds: at 180 here against Ergun's 150 below, mu*u/K came out
// exactly 20% above the viscous term it is supposed to reproduce. The Ergun
// value wins because dP is what the page reports and what experiments measure.
export const ERGUN_VISCOUS_CONSTANT = 150;
export function darcyPermeability(p) {
  const eps = p.voidFraction, solid = Math.max(1e-8, 1 - eps), constant = p.ergunViscousConstant ?? ERGUN_VISCOUS_CONSTANT;
  if (!(Number.isFinite(constant) && constant > 0)) throw new RangeError("ergunViscousConstant must be finite and positive");
  return Math.pow(eps, 3) * p.dp * p.dp / (constant * solid * solid);
}
export function packedBedTransport(tempC, p) {
  const g = 9.80665, eps = p.voidFraction, solid = Math.max(1e-8, 1 - eps), Cg = gasCapacityRate(p); let pressure = p.pressure, state, dP = 0, ReP = 0, NuP = 2, hgs = 0, permeability = 0;
  for (let pass = 0; pass < 3; pass++) {
    state = gasState(tempC, p, pressure); ReP = state.rho * state.u * p.dp / state.mu; NuP = 2 + 1.1 * Math.pow(Math.max(ReP, 0), .6) * Math.pow(state.Pr, 1 / 3); hgs = NuP * state.kg / p.dp;
    permeability = darcyPermeability(p);
    // Viscous term written as mu*u/K so the identity with Darcy's law holds by
    // construction instead of by two constants happening to agree.
    const dPdz = state.mu * state.u / permeability + 1.75 * state.rho * solid * state.u * state.u / (Math.pow(eps, 3) * p.dp);
    dP = dPdz * p.H; pressure = Math.max(5000, p.pressure - dP / 2);
  }
  const specificArea = 6 * solid / p.dp, UA = hgs * specificArea * p.volume * 1e-6, gasEffectiveness = Cg > 0 ? 1 - Math.exp(-UA / Cg) : 0, rhoCpEff = solid * p.rhoSolid * p.cpSolid + eps * state.rho * state.cp;
  return { ...state, Cg, ReP, NuP, hgs, specificArea, UA, gasEffectiveness, permeability, dP, outletPressure: Math.max(5000, p.pressure - dP), rhoCpEff };
}

// Constitutive data shared by a future coupled Maxwell/Darcy/LTE/species solve.
// Components follow the reactor-axis convention [longitudinal, transverse,
// transverse] used in Supplementary Note 3.
export function axisymmetricTensor(longitudinal, transverse) {
  if (![longitudinal, transverse].every(value => Number.isFinite(value) && value > 0)) throw new RangeError("tensor components must be finite and positive");
  return [[longitudinal, 0, 0], [0, transverse, 0], [0, 0, transverse]];
}

export function homogenizationValidity({ unitCellLength, macroLength, wavelength, maxRatio = 0.1 }) {
  if (![unitCellLength, macroLength, wavelength, maxRatio].every(value => Number.isFinite(value) && value > 0)) throw new RangeError("homogenization lengths and maxRatio must be finite and positive");
  const macroRatio = unitCellLength / macroLength, waveRatio = unitCellLength / wavelength;
  return { macroRatio, waveRatio, valid: macroRatio < maxRatio && waveRatio < maxRatio };
}

// Every anisotropy ratio is supplied independently and defaults to isotropic.
// The loss part in particular carries its own ratio: eps' and eps'' come from
// separate measurements, so scaling eps'' by the eps' ratio would silently
// assert a direction-independent loss tangent that no measurement supports.
export function porousContinuumClosures(tempC, p) {
  const e = dielectric(tempC, p), transport = packedBedTransport(tempC, p), kr = kBed(tempC, p);
  const epsLongitudinalRatio = p.epsLongitudinalRatio ?? 1, epsLossLongitudinalRatio = p.epsLossLongitudinalRatio ?? 1,
    permeabilityLongitudinalRatio = p.permeabilityLongitudinalRatio ?? 1, thermalLongitudinalRatio = p.kzRatio ?? 1;
  if (![epsLongitudinalRatio, epsLossLongitudinalRatio, permeabilityLongitudinalRatio, thermalLongitudinalRatio].every(value => Number.isFinite(value) && value > 0)) throw new RangeError("anisotropy ratios must be finite and positive");
  return {
    permittivityReal: axisymmetricTensor(e.ep * epsLongitudinalRatio, e.ep),
    permittivityLoss: axisymmetricTensor(e.epp * epsLossLongitudinalRatio, e.epp),
    permeability: axisymmetricTensor(transport.permeability * permeabilityLongitudinalRatio, transport.permeability),
    thermalDispersion: axisymmetricTensor(kr * thermalLongitudinalRatio, kr),
    rhoCpEffective: transport.rhoCpEff
  };
}

export function darcyVelocity(pressureGradient, permeabilityTensor, viscosity) {
  if (!Array.isArray(pressureGradient) || pressureGradient.length !== 3 || !pressureGradient.every(Number.isFinite)) throw new TypeError("pressureGradient must contain three finite components");
  if (!Number.isFinite(viscosity) || viscosity <= 0) throw new RangeError("viscosity must be finite and positive");
  if (!Array.isArray(permeabilityTensor) || permeabilityTensor.length !== 3 || permeabilityTensor.some(row => !Array.isArray(row) || row.length !== 3 || !row.every(Number.isFinite))) throw new TypeError("permeabilityTensor must be a finite 3 by 3 matrix");
  // A vanishing component must read as +0: negating a zero sum yields -0, which
  // fails strict comparisons and would surface as "-0.00 m/s" once a pressure
  // solver feeds these velocities to the page.
  return permeabilityTensor.map(row => {
    const component = -row.reduce((sum, value, column) => sum + value * pressureGradient[column], 0) / viscosity;
    return component === 0 ? 0 : component;
  });
}
export function naturalConvection(tempC, p) {
  const Rg = 8.314462618, g = 9.80665, Tfilm = (tempC + p.Ta) / 2 + 273.15, Mair = .0289652, mu = 1.81e-5 * Math.pow(Tfilm / 293.15, 1.5) * (293.15 + 111) / (Tfilm + 111), rho = p.pressure * Mair / (Rg * Tfilm), cp = 1007, k = .0263 * Math.pow(Tfilm / 293.15, .76), Pr = mu * cp / k, nu = mu / rho, beta = 1 / Tfilm, L = Math.max(p.domainHeight, 1e-6), dT = Math.abs(tempC - p.Ta), Gr = g * beta * dT * Math.pow(L, 3) / (nu * nu), Ra = Gr * Pr, Nu = Math.pow(.825 + .387 * Math.pow(Math.max(Ra, 0), 1 / 6) / Math.pow(1 + Math.pow(.492 / Pr, 9 / 16), 8 / 27), 2);
  return { Tfilm, Pr, Gr, Ra, Nu, h: Nu * k / L };
}

// Mesh, material map, and cell metrics for the axisymmetric domain. Extracted
// so the flow solver and its tests build the exact same grid solve2D uses:
// material 2 = packed bed, 1 = tube gas, 3 = quartz, 0 = outside air.
export function bedMesh(p) {
  const R = p.D / 2, Ro = R + p.tq, Rd = p.domainWidth / 2, Hd = p.domainHeight, dr = Rd / p.Nr, dz = Hd / p.Nz;
  const material = Array.from({ length: p.Nz }, (_, j) => Array.from({ length: p.Nr }, (_, i) => {
    const r = (i + .5) * dr, z = (j + .5) * dz - Hd / 2;
    if (r < R && Math.abs(z) < p.H / 2) return 2;
    if (r < R) return 1;
    if (r < Ro) return 3;
    return 0;
  }));
  const vols = [], areasZ = [];
  for (let i = 0; i < p.Nr; i++) { const rw = i * dr, re = (i + 1) * dr; vols[i] = Math.PI * (re * re - rw * rw) * dz; areasZ[i] = Math.PI * (re * re - rw * rw); }
  return { R, Ro, Rd, Hd, dr, dz, material, vols, areasZ };
}

// Superficial-velocity field on the packed bed from continuity plus Darcy's
// law, div(lambda grad P) = 0 with mobility lambda = rho*K/mu.
//
// Both bed end faces are isobars: the open tube above and below the bed adds no
// appreciable resistance, so flow splits between columns purely by local
// mobility. That radial maldistribution -- hot centre means higher mu, so less
// flow -- is the physics a single plug-flow stream cannot express, and it is
// exactly what a prescribed uniform inlet flux would erase.
//
// The system is linear in P, so this solves once with phi = 1 on the inlet face
// and phi = 0 on the outlet face, measures the unit-potential mass flow, then
// scales by dP = mdot / flowPerPascal. Face mass fluxes are read straight out of
// that solution, which makes the discrete mass balance exact rather than
// approximate -- a requirement for the energy equation Stage 2 builds on top.
//
// Inlet and outlet faces are found from the material map (a bed cell with no bed
// neighbour above or below) rather than assumed rectangular, so the geometry
// stays free to change. Gas enters at the top, matching solve2D's gas march.
export function darcyField({ p, T, material, dr, dz, areasZ, phi0 = null, maxIter = 20000, tol = 1e-13, omega = 1.7 }) {
  const Nz = material.length, Nr = material[0].length;
  const K = darcyPermeability(p), ratio = p.permeabilityLongitudinalRatio ?? 1;
  if (!(Number.isFinite(ratio) && ratio > 0)) throw new RangeError("permeabilityLongitudinalRatio must be finite and positive");
  if (!(Number.isFinite(omega) && omega > 0 && omega < 2)) throw new RangeError("omega must lie in (0, 2)");
  const lamZ = Array.from({ length: Nz }, () => Array(Nr).fill(0)), lamR = Array.from({ length: Nz }, () => Array(Nr).fill(0)), rho = Array.from({ length: Nz }, () => Array(Nr).fill(0));
  const bed = (j, i) => j >= 0 && j < Nz && i >= 0 && i < Nr && material[j][i] === 2;
  for (let j = 0; j < Nz; j++) for (let i = 0; i < Nr; i++) {
    if (!bed(j, i)) continue;
    const g = gasState(T[j][i], p);
    rho[j][i] = g.rho; lamZ[j][i] = g.rho * K * ratio / g.mu; lamR[j][i] = g.rho * K / g.mu;
  }
  // Harmonic-mean face mobilities, mirroring how interfaceG adds resistances.
  const gR = (j, i) => { const A = 2 * Math.PI * (i + 1) * dr * dz; return 1 / (dr / (2 * lamR[j][i] * A) + dr / (2 * lamR[j][i + 1] * A)); };
  const gZ = (j, i) => { const A = areasZ[i]; return 1 / (dz / (2 * lamZ[j][i] * A) + dz / (2 * lamZ[j + 1][i] * A)); };
  const gEnd = (j, i) => lamZ[j][i] * areasZ[i] / (dz / 2);
  const isInlet = (j, i) => bed(j, i) && !bed(j + 1, i), isOutlet = (j, i) => bed(j, i) && !bed(j - 1, i);
  // Seed with the linear profile, which is the exact answer for a uniform bed,
  // so a cold start costs almost nothing and a warm start costs less.
  const phi = phi0 ? phi0.map(row => row.slice()) : Array.from({ length: Nz }, () => Array(Nr).fill(0));
  if (!phi0) for (let i = 0; i < Nr; i++) {
    const rows = []; for (let j = 0; j < Nz; j++) if (bed(j, i)) rows.push(j);
    rows.forEach((j, n) => { phi[j][i] = (n + .5) / rows.length; });
  }
  let it = 0, converged = false, maxDelta = Infinity;
  for (it = 0; it < maxIter; it++) {
    maxDelta = 0;
    for (let j = 0; j < Nz; j++) for (let i = 0; i < Nr; i++) {
      if (!bed(j, i)) continue;
      let sum = 0, rhs = 0;
      if (bed(j, i - 1)) { const G = gR(j, i - 1); sum += G; rhs += G * phi[j][i - 1]; }
      if (bed(j, i + 1)) { const G = gR(j, i); sum += G; rhs += G * phi[j][i + 1]; }
      if (bed(j - 1, i)) { const G = gZ(j - 1, i); sum += G; rhs += G * phi[j - 1][i]; } else sum += gEnd(j, i);
      if (bed(j + 1, i)) { const G = gZ(j, i); sum += G; rhs += G * phi[j + 1][i]; } else { const G = gEnd(j, i); sum += G; rhs += G; }
      const next = rhs / Math.max(sum, 1e-300), updated = phi[j][i] + omega * (next - phi[j][i]);
      maxDelta = Math.max(maxDelta, Math.abs(updated - phi[j][i])); phi[j][i] = updated;
    }
    if (maxDelta < tol) { converged = true; break; }
  }
  const inlet = gasState(p.Ta, p), massFlow = inlet.ndot * inlet.M;
  let flowPerPascal = 0;
  for (let j = 0; j < Nz; j++) for (let i = 0; i < Nr; i++) if (isOutlet(j, i)) flowPerPascal += gEnd(j, i) * phi[j][i];
  const dP = flowPerPascal > 0 ? massFlow / flowPerPascal : 0;
  // Face mass flows in kg/s, signed along +z and +r. fluxZ[j][i] crosses the
  // face below cell j; fluxR[j][i] crosses the face at r = i*dr.
  const fluxZ = Array.from({ length: Nz + 1 }, () => Array(Nr).fill(0)), fluxR = Array.from({ length: Nz }, () => Array(Nr + 1).fill(0));
  for (let j = 0; j < Nz; j++) for (let i = 0; i < Nr; i++) {
    if (!bed(j, i)) continue;
    if (bed(j - 1, i)) fluxZ[j][i] = dP * gZ(j - 1, i) * (phi[j - 1][i] - phi[j][i]);
    else fluxZ[j][i] = dP * gEnd(j, i) * (0 - phi[j][i]);
    if (isInlet(j, i)) fluxZ[j + 1][i] = dP * gEnd(j, i) * (phi[j][i] - 1);
    if (bed(j, i - 1)) fluxR[j][i] = dP * gR(j, i - 1) * (phi[j][i - 1] - phi[j][i]);
  }
  let maxMassImbalance = 0, outletMassFlow = 0, uMax = 0, uMin = Infinity;
  const columnVelocity = Array(Nr).fill(0);
  for (let j = 0; j < Nz; j++) for (let i = 0; i < Nr; i++) {
    if (!bed(j, i)) continue;
    maxMassImbalance = Math.max(maxMassImbalance, Math.abs(fluxZ[j][i] - fluxZ[j + 1][i] + fluxR[j][i] - fluxR[j][i + 1]));
    if (!isOutlet(j, i)) continue;
    outletMassFlow += -fluxZ[j][i];
    const u = rho[j][i] > 0 ? -fluxZ[j][i] / (rho[j][i] * areasZ[i]) : 0;
    columnVelocity[i] = u; uMax = Math.max(uMax, u); uMin = Math.min(uMin, u);
  }
  if (!Number.isFinite(uMin)) uMin = 0;
  return { phi, dP, permeability: K, inletPressure: p.pressure + dP, outletPressure: p.pressure, massFlow, outletMassFlow, maxMassImbalance, fluxZ, fluxR, columnVelocity, uMax, uMin, maldistribution: uMin > 0 ? uMax / uMin : 1, it, converged, maxDelta };
}

export function solve2D(p) {
  const mesh = bedMesh(p), { R, Ro, Rd, Hd, dr, dz, material, vols, areasZ } = mesh;
  const T = Array.from({ length: p.Nz }, (_, j) => Array.from({ length: p.Nr }, (_, i) => {
    const r = (i + .5) * dr, z = (j + .5) * dz - Hd / 2, d = Math.hypot(Math.max(0, r - R), Math.max(0, Math.abs(z) - p.H / 2));
    return p.Ta + p.P * 22 * Math.exp(-d / .006);
  }));
  let V = 0;
  for (let j = 0; j < p.Nz; j++) for (let i = 0; i < p.Nr; i++) if (material[j][i] === 2) V += vols[i];
  const airK = temp => .0263 * Math.pow((temp + 273.15) / 293.15, .76) * p.airFactor;
  const conductivity = (m, temp, dir) => m === 2 ? (dir === "z" ? p.kzRatio : 1) * kBed(temp, p) : m === 1 ? gasConductivity(temp, p.gas) : m === 3 ? p.kq : airK(temp);
  const interfaceG = (i, j, ni, nj, area, spacing, dir) => {
    const m1 = material[j][i], m2 = material[nj][ni], k1 = conductivity(m1, T[j][i], dir), k2 = conductivity(m2, T[nj][ni], dir); let resistance = spacing / (2 * k1 * area) + spacing / (2 * k2 * area);
    if ((m1 === 2 && m2 === 3) || (m1 === 3 && m2 === 2)) resistance += 1 / (p.hContact * area);
    return 1 / resistance;
  };
  const heat = Array.from({ length: p.Nz }, () => Array(p.Nr).fill(0)), gasIn = Array(p.Nz).fill(p.Ta), gasG = Array(p.Nz).fill(0), Cg = gasCapacityRate(p), sampleRows = [];
  // Tube gas enters at the top (highest z, matching the "Gas flow" arrow drawn into the
  // bed from above) and exits at the bottom, so this must walk j from high to low: the
  // loops below start Tin at ambient and carry it row-to-row in sampleRows order.
  for (let j = p.Nz - 1; j >= 0; j--) if (material[j][0] === 2) sampleRows.push(j);
  let converged = false, maxDelta = Infinity, it = 0, lastTransport = packedBedTransport(p.Ta, p), gasEffectivenessUsed = p.gasTransferMode === "manual" ? p.gasEff : lastTransport.gasEffectiveness, hBoundaryEff = p.boundaryMode === "manual" ? p.hBoundary : naturalConvection(p.Ta, p).h;
  // p.fieldMode "helmholtz" replaces the fitted Gaussian-times-Beer-Lambert
  // shape with the solved frequency-domain field. Off by default: switching it
  // on changes every number the page reports, and the bed conductivities on the
  // Calibration tab were fitted against the old shape, so they have to be
  // refitted alongside it rather than carried over.
  //
  // eps''(T) moves with temperature, so the field is stale as soon as T
  // advances. It is refreshed every fieldEvery sweeps rather than every sweep:
  // the field costs a Krylov solve where a temperature sweep costs two
  // Gauss-Seidel passes, and eps'' varies slowly enough over 25 sweeps that the
  // final refresh below settles it.
  let solvedField = null, fieldSolve = null;
  const refreshField = () => {
    if (p.fieldMode !== "helmholtz") return;
    fieldSolve = solveField2D({ p, T, mesh });
    solvedField = fieldSolve;
  };
  refreshField();
  for (it = 0; it < p.maxIter; it++) {
    if (solvedField && it > 0 && it % (p.fieldEvery ?? 25) === 0) refreshField();
    let bedTemp = 0, bedTempV = 0; for (let j = 0; j < p.Nz; j++) for (let i = 0; i < p.Nr; i++) if (material[j][i] === 2) { bedTemp += T[j][i] * vols[i]; bedTempV += vols[i]; } bedTemp /= Math.max(bedTempV, 1e-30);
    let wallIndex = clamp(Math.floor(R / dr), 0, p.Nr - 1); for (let i = 0; i < p.Nr; i++) if (material[Math.floor(p.Nz / 2)][i] === 3) wallIndex = i; const wallGuess = T[Math.floor(p.Nz / 2)][wallIndex];
    lastTransport = packedBedTransport((bedTemp + p.Ta) / 2, p); gasEffectivenessUsed = p.gasTransferMode === "manual" ? p.gasEff : lastTransport.gasEffectiveness;
    const UAg = p.gasTransferMode === "manual" ? (Cg > 0 ? -Math.log(Math.max(1e-8, 1 - gasEffectivenessUsed)) * Cg : 0) : lastTransport.UA;
    hBoundaryEff = p.boundaryMode === "manual" ? p.hBoundary : naturalConvection(wallGuess, p).h;
    let norm = 0;
    for (let j = 0; j < p.Nz; j++) for (let i = 0; i < p.Nr; i++) {
      if (material[j][i] !== 2) { heat[j][i] = 0; continue; }
      if (solvedField) { heat[j][i] = solvedField.heat[j * p.Nr + i] * vols[i]; norm += heat[j][i]; continue; }
      const r = (i + .5) * dr, z = (j + .5) * dz - Hd / 2, e = dielectric(T[j][i], p), field = Math.exp(-Math.pow(r / (Math.max(.05, p.fieldWr) * R), 2) - Math.pow(z / (Math.max(.05, p.fieldWz) * p.H / 2), 2)), attenuation = Math.exp(-(R - r) / penetrationDepth(T[j][i], p));
      heat[j][i] = e.epp * field * attenuation * vols[i]; norm += heat[j][i];
    }
    for (let j = 0; j < p.Nz; j++) for (let i = 0; i < p.Nr; i++) heat[j][i] = p.P * heat[j][i] / Math.max(norm, 1e-30);
    let Tin = p.Ta;
    for (const j of sampleRows) {
      let sliceV = 0, sliceT = 0; for (let i = 0; i < p.Nr; i++) if (material[j][i] === 2) { sliceV += vols[i]; sliceT += vols[i] * T[j][i]; } sliceT /= sliceV;
      gasIn[j] = Tin; const eff = Cg > 0 ? 1 - Math.exp(-(UAg / sampleRows.length) / Cg) : 0; gasG[j] = Cg * eff; Tin += gasG[j] * (sliceT - Tin) / Math.max(Cg, 1e-30);
    }
    maxDelta = 0;
    for (let sweep = 0; sweep < 2; sweep++) for (let j = 0; j < p.Nz; j++) for (let i = 0; i < p.Nr; i++) {
      const rw = i * dr, re = (i + 1) * dr, Aw = 2 * Math.PI * rw * dz, Ae = 2 * Math.PI * re * dz, Az = areasZ[i], m = material[j][i]; let sum = 0, rhs = heat[j][i];
      if (i > 0) { const G = interfaceG(i, j, i - 1, j, Aw, dr, "r"); sum += G; rhs += G * T[j][i - 1]; }
      if (i < p.Nr - 1) { const G = interfaceG(i, j, i + 1, j, Ae, dr, "r"); sum += G; rhs += G * T[j][i + 1]; } else { const G = hBoundaryEff * Ae; sum += G; rhs += G * p.Ta; }
      if (j > 0) { const G = interfaceG(i, j, i, j - 1, Az, dz, "z"); sum += G; rhs += G * T[j - 1][i]; } else { const G = hBoundaryEff * Az; sum += G; rhs += G * p.Ta; }
      if (j < p.Nz - 1) { const G = interfaceG(i, j, i, j + 1, Az, dz, "z"); sum += G; rhs += G * T[j + 1][i]; } else { const G = hBoundaryEff * Az; sum += G; rhs += G * p.Ta; }
      if (m === 2) { let sliceV = 0; for (let ii = 0; ii < p.Nr; ii++) if (material[j][ii] === 2) sliceV += vols[ii]; const G = gasG[j] * vols[i] / sliceV; sum += G; rhs += G * gasIn[j]; }
      if (m === 3 && i < p.Nr - 1 && material[j][i + 1] === 0) { const Tk = T[j][i] + 273.15, Tak = p.Ta + 273.15, hrad = p.epsTube * p.radArea * sigma * (Tk + Tak) * (Tk * Tk + Tak * Tak), G = hrad * Ae; sum += G; rhs += G * p.Ta; }
      const next = clamp(rhs / Math.max(sum, 1e-30), p.Ta, 2500), updated = clamp(T[j][i] + p.omega * (next - T[j][i]), p.Ta, 2500); maxDelta = Math.max(maxDelta, Math.abs(updated - T[j][i])); T[j][i] = updated;
    }
    if (maxDelta < p.tol) { converged = true; break; }
  }
  let Tmax = -Infinity, Tmin = Infinity, Tavg = 0, bedV = 0;
  for (let j = 0; j < p.Nz; j++) for (let i = 0; i < p.Nr; i++) if (material[j][i] === 2) { Tmax = Math.max(Tmax, T[j][i]); Tmin = Math.min(Tmin, T[j][i]); Tavg += T[j][i] * vols[i]; bedV += vols[i]; } Tavg /= bedV;
  const jc = clamp(Math.floor((p.fbgZ + Hd / 2) / dz), 0, p.Nz - 1), ic = clamp(Math.floor(Math.abs(p.fbgR) / dr), 0, p.Nr - 1), jmid = Math.floor(p.Nz / 2); let iw = clamp(Math.floor(R / dr), 0, p.Nr - 1), is = 0; for (let i = 0; i < p.Nr; i++) { if (material[jmid][i] === 2) is = i; if (material[jmid][i] === 3) iw = i; } const center = T[jmid][0], fbg = T[jc][ic], wall = T[jmid][iw], surface = T[jmid][is];
  let gasTout = p.Ta;
  for (const j of sampleRows) { let sliceV = 0, sliceT = 0; for (let i = 0; i < p.Nr; i++) if (material[j][i] === 2) { sliceV += vols[i]; sliceT += vols[i] * T[j][i]; } sliceT /= sliceV; gasTout += gasG[j] * (sliceT - gasTout) / Math.max(Cg, 1e-30); }
  let qBoundary = 0, qrad = 0;
  for (let j = 0; j < p.Nz; j++) { const i = p.Nr - 1, A = 2 * Math.PI * Rd * dz; qBoundary += hBoundaryEff * A * (T[j][i] - p.Ta); }
  for (let i = 0; i < p.Nr; i++) { const A = areasZ[i]; qBoundary += hBoundaryEff * A * (T[0][i] - p.Ta) + hBoundaryEff * A * (T[p.Nz - 1][i] - p.Ta); }
  for (let j = 0; j < p.Nz; j++) for (let i = 0; i < p.Nr - 1; i++) if (material[j][i] === 3 && material[j][i + 1] === 0) { const A = 2 * Math.PI * (i + 1) * dr * dz, Tk = T[j][i] + 273.15, Tak = p.Ta + 273.15; qrad += p.epsTube * p.radArea * sigma * A * (Math.pow(Tk, 4) - Math.pow(Tak, 4)); }
  let qRadialBed = 0, qAxialBed = 0;
  for (let j = 0; j < p.Nz; j++) for (let i = 0; i < p.Nr; i++) if (material[j][i] === 2) {
    if (i < p.Nr - 1 && material[j][i + 1] !== 2) { const A = 2 * Math.PI * (i + 1) * dr * dz, G = interfaceG(i, j, i + 1, j, A, dr, "r"); qRadialBed += G * (T[j][i] - T[j][i + 1]); }
    if (j > 0 && material[j - 1][i] !== 2) { const G = interfaceG(i, j, i, j - 1, areasZ[i], dz, "z"); qAxialBed += G * (T[j][i] - T[j - 1][i]); }
    if (j < p.Nz - 1 && material[j + 1][i] !== 2) { const G = interfaceG(i, j, i, j + 1, areasZ[i], dz, "z"); qAxialBed += G * (T[j][i] - T[j + 1][i]); }
  }
  const qgas = Cg * (gasTout - p.Ta), balance = p.P - qBoundary - qrad - qgas;
  for (let j = 0; j < p.Nz; j++) for (let i = 0; i < p.Nr; i++) heat[j][i] /= vols[i];
  // Solved on the converged temperatures and reported only: the energy sweep
  // above still uses the legacy gas march, so nothing here feeds back into T.
  // Wiring these face mass flows into the energy equation is Stage 2's job.
  if (solvedField) refreshField();
  const darcy = p.flowMode === "off" ? null : darcyField({ p, T, material, dr, dz, areasZ });
  return { p, T, heat, material, R, Ro, Rd, Hd, dr, dz, V, center, fbg, wall, surface, Tavg, Tmax, Tmin, Tout: gasTout, dpCenter: penetrationDepth(center, p), epsCenter: dielectric(center, p), qgas, qRadialBed, qAxialBed, qBoundary, qrad, balance, it, converged, maxDelta, field: fieldSolve, transport: lastTransport, gasEffectivenessUsed, hBoundaryEff, darcy };
}

export function transportNumbers(sol) {
  const p = sol.p, Rg = 8.314462618, g = 9.80665, TgK = (p.Ta + sol.Tout) / 2 + 273.15, Mg = p.gas === "He" ? .004002602 : .0280134, Cpm = p.gas === "He" ? 20.786 : 29.124, mu0 = p.gas === "He" ? 1.96e-5 : 1.76e-5, S = p.gas === "He" ? 79.4 : 111;
  const mu = mu0 * Math.pow(TgK / 293.15, 1.5) * (293.15 + S) / (TgK + S), rho = p.pressure * Mg / (Rg * TgK), cp = Cpm / Mg, kg = (p.gas === "He" ? .1513 : .0258) * Math.pow(TgK / 293.15, .72), alpha = kg / (rho * cp), Pr = mu * cp / kg, ndot = p.flow * 1e-6 / 60 / .022414, Qactual = ndot * Rg * TgK / p.pressure, area = Math.PI * sol.R * sol.R, u = Qactual / Math.max(area, 1e-30);
  const packed = sol.transport || packedBedTransport((p.Ta + sol.Tavg) / 2, p), rhoBulk = p.rhoBulk, voidFraction = p.voidFraction, ReD = rho * u * p.D / mu, Rep = packed.ReP, ReInterstitial = Rep / voidFraction, PeD = ReD * Pr, Pep = Rep * Pr, Gz = ReD * Pr * p.D / p.H;
  const NuD = 3.66 + .0668 * Gz / (1 + .04 * Math.pow(Math.max(Gz, 0), 2 / 3)), hTube = NuD * kg / p.D, NuP = 2 + 1.1 * Math.pow(Math.max(Rep, 0), .6) * Math.pow(Pr, 1 / 3), hParticle = NuP * kg / p.dp;
  const permeability = packed.permeability, Da = permeability / (sol.R * sol.R), dP = packed.dP, Eu = u > 0 ? dP / (rho * u * u) : 0, Ar = g * Math.pow(p.dp, 3) * rho * Math.max(0, p.rhoSolid - rho) / (mu * mu);
  const TaK = (p.Ta + 273.15), TairK = (p.Ta + sol.wall) / 2 + 273.15, Mair = .0289652, muAir = 1.81e-5 * Math.pow(TairK / 293.15, 1.5) * (293.15 + 111) / (TairK + 111), rhoAir = p.pressure * Mair / (Rg * TairK), cpAir = 1007, kAir = .0263 * Math.pow(TairK / 293.15, .76), PrAir = muAir * cpAir / kAir, nuAir = muAir / rhoAir, alphaAir = kAir / (rhoAir * cpAir), beta = 1 / TairK, Lc = p.domainHeight, dT = Math.abs(sol.wall - p.Ta), Gr = g * beta * dT * Math.pow(Lc, 3) / (nuAir * nuAir), Ra = Gr * PrAir, NuAir = Math.pow(.825 + .387 * Math.pow(Math.max(Ra, 0), 1 / 6) / Math.pow(1 + Math.pow(.492 / PrAir, 9 / 16), 8 / 27), 2), hNatural = NuAir * kAir / Lc;
  const kb = kBed(sol.Tavg, p), TwK = sol.wall + 273.15, hrad = p.epsTube * p.radArea * sigma * (TwK + TaK) * (TwK * TwK + TaK * TaK), BiContact = p.hContact * sol.R / kb, BiRadiation = hrad * sol.R / kb, DpR = sol.dpCenter / sol.R;
  return { TgK, rho, mu, cp, kg, alpha, Pr, Qactual, u, rhoBulk, voidFraction, ReD, Rep, ReInterstitial, PeD, Pep, Gz, NuD, hTube, NuP: packed.NuP, hParticle: packed.hgs, permeability, Da, dP, Eu, Ar, TairK, PrAir, Gr, Ra, NuAir, hNatural, kb, hrad, BiContact, BiRadiation, DpR, specificArea: packed.specificArea, UA: packed.UA, gasEffectiveness: sol.gasEffectivenessUsed, rhoCpEff: packed.rhoCpEff, outletPressure: packed.outletPressure, hBoundaryUsed: sol.hBoundaryEff };
}

// ---- Frequency-domain field -------------------------------------------------
//
// The shipped source is a fitted shape: a Gaussian in (r,z) whose widths are two
// free inputs, multiplied by a Beer-Lambert skin measured inward from the bed
// surface, then renormalised so the total equals the absorbed power. At this
// geometry that is not merely approximate, it has the sign wrong. At 2.404 GHz
// the free-space wavelength is 124.7 mm and SiC at eps' = 7.96 brings it to
// 44.2 mm inside the bed, so a 10 mm bed spans D/lambda = 0.23 and |kR| = 0.71.
// The penetration depth is 140 mm against a 5 mm radius, so the skin term varies
// the source by 3.6% -- edge hot -- while the refraction it omits concentrates
// the field on axis by 14%, which is 30% in power density, centre hot.
//
// Solving for the field instead removes both fitted widths radially and lets
// eps''(T) feed back through the field rather than reweighting a fixed shape.
// The equation is the scalar Helmholtz problem for an axial E in axisymmetric
// coordinates,
//
//     (1/r) d/dr ( r dE/dr ) + d2E/dz2 + k0^2 eps(T,r,z) E = 0
//
// discretised by the same finite-volume pattern as the temperature field, on the
// same mesh. Two things make this cheap rather than a project of its own: the
// load is a quarter wavelength across, so it supports no internal cavity mode
// and nothing outside it has to be meshed; and the existing 0.5 mm cells already
// give 88 per wavelength.
//
// Scope, stated rather than implied. Scalar Helmholtz for E_z is exact when eps
// varies with r alone; axial variation adds a grad(eps) coupling that is dropped
// here, which is defensible only because the load is subwavelength. The incident
// field is taken uniform for the same reason and imposed as E = 1 on the domain
// boundary, so the solution is a shape and the absolute coupling efficiency is
// still not predicted -- the total is renormalised to the absorbed power exactly
// as before. Predicting how much power couples in needs the applicator.
export function permittivityAt(code, tempC, p) {
  if (code === 2) { const e = dielectric(tempC, p); return { re: e.ep, im: -e.epp }; }
  if (code === 3) return { re: p.epsQuartz ?? 3.8, im: -(p.epsQuartzLoss ?? 1e-4) };
  return { re: 1, im: 0 };
}

// BiCGSTAB on the real 2N block form of the complex system. The imaginary part
// lives only on the diagonal -- face coefficients are real -- so the block
// Jacobi preconditioner is an exact 2x2 inverse per cell.
export function solveField2D({ p, T, mesh, maxIterations = 2000, tolerance = 1e-11 }) {
  const { Nr, Nz } = { Nr: p.Nr, Nz: p.Nz }, { dr, dz, material, vols, areasZ } = mesh;
  const n = Nr * Nz, idx = (i, j) => j * Nr + i;
  const k0 = 2 * Math.PI * p.frequency / c0, k02 = k0 * k0;
  const ar = new Float64Array(n), ai = new Float64Array(n);
  const br = new Float64Array(n), bi = new Float64Array(n);
  const gW = new Float64Array(n), gE = new Float64Array(n), gS = new Float64Array(n), gN = new Float64Array(n);
  for (let j = 0; j < Nz; j++) for (let i = 0; i < Nr; i++) {
    const q = idx(i, j), code = material[j][i], eps = permittivityAt(code, T[j][i], p);
    // Inner radial face has zero area on the axis, which is the symmetry
    // condition -- no special case needed.
    const w = 2 * Math.PI * (i * dr) * dz / dr, e = 2 * Math.PI * ((i + 1) * dr) * dz / dr;
    const s = areasZ[i] / dz, nn = areasZ[i] / dz;
    gW[q] = w; gE[q] = e; gS[q] = s; gN[q] = nn;
    const sum = w + e + s + nn;
    ar[q] = sum - k02 * eps.re * vols[i];
    ai[q] = -k02 * eps.im * vols[i];
    // Uniform incident field on every outer face.
    if (i === Nr - 1) { br[q] += e; }
    if (j === 0) { br[q] += s; }
    if (j === Nz - 1) { br[q] += nn; }
  }
  const mul = (xr, xi, yr, yi) => {
    for (let j = 0; j < Nz; j++) for (let i = 0; i < Nr; i++) {
      const q = idx(i, j);
      let sr = ar[q] * xr[q] - ai[q] * xi[q], si = ai[q] * xr[q] + ar[q] * xi[q];
      if (i > 0) { const m = idx(i - 1, j); sr -= gW[q] * xr[m]; si -= gW[q] * xi[m]; }
      if (i < Nr - 1) { const m = idx(i + 1, j); sr -= gE[q] * xr[m]; si -= gE[q] * xi[m]; }
      if (j > 0) { const m = idx(i, j - 1); sr -= gS[q] * xr[m]; si -= gS[q] * xi[m]; }
      if (j < Nz - 1) { const m = idx(i, j + 1); sr -= gN[q] * xr[m]; si -= gN[q] * xi[m]; }
      yr[q] = sr; yi[q] = si;
    }
  };
  const precondition = (vr, vi, or_, oi) => {
    for (let q = 0; q < n; q++) {
      const d = ar[q] * ar[q] + ai[q] * ai[q] || 1e-30;
      or_[q] = (ar[q] * vr[q] + ai[q] * vi[q]) / d;
      oi[q] = (ar[q] * vi[q] - ai[q] * vr[q]) / d;
    }
  };
  const alloc = () => [new Float64Array(n), new Float64Array(n)];
  const [xr, xi] = alloc(), [rr, ri] = alloc(), [hr, hi] = alloc();
  const [pr, pi] = alloc(), [vr, vi] = alloc(), [sr_, si_] = alloc();
  const [tr, ti] = alloc(), [phr, phi] = alloc(), [shr, shi] = alloc();
  for (let q = 0; q < n; q++) xr[q] = 1;
  mul(xr, xi, vr, vi);
  let bNorm = 0, rNorm = 0;
  for (let q = 0; q < n; q++) {
    rr[q] = br[q] - vr[q]; ri[q] = bi[q] - vi[q]; hr[q] = rr[q]; hi[q] = ri[q];
    bNorm += br[q] * br[q] + bi[q] * bi[q]; rNorm += rr[q] * rr[q] + ri[q] * ri[q];
  }
  bNorm = Math.sqrt(Math.max(bNorm, 1e-30));
  const dot = (a1, a2, b1, b2) => { let s = 0; for (let q = 0; q < n; q++) s += a1[q] * b1[q] + a2[q] * b2[q]; return s; };
  let relative = Math.sqrt(rNorm) / bNorm, rhoOld = 1, alpha = 1, omega = 1, iteration = 0;
  for (; iteration < maxIterations && relative > tolerance; iteration++) {
    const rho = dot(hr, hi, rr, ri);
    if (Math.abs(rho) < 1e-300) break;
    const beta = (rho / rhoOld) * (alpha / (omega || 1e-30));
    for (let q = 0; q < n; q++) {
      pr[q] = rr[q] + beta * (pr[q] - omega * vr[q]);
      pi[q] = ri[q] + beta * (pi[q] - omega * vi[q]);
    }
    precondition(pr, pi, phr, phi); mul(phr, phi, vr, vi);
    const denominator = dot(hr, hi, vr, vi);
    if (Math.abs(denominator) < 1e-300) break;
    alpha = rho / denominator;
    for (let q = 0; q < n; q++) { sr_[q] = rr[q] - alpha * vr[q]; si_[q] = ri[q] - alpha * vi[q]; }
    precondition(sr_, si_, shr, shi); mul(shr, shi, tr, ti);
    const tt = dot(tr, ti, tr, ti);
    omega = tt > 1e-300 ? dot(tr, ti, sr_, si_) / tt : 0;
    let next = 0;
    for (let q = 0; q < n; q++) {
      xr[q] += alpha * phr[q] + omega * shr[q];
      xi[q] += alpha * phi[q] + omega * shi[q];
      rr[q] = sr_[q] - omega * tr[q]; ri[q] = si_[q] - omega * ti[q];
      next += rr[q] * rr[q] + ri[q] * ri[q];
    }
    relative = Math.sqrt(next) / bNorm;
    if (Math.abs(omega) < 1e-300) break;
    rhoOld = rho;
  }
  // Dissipation density follows |E|^2 weighted by the local loss, which is the
  // only place eps'' enters the heating. Normalisation to the absorbed power is
  // left to the caller, exactly as the fitted shape was.
  const magnitude = new Float64Array(n), heat = new Float64Array(n);
  for (let j = 0; j < Nz; j++) for (let i = 0; i < Nr; i++) {
    const q = idx(i, j);
    magnitude[q] = xr[q] * xr[q] + xi[q] * xi[q];
    heat[q] = material[j][i] === 2 ? -permittivityAt(2, T[j][i], p).im * magnitude[q] : 0;
  }
  return { Er: xr, Ei: xi, magnitude, heat, iterations: iteration, relativeResidual: relative, k0 };
}
