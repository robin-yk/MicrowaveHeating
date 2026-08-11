// Microwave heating numeric core: dielectric mixing, bed conductivity, gas
// transport, and the 2D axisymmetric FVM steady-state thermal solver behind
// microwave.html. Pure functions of plain input objects — no DOM access —
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
export function packedBedTransport(tempC, p) {
  const g = 9.80665, eps = p.voidFraction, solid = Math.max(1e-8, 1 - eps), Cg = gasCapacityRate(p); let pressure = p.pressure, state, dP = 0, ReP = 0, NuP = 2, hgs = 0, permeability = 0;
  for (let pass = 0; pass < 3; pass++) {
    state = gasState(tempC, p, pressure); ReP = state.rho * state.u * p.dp / state.mu; NuP = 2 + 1.1 * Math.pow(Math.max(ReP, 0), .6) * Math.pow(state.Pr, 1 / 3); hgs = NuP * state.kg / p.dp;
    permeability = Math.pow(eps, 3) * p.dp * p.dp / (180 * solid * solid);
    const dPdz = 150 * state.mu * solid * solid * state.u / (Math.pow(eps, 3) * p.dp * p.dp) + 1.75 * state.rho * solid * state.u * state.u / (Math.pow(eps, 3) * p.dp);
    dP = dPdz * p.H; pressure = Math.max(5000, p.pressure - dP / 2);
  }
  const specificArea = 6 * solid / p.dp, UA = hgs * specificArea * p.volume * 1e-6, gasEffectiveness = Cg > 0 ? 1 - Math.exp(-UA / Cg) : 0, rhoCpEff = solid * p.rhoSolid * p.cpSolid + eps * state.rho * state.cp;
  return { ...state, Cg, ReP, NuP, hgs, specificArea, UA, gasEffectiveness, permeability, dP, outletPressure: Math.max(5000, p.pressure - dP), rhoCpEff };
}
export function naturalConvection(tempC, p) {
  const Rg = 8.314462618, g = 9.80665, Tfilm = (tempC + p.Ta) / 2 + 273.15, Mair = .0289652, mu = 1.81e-5 * Math.pow(Tfilm / 293.15, 1.5) * (293.15 + 111) / (Tfilm + 111), rho = p.pressure * Mair / (Rg * Tfilm), cp = 1007, k = .0263 * Math.pow(Tfilm / 293.15, .76), Pr = mu * cp / k, nu = mu / rho, beta = 1 / Tfilm, L = Math.max(p.domainHeight, 1e-6), dT = Math.abs(tempC - p.Ta), Gr = g * beta * dT * Math.pow(L, 3) / (nu * nu), Ra = Gr * Pr, Nu = Math.pow(.825 + .387 * Math.pow(Math.max(Ra, 0), 1 / 6) / Math.pow(1 + Math.pow(.492 / Pr, 9 / 16), 8 / 27), 2);
  return { Tfilm, Pr, Gr, Ra, Nu, h: Nu * k / L };
}

export function solve2D(p) {
  const R = p.D / 2, Ro = R + p.tq, Rd = p.domainWidth / 2, Hd = p.domainHeight, dr = Rd / p.Nr, dz = Hd / p.Nz;
  const material = Array.from({ length: p.Nz }, (_, j) => Array.from({ length: p.Nr }, (_, i) => {
    const r = (i + .5) * dr, z = (j + .5) * dz - Hd / 2;
    if (r < R && Math.abs(z) < p.H / 2) return 2;
    if (r < R) return 1;
    if (r < Ro) return 3;
    return 0;
  }));
  const T = Array.from({ length: p.Nz }, (_, j) => Array.from({ length: p.Nr }, (_, i) => {
    const r = (i + .5) * dr, z = (j + .5) * dz - Hd / 2, d = Math.hypot(Math.max(0, r - R), Math.max(0, Math.abs(z) - p.H / 2));
    return p.Ta + p.P * 22 * Math.exp(-d / .006);
  }));
  const vols = [], areasZ = []; let V = 0;
  for (let i = 0; i < p.Nr; i++) { const rw = i * dr, re = (i + 1) * dr; vols[i] = Math.PI * (re * re - rw * rw) * dz; areasZ[i] = Math.PI * (re * re - rw * rw); }
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
  for (it = 0; it < p.maxIter; it++) {
    let bedTemp = 0, bedTempV = 0; for (let j = 0; j < p.Nz; j++) for (let i = 0; i < p.Nr; i++) if (material[j][i] === 2) { bedTemp += T[j][i] * vols[i]; bedTempV += vols[i]; } bedTemp /= Math.max(bedTempV, 1e-30);
    let wallIndex = clamp(Math.floor(R / dr), 0, p.Nr - 1); for (let i = 0; i < p.Nr; i++) if (material[Math.floor(p.Nz / 2)][i] === 3) wallIndex = i; const wallGuess = T[Math.floor(p.Nz / 2)][wallIndex];
    lastTransport = packedBedTransport((bedTemp + p.Ta) / 2, p); gasEffectivenessUsed = p.gasTransferMode === "manual" ? p.gasEff : lastTransport.gasEffectiveness;
    const UAg = p.gasTransferMode === "manual" ? (Cg > 0 ? -Math.log(Math.max(1e-8, 1 - gasEffectivenessUsed)) * Cg : 0) : lastTransport.UA;
    hBoundaryEff = p.boundaryMode === "manual" ? p.hBoundary : naturalConvection(wallGuess, p).h;
    let norm = 0;
    for (let j = 0; j < p.Nz; j++) for (let i = 0; i < p.Nr; i++) {
      if (material[j][i] !== 2) { heat[j][i] = 0; continue; }
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
  return { p, T, heat, material, R, Ro, Rd, Hd, dr, dz, V, center, fbg, wall, surface, Tavg, Tmax, Tmin, Tout: gasTout, dpCenter: penetrationDepth(center, p), epsCenter: dielectric(center, p), qgas, qRadialBed, qAxialBed, qBoundary, qrad, balance, it, converged, maxDelta, transport: lastTransport, gasEffectivenessUsed, hBoundaryEff };
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
