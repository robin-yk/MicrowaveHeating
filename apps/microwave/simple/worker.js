// Solver worker for the single-screen page. Builds the parameter object from a
// material profile's defaults (same derivation as ../index.html parameters())
// and runs solve2D off the main thread so dragging stays responsive.
import { materialProfiles, parseRows, clamp, solve2D } from "../solver.js";

function parameters(id, P, coarse) {
  const profile = materialProfiles[id], d = profile.defaults, D = d.diameter / 1000, H = d.length / 1000;
  const p = {
    materialId: id, materialLabel: profile.label, materialFormula: profile.formula, rhoSolid: profile.rhoSolid, cpSolid: profile.cpSolid,
    P: Math.max(0, P), frequency: d.frequency * 1e9, volume: Math.PI * D * D * H / 4 * 1e6, D, H, mass: d.mass, gas: d.gas, flow: d.flow,
    pressure: d["gas-pressure"] * 1e5, dp: d["particle-diameter"] * 1e-6, Ta: d.ambient,
    Nr: coarse ? 15 : 30, Nz: coarse ? 30 : 60, domainWidth: .03, domainHeight: .03,
    k200: d.k200, k500: d.k500, k800: d.k800, kzRatio: d["kz-ratio"], hContact: d["h-contact"], kq: d["k-quartz"], tq: d["tube-thickness"] / 1000,
    airFactor: d["air-factor"], boundaryMode: d["boundary-mode"], hBoundary: d["h-boundary"], epsTube: d.emissivity, radArea: d["rad-area"],
    gasTransferMode: d["gas-transfer-mode"], gasEff: d["gas-eff"], dielectricMode: d["dielectric-mode"], bedKMode: d["bed-k-mode"],
    fieldWr: d["field-wr"], fieldWz: d["field-wz"], fieldMode: profile.fieldMode, fbgR: d["fbg-r"] / 1000, fbgZ: d["fbg-z"] / 1000,
    diel: parseRows(profile.dielectric, 3), maxIter: coarse ? 3500 : 6000, tol: coarse ? 1.5e-3 : 3e-4, omega: coarse ? 1.08 : 1.05
  };
  const rhoBulk = (p.mass / 1000) / (p.volume * 1e-6), refMass = d.mass / 1000, refVolume = d.volume * 1e-6;
  p.rhoBulk = rhoBulk; p.voidFraction = clamp(1 - rhoBulk / p.rhoSolid, .01, .99); p.referenceVoidFraction = clamp(1 - (refMass / refVolume) / p.rhoSolid, .01, .99);
  return p;
}

function summary(s) {
  return { fbg: s.fbg, wall: s.wall, avg: s.Tavg, max: s.Tmax, converged: s.converged };
}

self.onmessage = event => {
  const m = event.data;
  if (m.type === "solve") {
    const s = solve2D(parameters(m.material, m.P, m.coarse));
    const T = Float32Array.from(s.T.flat());
    self.postMessage({ type: "solve", id: m.id, P: m.P, material: m.material, coarse: m.coarse, Nr: s.p.Nr, Nz: s.p.Nz, R: s.R, Ro: s.Ro, Rd: s.Rd, Hd: s.Hd, H: s.p.H, Ta: s.p.Ta, T, ...summary(s) }, [T.buffer]);
  } else if (m.type === "curve") {
    const points = m.powers.map(P => ({ P, ...summary(solve2D(parameters(m.material, P, m.coarse))) }));
    self.postMessage({ type: "curve", id: m.id, material: m.material, coarse: m.coarse, points });
  }
};
