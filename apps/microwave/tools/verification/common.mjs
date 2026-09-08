// Shared helpers for the numerical-verification runners. Everything here is
// about measuring discretization error: norms over cell-centered fields,
// observed convergence order between grid pairs, and Richardson extrapolation
// for cases with no exact solution.
"use strict";

// Volume-weighted L2 and pointwise Linf of (field - exact) over cells selected
// by `include(i, j)`. `field[j][i]`, `exact(i, j)`, `volume(i, j)`.
export function errorNorms({ nr, nz, field, exact, volume, include = () => true }) {
  let sum2 = 0, vol = 0, linf = 0;
  for (let j = 0; j < nz; j++) for (let i = 0; i < nr; i++) {
    if (!include(i, j)) continue;
    const e = field[j][i] - exact(i, j), v = volume(i, j);
    sum2 += e * e * v; vol += v; linf = Math.max(linf, Math.abs(e));
  }
  return { l2: Math.sqrt(sum2 / vol), linf };
}

// Observed order between successive grids with mesh ratio 2.
export const observedOrder = (coarse, fine) => Math.log2(coarse / fine);

// Richardson extrapolation for a scalar q computed on grids h, h/2, h/4.
// Returns the observed order, the extrapolated value, and the GCI-style
// relative error of the finest grid against that extrapolated value.
export function richardson(q1, q2, q3) { // q1 coarsest ... q3 finest
  const r = (q2 - q1) / (q3 - q2);
  const p = Math.log2(Math.abs(r));
  const qExtrap = q3 + (q3 - q2) / (Math.pow(2, p) - 1);
  return { p, qExtrap, fineError: Math.abs((q3 - qExtrap) / qExtrap) };
}

export function markdownTable(headers, rows) {
  const line = (cells) => "| " + cells.join(" | ") + " |";
  return [line(headers), line(headers.map(() => "---")), ...rows.map(line)].join("\n");
}

export const sci = (x, digits = 3) => Number(x).toExponential(digits);
export const fix = (x, digits = 3) => Number(x).toFixed(digits);
