# Verification record, 2026-09-08

Scope: the new `microwave3d/` implementation. Verification checks numerical
equations and conservation; it does not establish agreement with the supplied
reactor experiments. No measured-temperature curve was fitted to obtain these
results.

## Automated checks

`python -m pytest -q` exercises the following independent conditions:

| Check | Acceptance |
| --- | --- |
| Rectangular PEC TE101 vector mode | Discrete eigen-equation residual < 10⁻¹²; continuum-frequency error decreases at second order |
| 3D sinusoidal manufactured temperature, 6³/12³/24³ | Error reduction > 3.8 per grid halving; boundary heat closes the imposed source |
| Lossless closed cavity | All incident power reflected |
| Lossy cavity with port and metal wall | Passive reflection and electromagnetic power closure |
| Dielectric, ohmic, and aperture losses | All explicitly included in power balance |
| Power scaling at fixed properties | Quadrupling power doubles electric amplitude and quadruples deposition |
| Cell heat quadrature | Integrated heat equals the Maxwell dielectric mass loss |
| Incoherent frequency dwell | Weights normalize average power, not phasor amplitudes |
| Isothermal Darcy bed | Analytical pressure drop and per-cell mass conservation |
| Advection with heating | Outlet enthalpy equals volumetric source; independent streamlines do not mix artificially |
| Gray enclosure radiation | Matches analytical two-surface radiation resistance |
| Coupled calculation | Final temperatures pass a newly evaluated EM/flow/thermal audit |
| Zero drive | Ambient field; S11 remains a scattering property, not a fabricated zero-reflection claim |
| Failed calculation / invalid inputs | Explicit failure, no success flag or property extrapolation |
| Export | Numeric arrays, VTK vectors, self-contained viewer, and failed-run status retained |

The unchanged JavaScript regression suite also passes: 25 tests. Those tests
preserve the old model's behavior; they do not remedy its previously documented
physical or grid-convergence limitations.

## Empty cylindrical cavity

Command: `python tools/verification/microwave3d.py --grids 8,12,16,24`.
For a closed PEC cylinder of radius 52.46 mm and height 85 mm, the TE111
analytical frequency is 2.431901852 GHz, using the first zero of J₁′,
1.8411837813406593. The vector eigenproblem has the expected doublet.

| Uniform cells per axis | Electric unknowns | Calculated GHz | Frequency error | Cavity volume error |
| --- | ---: | ---: | ---: | ---: |
| 8 | 912 | 2.360775820 | −2.925% | +3.451% |
| 12 | 3268 | 2.402433228 | −1.212% | −0.970% |
| 16 | 8592 | 2.393122973 | −1.595% | +3.451% |
| 24 | 29128 | 2.423060087 | −0.364% | −0.970% |

The stepped circumference produces nonmonotone geometric error. A 0.36%
frequency error is about 8.8 MHz, which can exceed a narrow experimental sweep
bandwidth. This benchmark therefore **does not certify the default loaded
cavity's absorption or temperature**. Curved-boundary accuracy and frequency
resolution remain essential before experimental prediction.

## Illustrative coupled power sweep

Command: `python -m microwave3d.cli --powers 0,30,60,90,130 --output results3d/power-sweep`.
The default 14 × 14 × 14 clustered-grid case uses an assumed 8 mm lumped probe,
an illustrative SiC-bed table, and 2.404 GHz. It is deliberately not tuned to the
paper's experimental temperatures.

| Incident W | Sample absorption W | Centre °C | Quartz midplane °C |
| ---: | ---: | ---: | ---: |
| 0 | 0 | 20.000 | 20.000 |
| 30 | 0.060469 | 25.127 | 25.063 |
| 60 | 0.121209 | 30.186 | 30.057 |
| 90 | 0.182214 | 35.177 | 34.983 |
| 130 | 0.263956 | 41.730 | 41.448 |

All five states satisfy the configured convergence tests. These low absorptions
show why the assumed port, coarse geometry, and chosen frequency cannot stand in
for the tuned experiment. They are not evidence that the actual reactor fails
to heat. Do not use the example curve to infer experimental centre temperature.

## Remaining validation work

The channel extension adds an analytical flow benchmark and loaded channel-thermal
mesh sensitivity, recorded below. No loaded **EM** grid convergence, frequency-quadrature
convergence, experimental S11 fit, thermal-camera/FBG observation model, hot
variable-density flow validation, or uncertainty intervals are claimed.
The paper provides comparison targets, not automatic validation.

## Resolved channel extension

The Python suite now has 22 tests; the unchanged JavaScript suite has 25. All pass.
The added tests verify square-duct pressure drop, no-slip walls, weak continuity,
nonlinear momentum residual, phase-conservative intermesh heat transfer,
coupled channel enthalpy and power closure, invalid geometry rejection, and
agreement between iterative and direct thermal linear solves.

For a 1.3 mm square channel, L = 15 mm, U = 0.05 m/s, μ = 1.76 × 10⁻⁵ Pa s,
the Fourier-series pressure drop is 0.2222454614 Pa. With fully developed
velocities prescribed at both ends and six axial cells:

| Cross-section subdivisions | Centreline end-pressure difference error |
| ---: | ---: |
| 2 | +4.86786% |
| 3 | +1.17123% |
| 4 | +0.40507% |

The reported production pressure drop instead uses FEM area integrals over
the exact inlet/outlet planes. It does not use first/last cell centres, which
would omit half-cell lengths and introduce an artificial axial-grid dependence.

The 21-channel example at 130 W and 2.404 GHz converges with Re = 4.3353,
boundary-to-boundary pressure drop 0.2292874 Pa, nonlinear NS residual
3.26 × 10⁻¹¹, and conservative thermal-cell mass residual below 10⁻¹⁴.
The coarse thermal transfer corrects the FEM face flux by 0.5854% in relative L2
norm. This correction is explicitly reported and is not identified with the
FEM continuity residual.

This illustrative run absorbs only about 0.98 W in the sample because of the
assumed untuned port. It is not an experimental temperature prediction. At
the calculated gas temperatures, the ideal-gas density differs from the fixed
flow-reference density by approximately 18%; the output flags this limitation.

Targeted browser verification with agent-browser confirmed channel-only display,
Navier–Stokes speed selection, and XY cross sections without browser errors.

### Loaded channel-thermal refinement

At the same 130 W, 2.404 GHz setting, the EM mesh stays 14 × 14 × 14 and the
representative FEM flow stays 3 × 3 × 12 subdivisions. Thermal channel width
and length are refined together. The first series used
`python -m microwave3d.convergence --kind thermal --levels 3,4,5`;
level 6 additionally uses `thermal_cells_across: 6` and
`thermal_cells_axial: 36` in the same input case. The sequence can be reproduced
with `--levels 3,4,5,6`.

| Channel subdivisions across / axial | Full thermal grid | Sample maximum °C | Change from preceding level K |
| --- | --- | ---: | ---: |
| 3 / 18 | 34 × 34 × 26 | 82.2505 | — |
| 4 / 24 | 38 × 38 × 32 | 80.9962 | 1.2543 |
| 5 / 30 | 44 × 44 × 38 | 79.9353 | 1.0609 |
| 6 / 36 | 48 × 48 × 44 | 79.0076 | 0.9277 |

The final pair meets the 1 K temperature-change threshold; the preceding two
pairs do not. This is a pairwise sensitivity result, **not** a 1 K error bound.
The coarse cavity geometry, homogenized EM interfaces, radiation approximation,
and fixed-density flow are unchanged. Their uncertainties are not reduced by
this thermal refinement. These temperatures are illustrative model outputs.
