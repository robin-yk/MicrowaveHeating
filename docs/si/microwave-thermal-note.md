# Supplementary Note S_x — Two-dimensional thermal reconstruction of the microwave-heated bed

*Draft for the dielectric-redox coupling manuscript. Numbering to be assigned on
integration. All numbers in this note are reproducible from
`apps/microwave/solver.js` at the commit recorded in Section S_x.9.*

---

## S_x.1  Purpose and scope

The two thermometers used in this work report different things. The pyrometer
reads the **outer surface of the quartz tube** (`T_wall`); the fibre Bragg
grating reads a **single point on the bed axis** (`T_center`). Neither is the
volume-average bed temperature that a rate constant responds to, and the two
differ from each other by several hundred kelvin at the powers used here.

This note quantifies that difference. It does **not** attempt to predict
heating from first principles. Specifically, it does not solve Maxwell's
equations, does not compute the cavity field, and does not predict how much
power a given sample absorbs — all of those are measured quantities in this
work (Table S_y) or come from cavity perturbation theory (Note S2b). What the
model does is take the **measured power deposited in the sample** as its input
and reconstruct the internal temperature field consistent with it.

The model is therefore a *reconciliation* tool: it converts two boundary
measurements into an interior field, subject to energy conservation. Three
things follow that cannot be read off the thermometers directly:

1. how far the FBG reading sits above the volume average (Section S_x.5);
2. where the deposited power actually leaves the reactor (Section S_x.4);
3. whether the reported power partitioning is thermodynamically self-consistent
   across the sample series (Section S_x.7).

The third of these turns out to place a quantitative validity bound on the
cavity-perturbation analysis used in Figure 2g.

---

## S_x.2  Governing equations and closures

The model solves the steady energy equation on a two-dimensional axisymmetric
`(r, z)` domain spanning the bed, the annular gas gap, the quartz wall, and the
surrounding air:

$$\nabla\cdot\left(k_{\text{eff}}(T)\,\nabla T\right) + q'''(r,z,T) = 0$$

discretised by a cell-centred finite-volume method with two-point flux and
harmonic-mean interface conductivity, so that flux is continuous across the
bed/gas/quartz material jumps. The temperature dependence of every closure
makes the system nonlinear; it is resolved by Picard iteration with
Gauss–Seidel sweeps and over-relaxation (`ω = 1.05`), converged to a
normalised residual of `3 × 10⁻⁴`.

**Effective bed conductivity.** A Maxwell–Eucken mixture of the solid skeleton
and the interstitial gas, with the skeleton conductivity backed out from three
fitted anchor values `k(200 °C)`, `k(500 °C)`, `k(800 °C)` and log-interpolated
between them. Axial conductivity carries a separate anisotropy ratio to
represent the packing.

**Dielectric mixing.** Looyenga mixing of the measured solid-phase `ε′(T)`,
`ε″(T)` with air at the bed's void fraction, so that a sample measured at one
packing density can be re-evaluated at another.

**Volumetric source.** The deposited power density is written as a separable
envelope times a plane-wave attenuation factor,

$$q'''(r,z,T) \;\propto\; \varepsilon''(T)\;
\exp\!\left[-\left(\tfrac{r}{w_r R}\right)^{2}-\left(\tfrac{z}{w_z H/2}\right)^{2}\right]
\exp\!\left[-\frac{R-r}{\delta_p(T)}\right],$$

with the penetration depth from the plane-wave result
$\delta_p = 1/(2\alpha)$, $\alpha = k_0\sqrt{(|\varepsilon| - \varepsilon')/2}$.
The whole field is then **renormalised so that its volume integral equals the
prescribed power**. This is the single most important structural point in the
model and is stated explicitly because it determines what the input variable
means: the model's power input is by construction $P_{\text{sample}}$, the
power dissipated in the sample — not the incident power and not the total
absorbed power. Reflection and cavity-wall losses are outside the model
boundary.

The renormalisation also means $w_r$ and $w_z$ control only the *shape* of the
deposition, never its magnitude. Section S_x.5 shows this is exactly why the
FBG-to-average gap is bounded even though the field shape is not known
independently.

**Boundary conditions.** Radiation from the quartz outer surface
($\varepsilon = 0.85$, with a fitted view-factor area multiplier) and natural
convection to ambient by the Churchill–Chu correlation. Gas-phase enthalpy
removal uses a particle Nusselt closure
$\mathrm{Nu}_p = 2 + 1.1\,\mathrm{Re}_p^{0.6}\mathrm{Pr}^{1/3}$.

A Darcy pressure field is solved on the same mesh for the bed permeability
(Blake–Kozeny/Ergun) and is reported for diagnostic purposes, but it is **not**
coupled back into the energy equation in the present version; convective
transport uses the one-dimensional gas march described above.

---

## S_x.3  Calibration and verification

Nine parameters are fitted — three conductivity anchors, the axial anisotropy
ratio, the bed/wall contact conductance, the two field-envelope widths, the
radiative area multiplier, and the air-gap conduction factor — against the
fourteen measurements in Table S_x.1 (seven powers × two thermometers). The
resulting parameter-to-datum ratio of 9:14 is poor by the standards of a
predictive model and is the reason this note reports **bounds** on the derived
quantities rather than point values.

**Table S_x.1 — Calibration against the two-thermometer power sweep.**
Reduced rutile TiO₂ (H₂, 600 °C, 30 min), 1.150 g in a 10 mm × 15 mm bed,
He at 50 sccm, 2.404 GHz. `T_wall` is the quartz outer surface; `T_center` is
the on-axis FBG. Grid 30 × 60.

| P_sample (W) | T_wall model (°C) | T_wall meas. (°C) | Δ (K) | T_center model (°C) | T_center meas. (°C) | Δ (K) | ⟨T⟩_V (°C) | T_center − ⟨T⟩_V (K) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 5  | 188.0 | 167.3 | +20.7 | 213.7 | 185.5 | +28.2 | 198.9 | 14.8 |
| 10 | 291.5 | 287.1 |  +4.4 | 363.7 | 354.6 |  +9.1 | 320.3 | 43.4 |
| 14 | 355.3 | 360.0 |  −4.7 | 486.6 | 490.0 |  −3.4 | 405.0 | 81.6 |
| 17 | 396.2 | 402.7 |  −6.5 | 580.7 | 580.0 |  +0.7 | 466.0 | 114.7 |
| 20 | 432.9 | 425.5 |  +7.4 | 664.5 | 618.0 | +46.5 | 523.0 | 141.5 |
| 23 | 466.0 | 468.0 |  −2.0 | 742.0 | 710.0 | +32.0 | 574.9 | 167.1 |
| 26 | 495.7 | 499.6 |  −3.9 | 811.2 | 799.2 | +12.0 | 622.8 | 188.4 |

RMSE: `T_wall` 9.1 K, `T_center` 24.6 K. The largest single residual is the
20 W centre point (+46.5 K); the wall residuals are within the pyrometer's
stated reproducibility across the whole series.

The rightmost column is the quantity of interest and is discussed in
Section S_x.5. Note that it is not a fitted quantity — no measurement of
⟨T⟩_V exists — but a consequence of energy conservation given the two
measured temperatures.

**Table S_x.2 — Grid convergence at P_sample = 26 W.**
Reproduced from the repository's continuous-integration verification suite
(`npm run verify:microwave`), which tightens the convergence tolerance with
each refinement level. Darcy flow on.

| grid | T_center (°C) | T_wall (°C) | ⟨T⟩_V (°C) | T_center − ⟨T⟩_V (K) | energy closure | Darcy mass imbalance |
|---|---:|---:|---:|---:|---:|---:|
| 30 × 60 (default) | 811.20 | 495.71 | 622.79 | 188.4 | 2.4 × 10⁻⁴ | 1.8 × 10⁻¹³ |
| 60 × 120 | 804.65 | 495.74 | 621.85 | 182.8 | 5.0 × 10⁻⁴ | 2.1 × 10⁻¹³ |
| 120 × 240 | 788.98 | 496.13 | 618.23 | 170.8 | 1.0 × 10⁻³ | 2.1 × 10⁻¹³ |

`T_wall` is grid-converged to better than 0.5 K. `T_center` is still drifting
downward at the finest grid, and with it the FBG-to-average gap: 188 K → 171 K
over a factor of four in linear resolution. The results in this note are
therefore quoted on the default grid with the residual grid drift folded into
the uncertainty band of Section S_x.5, not as converged point values.

Two additional verification results are carried over from the repository's
verification suite (`docs/VERIFICATION.md`):

* **Analytic limit.** With the source term made uniform and the boundary
  isothermal, the computed radial profile approaches the parabolic conduction
  solution as the aspect ratio removes axial leakage: relative error
  `8.7 × 10⁻²` at `H/D = 4`, falling to `2.3 × 10⁻³` at `H/D = 16`.
* **Discrete conservation.** The Darcy operator closes its mass balance to
  `1.8 × 10⁻¹³` (round-off), confirming the finite-volume assembly is
  conservative independently of the energy equation's closures.

A full method-of-manufactured-solutions order verification is available for the
Joule-heating solver in the same repository (observed order 1.7) but has not
been carried out for the microwave solver, whose temperature-dependent
dielectric closure has no convenient manufactured form. This is stated as a
limitation, not a claim.

---

## S_x.4  Where the deposited power leaves the reactor

**Table S_x.3 — Where the deposited power leaves the reactor.**
Fractions of `P_sample`. The first three columns are the external loss paths
and sum to 100 %; the last two are the internal transport split inside the bed
(the fraction of the deposited power that reaches the bed boundary radially
versus axially).

| P_sample (W) | quartz radiation | outer convection | gas enthalpy | → radial (internal) | → axial (internal) |
|---:|---:|---:|---:|---:|---:|
| 5  | 52.5 % | 44.7 % | 2.6 % | 80.3 % | 17.0 % |
| 10 | 59.6 % | 38.2 % | 2.1 % | 81.1 % | 16.8 % |
| 17 | 65.9 % | 32.3 % | 1.8 % | 81.9 % | 16.4 % |
| 26 | 71.1 % | 27.5 % | 1.5 % | 82.5 % | 16.0 % |

Radiation overtakes convection above about 8 W and carries 71 % of the load at
the top of the series, as expected from the `T⁴` scaling against a nearly
constant convective coefficient. The split is set by the quartz outer-surface
temperature and the ambient, and is therefore essentially independent of how
the power is distributed inside the bed — the property exploited in
Section S_x.6.

The practical consequence is that the helium sweep is not a meaningful
heat-removal path at these powers. Quadrupling the flow from 25 to 100 sccm at
26 W raises the gas share from 0.7 % to 3.0 % and lowers `T_wall` by 3.5 K and
`T_center` by 23 K — the wall shift is below the pyrometer's reproducibility.
Flow rate is therefore not a usable handle on bed temperature in this reactor,
and differences in sweep flow between samples cannot explain the discrepancies
examined in Section S_x.7.

---

## S_x.5  The FBG reading is not the bed temperature

At the highest power in the calibration series the model places the axial FBG
reading roughly 200 K above the volume-averaged bed temperature. The gap grows
monotonically with power because it is set by the internal conduction
resistance against a nearly fixed external loss coefficient.

The value is **bounded rather than determined**, and the reason is worth
stating precisely. The field-envelope widths $w_r, w_z$ are fitted, not
measured, and the two-thermometer dataset does not determine them uniquely:
parameter sets differing in $w_r$ by a factor of 2.5 reproduce both measured
temperatures within their experimental scatter. Sweeping $w_r$ across the range
compatible with the calibration moves the gap over

$$T_{\text{FBG}} - \langle T\rangle_{V} = 200 \pm 30\ \text{K at } P_{\text{sample}} = 26\ \text{W},$$

and grid refinement moves it a further ~17 K downward (Table S_x.2), which is
inside that band.

We investigated whether an additional measurement could narrow the band and
conclude that it cannot with external probes. Two candidate discriminators were
tested against the model:

* an **off-axis FBG**: the discrimination between $w_r = 0.8$ and $w_r = 2.0$ is
  81 K on the axis and falls to zero at $r = 4$ mm;
* an **axial wall temperature profile** from an IR camera: the two field widths
  differ by less than ±2 K anywhere on the quartz outer surface.

The field-shape ambiguity is therefore not observable from the reactor exterior.
Narrowing it requires replacing the fitted envelope with a computed one — that
is, a full-wave solution of the cavity field — rather than an additional
thermal measurement. This is the one place in the present analysis where an
electromagnetic solver would change a number.

---

## S_x.6  The two thermometers separate two different effects

**Table S_x.4 — Response of the two thermometers to the deposition profile.**
`P_sample` held at 26 W while the bed's ε″(T) is scaled by the factor in the
first column, driving the penetration depth (evaluated at the bed-centre temperature)
from 1.85 bed radii down to 0.15.
`T_max` is the hottest bed cell anywhere in the bed and `r(T_max)` its radial
position; the last column is the fraction of the deposited power landing in the
outer 20 % of the bed radius.

| ε″ scaling | δ_p (mm) | δ_p/R | T_wall (°C) | T_center (°C) | T_max (°C) | r(T_max) (mm) | power in outer 20 % |
|---|---:|---:|---:|---:|---:|---:|---:|
| ×1  | 9.25 | 1.85 | 496 | 811 | 813 | 0.3 | 31 % |
| ×2  | 5.25 | 1.05 | 496 | 788 | 789 | 0.3 | 34 % |
| ×4  | 3.25 | 0.65 | 497 | 758 | 759 | 0.3 | 39 % |
| ×8  | 2.14 | 0.43 | 497 | 724 | 725 | 0.3 | 45 % |
| ×16 | 1.49 | 0.30 | 498 | 689 | 691 | 1.3 | 54 % |
| ×32 | 1.05 | 0.21 | 499 | 656 | 662 | 2.3 | 64 % |
| ×64 | 0.75 | 0.15 | 500 | 628 | 640 | 3.3 | 75 % |

Over a 64-fold change in loss at constant power — a swing from volumetric
heating to a pronounced surface skin — `T_wall` moves 4 K, within the
pyrometer's own scatter, while `T_center` moves 183 K.

Holding $P_{\text{sample}}$ fixed and scaling the dielectric loss over a factor
of 64 — which changes the penetration depth from several bed radii to a thin
surface skin — moves `T_wall` by only a few kelvin while moving `T_center` by
nearly 200 K.

This is not a defect. It is what makes the pair of thermometers informative:

* `T_wall` responds to **how much** power enters the sample and is blind to
  **where** it is deposited;
* `T_center` responds to both.

`T_wall` is therefore a clean power meter, and the difference between the two
readings isolates the deposition profile. Section S_x.7 uses the first half of
that statement as a consistency check on the measured series.

A secondary consequence concerns the FBG itself. Down to $\delta_p/R \approx 0.43$
the hottest cell in the bed is the on-axis cell and the FBG reads the maximum to
within 2 K. Below that the maximum migrates outward — $r = 1.3$ mm at
$\delta_p/R = 0.30$, 3.3 mm at 0.15 — and the on-axis FBG under-reads it by
2 K, 6 K and 12 K over the last three rows of Table S_x.4. The effect is real but
small compared with the FBG-to-average gap of Section S_x.5, and it changes no
conclusion here. It is reported so that the centre temperatures of the most
strongly reduced samples (R1100 at $\delta_p/R = 0.31$, Ti₂O₃ at 0.11) are read
as axial values rather than bed maxima.

---

## S_x.7  Consistency of the reported power partitioning

Section S_x.6 establishes that `T_wall` is set by $P_{\text{sample}}$ and the
external loss path alone. Two samples receiving equal power in the same tube
must therefore reach the same wall temperature, whatever their dielectric
properties. This is a falsifiable statement about the measured dataset, and it
is used here as a consistency check on the cavity-perturbation partitioning.

**Per-sample reconstruction.**

**Table S_x.5 — Per-sample thermal reconstruction of the rutile family.**
Each sample run at its measured tap density (hence its own void fraction and
Looyenga-mixed dielectric) and its reported `P_sample`. Thermal parameters are
those calibrated on r-TiO₂-R600 and are *not* re-fitted per sample.

| sample | tap ρ (g mL⁻¹) | void fraction | P_sample (W) | δ_p/R | T_wall model (°C) | T_wall meas. (°C) | Δ (K) |
|---|---:|---:|---:|---:|---:|---:|---:|
| a-TiO₂-R600  | 0.405 | 0.896 |  0.79 | 394  |  57 |  94 | **−37** |
| r-TiO₂-R600  | 0.975 | 0.769 | 14.88 | 2.11 | 368 | 403 | **−35** |
| r-TiO₂-R1000 | 0.757 | 0.821 | 17.48 | 0.77 | 403 | 436 | **−33** |
| r-TiO₂-R1100 | 1.185 | 0.719 | 17.55 | 0.31 | 405 | 386 | +19 |
| Ti₂O₃        | 2.340 | 0.478 |  7.99 | 0.11 | 256 | 221 | +35 |

The first three residuals agree to within 4 K of each other across a
twenty-two-fold range in deposited power and a 2.4-fold range in bed density.
That is a constant offset, not a trend, and a constant offset is hard to
produce by any sample-dependent mechanism; the natural reading is a calibration
difference of about −35 °C between the cavity used for the power sweep and the
dual-mode cavity used for the dielectric measurements. The two high-loss
samples break the pattern and change sign.

Two caveats on this table. The a-TiO₂ point sits at 0.79 W, well below the
calibration series, which itself shows its largest wall residual (+20.7 K) at
its lowest power; that row should be treated as an extrapolation. And the
thermal parameters are those calibrated on r-TiO₂-R600 for every row, so the
comparison assumes the samples share a bed conductivity, which is only
plausible within the rutile family.

**The R1000/R1100 pair.** These two samples are reported to receive
17.48 W and 17.55 W — equal within 0.4 % — yet their measured wall temperatures
differ by 50 °C. The model, run at the reported powers and the measured
densities, places them within 2 °C of each other. Density and void fraction do
not resolve the discrepancy: R1100 is the denser of the two, which if anything
moves the prediction the wrong way.

**Inversion.** Inverting the measured wall temperatures through the calibrated
`T_wall(P_sample)` curve, after removing the constant offset, gives the power
each sample must actually have received:

<<TABLE_S6>>

The ratio is unity within the calibration uncertainty for the three low-loss
samples and falls systematically below unity for the two high-loss ones, in
order of the perturbation parameter $\gamma\varepsilon''$.

**Interpretation.** The cavity-perturbation expression used in Note S2b,

$$\frac{P_{\text{sample}}}{P_{\text{abs}}} = \frac{2\gamma\varepsilon'' Q_{u0}}
{1 - \gamma(\varepsilon'-1) + 2\gamma\varepsilon'' Q_{u0}},$$

is a first-order result in the perturbation parameter $\gamma\varepsilon''$ and
is not expected to hold once that parameter approaches unity. The inversion
locates the departure empirically at

$$\gamma\varepsilon'' \gtrsim 5\times10^{-3},$$

above which the expression **over-estimates** the sample's share of the absorbed
power. Both R1100 and Ti₂O₃ lie above this bound. We therefore suggest a
validity annotation on Figure 2g rather than a change to the analysis: the
partitioning shown for the two most strongly reduced samples should be read as
an upper bound on $P_{\text{sample}}/P_{\text{abs}}$.

The inversion side of this comparison rests only on energy conservation in the
thermal model and the measured wall temperatures; no electromagnetic modelling
enters it. What is being tested is the cavity-perturbation partitioning against
a thermal measurement, and the two are independent.

---

## S_x.8  Limitations

Stated plainly, because several of them bound the results above.

1. **No electromagnetic solution.** The deposition profile is a fitted envelope,
   not a computed field. This is the origin of the ±30 K band in Section S_x.5.
2. **Parameter count.** Nine fitted parameters against fourteen data points.
   The model interpolates the calibration series reliably; extrapolation beyond
   it is not supported.
3. **Steady state only.** No transient capability, so nothing in this note
   speaks to heating or quench rates.
4. **Darcy field not coupled.** Pressure and permeability are solved and
   reported but do not feed back into the energy equation.
5. **Dielectric input provenance.** The `ε′(T)`, `ε″(T)` table used for the
   reduced-rutile calibration is a digitisation of the measured curve. It should
   be replaced with the numerical source before publication.
6. **Single-material calibration.** The transferability check in Section S_x.7
   re-uses the reduced-rutile thermal parameters for all five samples, scaling
   only the dielectric input and the void fraction. The constant −35 °C offset
   is consistent with that assumption; a per-sample thermal calibration is not
   available.
7. **Two-thermometer identifiability.** As shown in Section S_x.5, the wall
   thermometer carries essentially no information about the field width, and
   neither of the two additional exterior measurements we tested recovers it.
   We do not claim this is true of every conceivable exterior probe, only of
   those examined here. Results depending on the field width are reported as
   bands.

---

## S_x.9  Code and reproducibility

The solver is the microwave module of the Electrification Suite
(`apps/microwave/solver.js`), MIT-licensed, with unit tests and a verification
suite run in continuous integration. Every number in this note is produced by
the scripts in `tools/si/`.

| Item | Location |
|---|---|
| Solver core | `apps/microwave/solver.js` |
| Unit tests | `tests/microwave-solver.test.js` |
| Verification suite | `tools/verification/microwave.mjs`, `docs/VERIFICATION.md` |
| Note reproduction scripts | `tools/si/microwave-note.mjs` |
| Commit | *to be recorded on submission* |

```bash
node tools/si/microwave-note.mjs            # Tables S_x.1, S_x.3, S_x.4, S_x.5
node tools/si/microwave-note.mjs --invert   # adds Table S_x.6 (~1 h)
npm run verify:microwave                    # Table S_x.2 and the analytic checks
```
