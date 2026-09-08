# Numerical verification of the microwave 2D solver

Extracted from the Electrification Suite verification report; the studies below
are the ones that exercise `solver.js`.

## Methods

Four standard techniques are used:

1. **Discrete exactness**: a uniform volumetric source with constant
   conductivity has an exactly parabolic radial temperature profile, and a
   conservative two-point-flux FV scheme reproduces quadratics exactly at cell
   centers. Any mismatch is therefore *physics left in the setup* (finite-rod
   axial leakage), not discretization error, and must vanish as the element is
   made longer.
2. **Analytic multi-layer benchmark**: mid-plane temperature drops across the
   element / wall / surrounding-air layers versus the ln-resistance solution
   of an infinite cylinder.
3. **Manufactured solution (MMS)**: a smooth exact solution
   T\*(r,z) = T_a + A·(1−(r/R_d)²)²·(1−(2z/H_d)²)² is imposed by injecting the
   analytically derived source q = −k∇²T\* on a uniform-conductivity domain
   (element k = wall k = gap k = outside-air k) with radiation and convection
   off. The quartic bump vanishes with zero slope on every outer boundary, so
   the solver's ambient boundary handling is exactly consistent with T\*. L2
   and L∞ errors against T\* measure the observed convergence order.
4. **Grid sensitivity of the shipped default cases**: the full nonlinear
   models (radiation, He purge flow, temperature-dependent properties) on a
   doubling grid sequence.

## Microwave solver (`solver.js`, `solve2D`)

### Scope note

A full-domain manufactured solution is deliberately out of scope for this
solver: the tube-gas and outside-air conductivities are temperature-dependent
by design and not overridable through public parameters, so a
uniform-conductivity domain cannot be produced without modifying shipped
code. The microwave implementation is instead verified against analytic
conduction (below) and grid-refinement behavior; the shared FV pattern
(two-point flux, harmonic interface resistance) is MMS-verified in the Joule
solver.

### 1. Radial parabola: pure-conduction reduction

Public knobs alone reduce solve2D to pure conduction with a uniform bed
source: manual constant bed conductivity, radiation off, gas exchange off,
flat field shape, and a near-infinite penetration depth (the source
normalization keeps total power exact). Mid-plane bed profile vs the exact
parabola:

| case | worst relative mismatch |
| --- | --- |
| H/D = 4 | 8.7e-2 |
| H/D = 8 | 2.4e-2 |
| H/D = 16 | 2.3e-3 |

Same signature as the Joule study: the residual is finite-bed axial leakage
and collapses as the bed is made longer; the bed conduction discretization
matches the analytic profile.

### 2. Default calibrated case: grid sensitivity

The app's default reduced-rutile TiO₂ calibration case (26 W absorbed, He
flow, Looyenga dielectric mixing, wall radiation, Darcy flow on):

| grid | center T (°C) | wall T (°C) | avg bed T (°C) | energy closure | Darcy mass imbalance |
| --- | --- | --- | --- | --- | --- |
| 30×60 (default) | 811.20 | 495.71 | 622.79 | 2.4e-4 | 1.8e-13 |
| 60×120 | 804.65 | 495.74 | 621.85 | 5.0e-4 | 2.1e-13 |
| 120×240 | 788.98 | 496.13 | 618.23 | 1.0e-3 | 2.1e-13 |

Energy closure holds to ≤1e-3 of the absorbed power and the Darcy flow field
conserves mass to machine precision (~1e-13 relative) on every grid. As in
the Joule case, the sequence is not yet asymptotic (the near-wall exponential
source deposition and T⁴ wall radiation are resolved progressively), so a
sensitivity bound is reported instead of an order: **the default 30×60 grid
differs from the 120×240 grid by 22.2 K at the bed center (2.9% of the rise),
4.6 K in average bed temperature (0.8%), and 0.4 K at the wall**. The wall
temperature, which is what the experimental calibration constrains, is
grid-insensitive.

## Microwave field solve (`solver.js`, `solveField2D`)

Opt-in through `p.fieldMode = "helmholtz"`; the shipped default still uses the
fitted source, because switching over changes every number the page reports and
the bed conductivities on the Calibration tab were fitted against the old shape.

### Why a field solve is tractable here

At 2.404 GHz the free-space wavelength is 124.7 mm, and SiC at eps' = 7.96 brings
it to 44.2 mm inside the bed. A 10 mm bed therefore spans **D/lambda = 0.23**, so
the load supports no internal cavity mode and nothing outside it needs meshing,
and the existing 0.5 mm cells already give **88 per wavelength**. What is solved
is the scalar Helmholtz problem for an axial E,

    (1/r) d/dr ( r dE/dr ) + d2E/dz2 + k0^2 eps(T,r,z) E = 0

on the temperature mesh, as a complex system in real 2N block form with an exact
per-cell 2x2 Jacobi preconditioner. Stated limits: the scalar form is exact for
eps varying with r alone and drops a grad(eps) coupling; the incident field is
taken uniform and imposed on the domain boundary; and the absolute coupling
efficiency is still not predicted, so the total is renormalised to the absorbed
power exactly as the fitted shape was.

### Control

An infinite lossy dielectric cylinder in a uniform axial field has a closed form,
`E(r)/E(R) = J0(kr)/J0(kR)` with `k = k0 sqrt(eps)` and J0 of a complex argument.
At |kR| = 0.711 the solved profile reproduces it and converges cleanly:

| grid | axis \|E\|/\|E(R)\| | exact | error |
| --- | --- | --- | --- |
| 15x30 | 1.109556 | 1.109332 | 2.02e-4 |
| 30x60 | 1.124257 | 1.124193 | **5.67e-5** |
| 60x120 | 1.131765 | 1.131749 | 1.47e-5 |

Error falls by 3.6x then 3.9x per halving — second order.

### What it changes, and a correction

Power density on the axis relative to the bed edge:

| | axis / edge |
| --- | --- |
| exact (refraction) | 1.264 |
| shipped fitted source | 1.809 — **43% too peaked** |
| its Beer-Lambert skin alone | 0.969 |

The centre peaking is real physics: refraction into a subwavelength load. An
earlier reading of this comparison quoted the skin factor alone, concluded the
shipped source was edge-peaked and therefore had the sign wrong, and was mistaken
— the fitted Gaussian in front of it dominates and peaks the source hard on the
axis. The skin term does have the sign backwards, but at delta = 140 mm against a
5 mm radius it moves the source by 3.2% and nothing rests on it.

The real finding is that the shipped model obtains a genuine physical effect from
a **fitted width** rather than from the field, and overshoots it by half again.
The temperature consequence is modest — centre 473.8 -> 469.6 C, spread 35.7 ->
29.7 K on the SiC default — but `fieldWr` stops being a free parameter, and the
bed conductivities fitted against a 43%-too-peaked source were absorbing that
error, so refitting them against the solved field is what makes those numbers
mean what they claim.

## Microwave grid convergence, and the gas-exchange term that breaks it

With the Krylov solve in place the grid sequence can be measured honestly. On the
reduced-rutile default case it does not converge:

| grid | centre (°C) | avg bed (°C) | gas outlet (°C) | energy closure |
| --- | --- | --- | --- | --- |
| 30×60 | 811.16 | 622.75 | 520.39 | 4.25e-7 |
| 60×120 | 804.56 | 621.76 | 514.61 | 1.31e-7 |
| 120×240 | 788.80 | 618.05 | 509.19 | 8.07e-8 |

Centre differences **grow**, −6.60 then −15.76, and the gas outlet drifts by a
near-constant 5.8 K then 5.4 K rather than halving. Energy closure is excellent
throughout, so nothing is being lost — the answer simply keeps moving.

Switching the bed-to-gas exchange off isolates it completely:

| grid | centre (°C), no gas exchange |
| --- | --- |
| 30×60 | 1061.28 |
| 60×120 | 1066.61 (+5.34) |
| 120×240 | 1069.30 (+2.68) |

Monotone, ratio 1.99, **first order**. The rest of the solver converges cleanly.

### Why

The bed-to-gas march gives each bed row a finite-effectiveness stage,

    eff = 1 − exp( −(UA_total / rowCount) / C_gas )

and `rowCount` is the number of bed rows, which is the grid. Refining the mesh
therefore changes the number of stages in the exchanger, so it is not one model
being solved more accurately, it is a different model on every grid. The total UA
is preserved, but a cascade of N finite stages is not the same object as the
continuous exchanger it is standing in for, and the temperature field feeds back
into each stage's driving difference, so the drift does not settle at the O(1/N)
a decoupled cascade would give. This is the same species of defect as the He
purge stream in the Joule solver, whose ambient clamp scaled as 1/(dz/2) and
drove that solver's observed order negative.

### What it does and does not invalidate

The defect bites hardest where axial gas transport competes with radial
conduction. Reduced rutile has a bed conductivity of 0.6–1.7 W/m·K and is
dominated by it; SiC at 4–18 W/m·K is nearly isothermal across the bed and is
not, which is why the SiC sequence used for the recalibration converges cleanly
(473.56 → 476.04 → 477.03, differences +2.48 then +0.99, p ≈ 1.3) with a 30×60
grid error near 3.5 K, under the 5.5 K residual being fitted.

So the SiC refit stands, and **the reduced-rutile profile must not be recalibrated
until the exchange term is made grid-independent.**

### Refinement has a floor, and it is material-specific

`bedHomogenization()` is now evaluated on every solve and the grid study refuses
levels that cross it. The mesh cell must stay coarser than the packing unit cell,
or the grid resolves structure the continuum model does not carry:

| grid | SiC, d_p = 194 µm | TiO₂, d_p = 50 µm |
| --- | --- | --- |
| 30×60 | h/d_p = 2.58 | 10.0 |
| 60×120 | 1.29 | 5.0 |
| 120×240 | **0.64 — below** | 2.5 |

A SiC grid study may reach 60×120 and no further. The sequences above are the
TiO₂ profile, which has room at 120×240.

## Microwave calibration is fitting mesh error (`tools/verification/microwave-calibrate.mjs`)

The page fits four numbers -- the bed conductivity anchors k200, k500, k800 and
the radiation-area multiplier radArea -- by bounded coordinate search against a
measured [P_abs, T_wall, T_centre] sweep. The search runs on a 10x30 mesh and the
result is reported on 30x60. At fixed parameters those two meshes disagree by far
more than the residual being minimised:

| P (W) | centre 10x30 -> 30x60 | wall 10x30 -> 30x60 |
| --- | --- | --- |
| 10 | 225.0 -> 217.2 (−7.8) | 188.4 -> 185.5 (−2.9) |
| 20 | 340.0 -> 325.5 (−14.5) | 277.0 -> 271.3 (−5.7) |
| 30 | 427.4 -> 406.4 (−21.0) | 339.7 -> 331.8 (−7.9) |
| 40 | 501.0 -> 473.8 (−27.3) | 389.0 -> 379.3 (−9.7) |

**The discretisation error reaches 27 K where the RMSE being minimised is about
10 K**, and it grows systematically with power rather than scattering. A search
under those conditions does not find parameters that describe the physics; it
finds parameters that cancel coarse-mesh error. The reproduction shows exactly
that: driving the coarse objective from 10.21 down to 7.87 leaves the 30x60
result at 10.83, *worse* than the 8.72 the unsearched defaults give.

### What this invalidates, including a claim made here

Comparing the fitted-Gaussian source against the solved field on this objective
does not work yet, and an intermediate reading of it was wrong. On the 10x30
search mesh the solved field scored 8.80 against 10.21 and looked like a 14%
improvement; on 30x60 the order reverses, 11.23 against 8.72. Neither number
means anything, for two reasons that have to be fixed in order:

1. the objective is dominated by discretisation error, and
2. the shipped conductivities were fitted *with* the Gaussian source, so they
   flatter it. Each source model has to be refit independently before either
   can be scored.

### Resolved: the search mesh was the problem, not the report mesh

`solve2D` now assembles its linearised system and hands it to a preconditioned
conjugate-gradient solve, ported from the Joule solver. The matrix is symmetric
by construction — every interior face contributes the same conductance to both
its cells — and diagonally dominant once boundary, gas-exchange and radiation
terms are added, so CG applies directly. It reproduces the relaxation answer to
0.2 K (10x30 centre 501.0 -> 500.88; 30x60 473.8 -> 473.56), and makes the fine
meshes affordable:

| grid | centre (°C) | Picard steps | wall clock |
| --- | --- | --- | --- |
| 10x30 | 500.88 | 80 | 0.26 s |
| 30x60 | 473.56 | 36 | 1.1 s |
| 60x120 | 476.04 | 21 | 3.6 s |
| 120x240 | 477.03 | 16 | 20 s |

That locates the fault precisely. The converged centre is about 477 °C, so the
**10x30 search mesh was 24 K away while 30x60 is 3.5 K away** — the report mesh
was always adequate and only the search mesh was not. Moving the search to 30x60
puts discretisation error under the residual, and at about a second per solve a
full coordinate search costs minutes.

### The refit, and what the fitted conductivities were carrying

With both source models independently refit on the 30x60 mesh under an identical
search schedule:

| | k200 | k500 | k800 | radArea | combined RMSE |
| --- | --- | --- | --- | --- | --- |
| shipped defaults | 4.000 | 18.00 | 18.00 | 6.000 | 8.85 |
| refit, fitted Gaussian source | 2.761 | 12.45 | **27.50** | 5.868 | 5.67 |
| refit, solved Helmholtz field | 2.879 | 11.57 | **13.93** | 5.802 | **5.50** |

Two results. The shipped calibration was **60% worse than achievable** (8.85
against 5.50) purely because it had been searched on a mesh whose error exceeded
the residual. And the hypothesis that motivated the field solve is confirmed and
localised: k200 and k500 barely distinguish the two source models, but **k800
differs by a factor of two**. The high-temperature anchor is where eps''(T) is
largest and the source shape matters most, and fitting against a 43%-too-peaked
source drives k800 to 27.5 — near its upper bound of 30 — to spread heat the
source should never have concentrated. Against the solved field the same data
asks for 13.9.

The solved field also fits marginally better, 5.50 against 5.67, now that the
comparison is like-for-like.
