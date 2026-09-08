# Microwave3D research solver

This addition solves a vector electromagnetic field and a coupled three-dimensional
thermal field locally using NumPy and SciPy. It also solves optional three-dimensional
Darcy flow in a homogenized packed bed and puts its conservative face mass fluxes
into the energy equation. It does **not** implement channel-resolved Navier–Stokes
CFD of the monolith in the supplied paper. Existing microwave2D files are retained.

## Run locally

Python 3.11 or newer is required. From the repository root:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -e '.[test]'
.venv/bin/python -m pytest -q
.venv/bin/python -m microwave3d.cli --power 130 --output results3d/itaca
```

Open `results3d/itaca/viewer.html` to inspect the actual calculated volume and
XY/XZ/YZ slices. This viewer is self-contained and works offline; it does not
recalculate the physics. `result.json` records inputs, provenance, convergence,
and power accounting. `fields.npz` retains numeric arrays, including complex
electric components per frequency. `fields.vtk` opens in a scientific field viewer.

The default is a coarse demonstration mesh. Its temperature is not a validated
prediction for the experimental apparatus. The same warning remains in every
export, even if the nonlinear equations converge.

```sh
# Incident-power continuation; input order is preserved, including descending runs.
.venv/bin/python -m microwave3d.cli --powers 0,30,60,90,130 --output results3d/power

# Cold-cavity frequency scan, not a temperature calculation.
.venv/bin/python -m microwave3d.cli --spectrum-only --power 130 \
  --frequencies-mhz 2350,2375,2400,2425,2450,2475,2500 --output results3d/cold

# Fast-sweep, equal-dwell average at each nonlinear thermal iteration.
.venv/bin/python -m microwave3d.cli --power 130 \
  --frequencies-mhz 2400,2402,2404,2406,2408 --output results3d/band

# Independent geometric benchmark for an empty closed PEC cylinder.
.venv/bin/python tools/verification/microwave3d.py --grids 8,12,16,24
```

Power sweeps also write
`power-sweep.csv` and `power-sweep.html`, with separate centre, sample-average,
and quartz temperature curves against incident power. The table retains actual
calculated sample absorption at each point.

A failed solve exits with status 2. A finite, unconverged result is saved with
`converged: false`, and power continuation stops there. No temperature clipping
converts an out-of-range state into apparent success. Property extrapolation is
rejected. Physical stability and missing steady branches are not diagnosed by
fixed-point convergence.

## Electromagnetic model

The peak electric phasor E uses exp(−iωt). The constitutive convention is
ε* = ε′ + iε″ for passive dielectric loss. The loss sign in the supplied paper
uses the opposite phasor convention; the physical dissipated power is unchanged.

$$
\nabla\times\mu_0^{-1}\nabla\times\mathbf E
-\omega^2\epsilon_0\epsilon'\mathbf E
-i\omega(\omega\epsilon_0\epsilon''+\sigma_e)\mathbf E=0.
$$

Here ω = 2πf, f is frequency in Hz, μ₀ and ε₀ are vacuum permeability and
permittivity, ε′ and ε″ are relative dielectric properties, and σₑ is optional
electrical conductivity. The implementation uses all three electric components
on an orthogonal Yee grid. It is not three independent scalar Helmholtz solves.
Edge voltages and oriented face circulations form a discrete curl-curl operator.
The nonuniform grid has a refined central region and a stepped cylindrical wall.

The lumped port voltage V is a weighted edge-line integral. The configured probe
length is preserved with a fractional terminal edge. With port conductance G,
incident peak voltage V⁺ = √(2Pᵢₙ/G), and reflected voltage V⁻ = V − V⁺:

$$
S_{11}=V^-/V^+,
\qquad P_{\mathrm{ref}}=P_{\mathrm{in}}|S_{11}|^2.
$$

The passive port contributes −iωGbbᵀ to the matrix and
−2iωGbV⁺ to the right-hand side, where b extracts V from the edge voltages.
This is an ideal circuit port coupled to the field. The measured coax matching
network, feed metal, and viewing windows are not represented.

For a finite-conductivity cavity wall, the surface impedance is
Zₛ = (1−i)√(ωμ₀/(2σwall)). Surface absorption is integrated separately. Optional
first-order scattering at the tube end apertures uses the vacuum wave admittance.
It approximates an outgoing boundary; it is not a meshed external waveguide or PML.

$$
q'''_s=\tfrac12(\omega\epsilon_0\epsilon''_s+\sigma_{e,s})|\mathbf E|^2,
\qquad
P_{\mathrm{in}}=P_{\mathrm{ref}}+P_s+P_q+P_{\mathrm{wall}}+P_{\mathrm{aperture}}.
$$

The symbols Pₛ and Pq denote sample and quartz absorption. Power escape through
the tube apertures is P_aperture. The dielectric mass and thermal deposition use
the same positive cell quadrature, so their integrated losses agree to round-off.
No prescribed sample absorption or rescaling to a temperature target is applied.

For frequency dwell weights aⱼ ≥ 0 with Σaⱼ = 1, solve each frequency at the
current temperature and use Σaⱼqⱼ as the average heat source. Complex fields at
different frequencies are never added coherently. This is appropriate only when
temperature changes little during a frequency sweep. Repeating four identical
sweeps does not multiply average power by four. Frequency and dwell resolution
remain user-controlled and must be converged around narrow resonances.

## Thermal and porous-flow model

The thermal unknown is one continuum cell temperature. In the sample it is a
local thermal-equilibrium (LTE) approximation. It cannot predict a pointwise
solid–gas temperature difference in a channel or particle.

$$
\nabla\cdot(\rho_g c_{p,g}\mathbf u T)
-\nabla\cdot(k\nabla T)=q'''_s+q'''_q.
$$

Gas density ρg, specific heat cp,g, superficial velocity u, and conductivity k
are evaluated under the chosen material/gas assumptions. Conservative harmonic
face resistances give conduction. The exact same mass flux from the Darcy solve
gives upwind advection, avoiding the grid-dependent radial mixing in the old
slice-effectiveness closure. Outside the bed, quartz, stagnant tube gas, and
cavity air conduct heat. Open-tube and external buoyancy flows are not solved.

$$
\nabla\cdot(\rho_g\mathbf u)=0,
\qquad \mathbf u=-\mathbf K\nabla p/\mu_g.
$$

The diagonal permeability tensor K uses the packed-bed viscous Ergun expression
ε³dₚ²/[150(1−ε)²] and independent directional multipliers. Here ε is void fraction,
dₚ particle diameter, p gas pressure, and μg viscosity. Isobaric bed end faces and
the measured total mass flow determine the pressure drop and distribution. The
implementation rejects Δp/p > 1% and particle Reynolds number > 1. These guards
bound the reduced flow assumptions; they do not validate LTE or homogenization.

The outer quartz surface exchanges diffuse-gray radiation with an isothermal
gray cavity enclosure. A uniform enclosure radiosity closes multiple reflection
between a convex tube and cavity. Total tube area is geometric, with no fitted
area multiplier. Quartz transmission, internal sample radiation, and gas radiation
are excluded. Those omissions are particularly relevant when reproducing the
paper's spectrally transmitting quartz measurements. Cavity metal absorption is
rejected to the prescribed cooled wall and is not counted again as sample heating.

$$
P_s+P_q=Q_{\mathrm{boundary}}+Q_{\mathrm{radiation}}+Q_{\mathrm{gas}}.
$$

All terms are recomputed at the final temperature, including a fresh Maxwell
solve, material properties, Darcy fluxes, and radiation. Acceptance requires an
undamped temperature-step tolerance, a global heat residual, and a cell residual.
Energy balance alone does not establish mesh convergence or measurement agreement.

Subcell material fractions prevent thin quartz from disappearing on coarse grids.
They use arithmetic property mixing within a voxel and therefore introduce an
interface approximation. The optional contact conductance applies only to faces
separating predominantly sample and quartz cells; it is off by default and is not
a subcell interface model. Sample-average temperature is volume weighted. The
reported centre is the closest sample cell, not an FBG-length average. Wall
temperature is an area average around the quartz midplane, not a pyrometer model.

## Source accounting and validation boundary

Source: Malhotra et al., *Temperature Homogeneity under Selective and Localized
Microwave Heating in Structured Flow Reactors*, Ind. Eng. Chem. Res. 2021, 60,
6835–6847, DOI [10.1021/acs.iecr.0c05580](https://doi.org/10.1021/acs.iecr.0c05580),
main paper and `ie0c05580_si_001.pdf` supplied by the user. The PDFs are not
redistributed in this repository.

The paper supplies the cavity envelope, impedance-wall/scattering-boundary
approach, a 10 × 15 mm solid-cylinder benchmark, and distinct channelled monolith
cases. It reports empty-tube accepted power of 7–11 W at 20 MHz bandwidth and an
8 ± 1 W offset for one experimental comparison. This does not establish a
universal cavity loss at every loading, frequency, and temperature. The code
therefore exposes the power partition instead of subtracting 8 W automatically.

The quoted P_input(1−S11) is ambiguous unless S11 denotes a power fraction. The
code defines S11 explicitly as an amplitude ratio and consistently uses |S11|².

Before use as an experimental prediction, supply the coax/probe and window
geometry or a validated port equivalent, actual metal-wall loss, measured complex
dielectric curves with a clear conductivity convention, specimen-specific thermal
properties, thermal enclosure conditions, and spatial/temperature reference data.
Validate empty and loaded S11, frequency-averaged absorption, and temperatures on
independent runs. Resolve monolith channels and laminar fluid mechanics before
claiming reproduction of the paper's channel gas temperatures.
