# Resolved square-channel Navier–Stokes and conjugate heat transfer

## Run

```sh
.venv/bin/python -m pip install -e '.[test]'
.venv/bin/python -m microwave3d.cli --config examples3d/itaca-monolith-ns.json --output results3d/monolith-ns
.venv/bin/python -m microwave3d.convergence --kind thermal --levels 2,3,4 --output results3d/channel-grid
.venv/bin/python -m microwave3d.convergence --kind flow --levels 2,3,4 --output results3d/flow-grid
# A separate, more expensive cavity-grid study:
.venv/bin/python -m microwave3d.convergence --kind em --levels 6,8,10 --output results3d/em-grid
```

Open the generated `viewer.html`, not the unpopulated template in `microwave3d/`.
Select **Channel gas** and **Navier–Stokes speed**, **Gas pressure**, or **Temperature**.
The numeric export includes sampled FEM velocity, conservative face-reconstructed
velocity, pressure, channel fraction, temperature, and deposited power.

## Geometry and evidence

The supplied Malhotra paper, DOI 10.1021/acs.iecr.0c05580, p. 6838 reports 21 square
channels, each 1.3 mm wide, in a 15 mm long monolith inside a 10 mm ID quartz tube.
The example uses these dimensions. The arrangement is an **assumed** 5 × 5 pattern
with four corners removed, at 1.65 mm pitch. Actual centre coordinates were not
tabulated. `channels.centres_m` accepts explicit coordinates relative to the tube
centre. Overlapping channels or corners outside the monolith are rejected.

The material example uses k = 50 W/(m K) from the SI. Its constant dielectric
value is an illustrative extension of the room-temperature value in Figure 9,
which describes a different 19-round-channel case. It is **not** a measured
temperature-dependent material table for the 21-channel experiment. Additional
electrical conductivity is zero to avoid adding a second unverified loss term.
Packed-bed porosity and particle diameter are ignored by the channel flow model;
void volume is calculated from the actual channel geometry.

## Momentum and continuity

The channel solver solves all three velocity components and pressure on a
tetrahedral mesh using quadratic velocity and linear pressure (Taylor–Hood P2/P1):

$$
\rho(\mathbf u\cdot\nabla)\mathbf u-\mu\nabla^2\mathbf u+\nabla p=0,
\qquad \nabla\cdot\mathbf u=0.
$$

Here ρ is reference gas density, μ is dynamic viscosity, u is velocity, and p is
gauge pressure. Length is scaled by channel width a, speed by the prescribed mean
speed U, and pressure by μU/a. The dimensionless convection coefficient is
Re = ρUa/μ. No permeability or Darcy resistance appears in these equations.

Side walls are no-slip. The inlet has a smooth biquadratic axial profile,
normalized using its FEM surface integral to give the specified mass flow.
The outlet uses zero natural traction for the grad-u weak formulation,
μ∂u/∂n − pn = 0. This is a do-nothing outlet, not a pointwise p = 0 condition.
The benchmark instead prescribes a fully developed profile at both ends and
fixes one pressure gauge degree of freedom. The fully developed profile is
computed from the square-duct Fourier series, not fitted to the numerical result.

Picard/Oseen iteration retains convective inertia and checks the nonlinear
momentum residual and velocity change. A failed linear solve or iteration raises
an error. The current unstabilized implementation rejects Re > 200. This limit
is a numerical scope guard, not a universal transition criterion.

Identical isolated channels receive equal mass flow. A representative 3D solve
is reused by translational symmetry; this does not replace it with a 1D pressure
law. Channel-specific inlet maldistribution and shared manifolds are not solved.
Gas density and viscosity are fixed at `gas.reference_temperature_c`. Heat changes
the thermal conductivity and configured dielectric properties, but does not
change flow density or viscosity. Consequently, large temperature rises are
outside a justified constant-density approximation. Variable-density low-Mach
flow remains an explicit extension, not an implemented capability.

## Heat transfer and conservative coupling

A separate Cartesian thermal mesh includes every EM cell boundary and every
square-channel boundary. Channel walls are therefore exact grid planes. The
circular outer monolith and quartz still use subcell volume quadrature. Solid
and gas occupy separate cells with separate temperatures; there is no packed-bed
local-thermal-equilibrium assumption in this option.

Solid/gas interfaces use harmonic normal conduction resistance. This imposes
heat-flux continuity without a fitted wall heat-transfer coefficient. The inlet
has prescribed temperature for diffusion and incoming enthalpy. The outlet has
zero diffusive flux and outgoing enthalpy. Flow stops at the monolith end planes;
the surrounding tube gas conducts heat but is not a resolved flowing manifold.
Radiation retains the original outer-quartz gray-enclosure approximation.

FEM velocity is integrated on thermal faces with 2 × 2 Gaussian points.
Continuous P2/P1 velocity is weakly divergence-free, not exactly conservative on
arbitrary finite-volume cells. A graph-based flux projection enforces each
thermal cell's mass balance while preserving each prescribed end flow. Both
the correction's relative L2 norm and the final cell mass residual are reported.
Large transfer corrections indicate a need to refine the flow/thermal grids.
The pressure field is the FEM pressure, not the projection multiplier.

Thermal material properties are volume-averaged onto their parent EM cells.
EM sample/quartz absorption is transferred back separately in **watts**, weighted
by each fine cell's dielectric-plus-ohmic loss coefficient and volume at that
frequency and temperature. Sample heat is never deposited in pure gas
cells. Each parent cell's absorption is conserved. Coarse EM voxels still use
arithmetic dielectric homogenization; thermal channel resolution does not imply
channel-resolved electromagnetics. Refining the channel thermal mesh alone cannot
remove the coarse cavity resonance or port error.

## Verification and uncertainty

The test suite compares fully developed pressure drop with the square-duct series:

$$
\Delta p=28.45415377\,\mu L U/a^2.
$$

Here L is channel length. Tests also check no slip, weak continuity, nonlinear
momentum residual, exact phase-wise heat transfer, zero-power equilibrium,
positive gas enthalpy rise, and final coupled power/mass closure.

The convergence command saves every run and reports changes in maximum
temperature, absorption, and pressure drop. At least three increasing levels
are required. `last_pair_passed` means only that the final tested pair meets
the specified tolerances; it is neither an error bound nor validation. An EM
study changes the nested thermal mesh too, as stated in its report.

The example selects an equilibrated GMRES thermal solve with an incomplete-LU
preconditioner. Reuse affects only the preconditioner: the current thermal matrix
is still assembled and its residual checked on every iteration. A regression
test compares the solution with direct sparse factorization. Set
`solver.thermal_linear_solver` to `direct` to use the original factorization.

The feed remains an ideal lumped probe. Real antenna metal, coax geometry,
windows, measured dielectric curves, thermal expansion, and experimental
validation have not been added by this channel extension.
