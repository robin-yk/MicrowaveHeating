# Input cases

`itaca-cylinder.json` uses the supplied ITACA cavity envelope, 104.92 mm diameter
and 85 mm height. The sample is an **illustrative effective SiC packed bed** in a
10 mm ID tube, 15 mm long. Its material table, probe, wall metal conductivity,
quartz properties, and thermal enclosure properties are assumptions. They are
not fitted to the Malhotra paper and are not the paper's channelled monolith.

All lengths are metres, powers watts, temperatures Celsius, frequencies Hz,
pressure Pa, conductivity W/(m K) or S/m as explicitly named, gas flow sccm at
the 0 °C, 1 atm molar volume 0.022414 m³/mol.

The four material-table columns are temperature, real relative permittivity,
relative dielectric loss, and effective sample thermal conductivity. Values
interpolate linearly only within the supplied range. An optional
`material.electrical_conductivity_s_m` adds ohmic loss consistently in both
Maxwell and heat deposition. Set `loss_definition` to `total_measured` for a
dielectric table that already includes conductivity; adding a nonzero separate
conductivity is then rejected to prevent double counting.

The supplied paper provides several distinct cases, which must not be merged:

| Case | Sourced information | Missing input for reproduction |
| --- | --- | --- |
| Solid SiC benchmark | D = 10 mm, H = 15 mm; Figure 2 power–temperature comparison; 8 ± 1 W empty-tube offset | Numerical property curves, actual port geometry, digitized/raw comparison data |
| Experimental honeycomb | 21 square 1.3 × 1.3 mm channels; 15 mm length; 10 mm ID tube; Figure 3 at 2410 MHz and 125 W input | Exact channel placement, wall dimensions, channel-resolved flow mesh and boundary data |
| Figure 9 parameter study | 19 round channels, radius 0.5 mm; ε* = 9.72 − 2.3j at 300 K for SiC; caption gives k = 100 W/(m K), σ = 50 S/m | Full temperature-dependent dielectric table and matched numerical geometry |

SI Section III separately discusses a conservative k = 50 W/(m K), while the
Figure 9 caption uses 100 W/(m K). These are different cases, not interchangeable
calibration values. No monolith preset is shipped that silently fills channels
with a porous powder model.

The first-order tube-aperture boundary is an outgoing impedance approximation.
The ideal lumped radial port has a defined 50-ohm reference; it does not reconstruct
the coax antenna or its matching network. Probe length and position are editable.

Use denser frequency sampling around an observed resonance. A handful of equally
spaced points can miss a narrow resonance, so averaged absorption must be checked
for convergence in frequency as well as in space. Dwell weights describe time
spent at each frequency, not coherent field amplitudes.
