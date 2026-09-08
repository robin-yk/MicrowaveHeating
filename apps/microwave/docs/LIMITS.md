# What this model does not do

## Microwave 2D field

### Verification is by analytic reduction, not a manufactured solution

A full-domain manufactured solution is deliberately out of scope: the tube-gas
and outside-air conductivities are temperature dependent by design and not
overridable through public parameters, so a uniform-conductivity domain cannot
be produced without editing shipped code. The implementation is verified
against analytic pure conduction and grid refinement instead, and the shared
finite-volume pattern is manufactured-solution verified in the Joule solver.

### The dielectric response is fitted, not predicted

ε′(T) and ε″(T) come from tabulated measurements, mixed to the bed with
Looyenga or Maxwell-Garnett. Bed conductivity is fitted to FBG measurements
through `k200`, `k500`, `k800` and a radiating-area factor.

**Do not** transfer a calibration to a different powder, packing, or applicator
geometry. The fit absorbs whatever the model does not represent.

### The homogenization has a stated validity window

`homogenizationValidity()` checks the unit cell against the macroscopic length
and the wavelength. Outside that window the continuum treatment of the packed
bed is not justified, and the check exists so that condition is visible rather
than assumed.
