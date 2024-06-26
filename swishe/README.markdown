# SWISHE: (S)uppression of (w)ind-(i)nduced (s)urface (h)eat (e)xchange

Last updated: 2023-06-14

## What is this?
SWISHE is an attempt to suppress tropical cyclone (TC) formation by limiting evaporation in global climate models (GCMs) for cells that correspond to TCs. The theoretical approach here is:
- mesoscale: to investigate how TCs can intensify with a suppression of vertical moist enthalpy transfer
- large-scale: investigate what role TCs play in climatic processes by suppressing them (mechanism denial)i

## Directory information
The directory herein is a (probably feeble and naive) attempt at modifying surface flux calculations to incorporate TC-specific criteria, such as vorticity thresholds and warm-core identification, to locate where to suppress evaporative winds (SWISHE, in other words).

## Notes to self
- Latest commit features vorticity data being passed from `fv_diagnostics.F90` to `surface_flux.F90`. Climate data appears to be equivalent to climate data from Wenchang's original KillTC runs, with TC suppression occurring (64 to 25 tropical-storm strength vortices, 25 to 9 hurricane-strength vortices) despite potential Findlater jet suppression. Note that vorticity criteria are currently only being passed for the NH and latitude constraints have not been imposed yet.
