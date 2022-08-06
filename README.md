# AstrodynamicsBase

A set of base routines for astrodynamics. 

## Installation

```julia-repl
(@v1.7) pkg> dev git@github.com:Yuricst/AstrodynamicsBase.jl.git
```

## Capabilities
- [x] Definition of canonical scales
- [x] Solar system constants calls (GM, SMA, SOI)
- [x] Cartesian/Keplerian/MEE/Poincare representation of two-body states
- [ ] Keplerian orbit parameter calls (period, fpa...)
- [x] Transformation between inertial and rotating systems
- [x] Wrap to SPICE functions for ephemeris handling & `spkssb`
- [x] Kepler's problem
