# Viscosity Models

## GRD08 (Giordano, Russell & Dingwell, 2008)
Viscosity is computed using:

    log10(η) = A + B / (T_K - C)

Parameters A, B, C are computed from oxide mole fractions.

## Mueller et al. (2011)
Relative viscosity:

    η_rel = (1 - φc/φm)^(-2)

Maximum packing:

    φm = φm1 * exp(-(log10(rp))^2 / (2 b^2))

Aspect ratios:
- plagioclase: 4.0
- clinopyroxene: 1.5
- olivine: 1.5
- others: 1.0

Suspension viscosity:

    η_susp = η_rel × η_liq
