# MELTS Streamlit GUI — Workflow Overview

The GUI provides three major modes:

## 1. Liquidus Temperature
- MELTSdynamic is initialized
- Bulk composition is loaded
- Fe redox is normalized using Fe3+/FeT slider
- Pressure and fO2 buffer are applied
- A coarse + fine temperature descent finds the first appearance of solids

## 2. Sub-liquidus Equilibrium
- MELTS is run at user-defined temperature
- Melt fraction and solid fraction are computed
- Solid phases are listed with wt%
- Near-liquidus smoothing is applied when ΔT < 8 °C
- Optional Mueller viscosity:
  - Aspect ratios assigned by phase
  - φm computed from Mueller Eq. 4
  - Relative viscosity computed from Eq. 3
  - Suspension viscosity = η_rel × η_liq

## 3. Liquid Viscosity (GRD08)
- GRD mole fractions computed from wt%
- VFT parameters (A, B, C, Tg, F) computed
- Viscosity evaluated at:
  - Liquidus
  - User temperature
  - Temperature range
- Results plotted using Streamlit charts

## Data Tables
The following tables are included:
- Bulk_comp_tbl.txt
- Liquid_comp_tbl.txt
- Solid_comp_tbl.txt
- Phase_main_tbl.txt
- System_main_tbl.txt

These are used internally by MELTSdynamic.
