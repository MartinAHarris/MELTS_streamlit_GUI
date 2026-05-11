# -*- coding: utf-8 -*-
"""
Revised MELTS Liquidus + Sub-liquidus + GRD Viscosity + Mueller 2011 GUI
Updated: May 2026
Author: Martin Harris, PhD
"""

import math
import pandas as pd
import streamlit as st
from meltsdynamic import MELTSdynamic

# ---------------------------------------------------------
# GRD08 (Giordano, Russell & Dingwell, 2008) implementation
# ---------------------------------------------------------

def _mole_fractions_grd(wt):
    """
    Convert wt% oxides to mole percent in GRD order:
    [SiO2 TiO2 Al2O3 FeO(T) MnO MgO CaO Na2O K2O P2O5 H2O F2O-1]
    """
    M = {
        'SiO2': 60.0843,
        'TiO2': 79.866,
        'Al2O3': 101.961,
        'FeO': 71.844,      # FeO(T)
        'MnO': 70.937,
        'MgO': 40.304,
        'CaO': 56.077,
        'Na2O': 61.979,
        'K2O': 94.196,
        'P2O5': 141.944,
        'H2O': 18.015,
        'F2O-1': 37.0       # dummy; wt% usually 0
    }

    order = ['SiO2', 'TiO2', 'Al2O3', 'FeO', 'MnO',
             'MgO', 'CaO', 'Na2O', 'K2O', 'P2O5', 'H2O', 'F2O-1']

    moles = []
    for ox in order:
        w = wt.get(ox, 0.0)
        moles.append(w / M[ox] if M[ox] > 0 else 0.0)

    tot = sum(moles)
    if tot <= 0:
        return [0.0] * len(order)

    return [m / tot * 100.0 for m in moles]


def grd_params_from_wt(wt):
    """
    Compute GRD08 VFT parameters A, B, C, Tg, F from wt% oxides.
    wt: dict with wt% of SiO2, TiO2, Al2O3, FeO(T), MnO, MgO,
        CaO, Na2O, K2O, P2O5, H2O, (optional) F2O-1.
    Returns (A, B, C, Tg, F)
    """
    AT = -4.55
    bb = [159.56, -173.34, 72.13, 75.69, -38.98,
          -84.08, 141.54, -2.43, -0.91, 17.62]
    cc = [2.75, 15.72, 8.32, 10.2, -12.29, -99.54, 0.3]

    xmf = _mole_fractions_grd(wt)
    # indices: 0 SiO2, 1 TiO2, 2 Al2O3, 3 FeO(T), 4 MnO,
    #          5 MgO, 6 CaO, 7 Na2O, 8 K2O, 9 P2O5, 10 H2O, 11 F2O-1

    siti = xmf[0] + xmf[1]
    tial = xmf[1] + xmf[2]
    fmm  = xmf[3] + xmf[4] + xmf[5]
    nak  = xmf[7] + xmf[8]

    b1  = siti
    b2  = xmf[2]
    b3  = xmf[3] + xmf[4] + xmf[9]
    b4  = xmf[5]
    b5  = xmf[6]
    b6  = xmf[7] + xmf[10] + xmf[11]
    b7  = xmf[10] + xmf[11] + math.log(1.0 + xmf[10])
    b12 = siti * fmm
    b13 = (siti + xmf[2] + xmf[9]) * (nak + xmf[10])
    b14 = xmf[2] * nak

    c1  = xmf[0]
    c2  = tial
    c3  = fmm
    c4  = xmf[6]
    c5  = nak
    c6  = math.log(1.0 + xmf[10] + xmf[11])
    c11 = xmf[2] + fmm + xmf[6] - xmf[9]
    c11 = c11 * (nak + xmf[10] + xmf[11])

    bcf = [b1, b2, b3, b4, b5, b6, b7, b12, b13, b14]
    ccf = [c1, c2, c3, c4, c5, c6, c11]

    BT = sum(b * c for b, c in zip(bb, bcf))
    CT = sum(c * d for c, d in zip(cc, ccf))

    TG = BT / (12.0 - AT) + CT
    F  = BT / (TG * (1.0 - CT / TG) * (1.0 - CT / TG))

    return AT, BT, CT, TG, F


def viscosity_grd(wt, T_C):
    """
    GRD08 viscosity (Pa·s) at temperature T_C (°C).
    wt: wt% dict for GRD (FeO as FeO(T), no Fe2O3 split).
    """
    A, B, C, _, _ = grd_params_from_wt(wt)
    T_K = T_C + 273.15
    log_eta = A + B / (T_K - C)
    return 10.0 ** log_eta


# ---------------------------------------------------------
# Mueller et al. 2011 relative viscosity
# ---------------------------------------------------------

def phi_m_from_rp_mueller(rp, phi_m1=0.55, b=1.0):
    """
    Mueller et al. (2011) Eq. 4
    φm = φm1 * exp( - (log10(rp))^2 / (2 b^2) )
    """
    if rp <= 0:
        return None
    return phi_m1 * math.exp(-(math.log10(rp)**2) / (2 * b**2))


def relative_viscosity_mueller(phi_c, phi_m):
    """
    Mueller et al. (2011) Eq. 3 (Maron–Pierce form)
    η_rel = (1 - φc/φm)^(-2)
    """
    if phi_m is None or phi_m <= 0:
        return None
    if phi_c >= phi_m:
        return float("inf")
    return (1.0 - phi_c / phi_m) ** -2


ASPECT_RATIOS = {
    "plag": 4.0,
    "plagioclase": 4.0,
    "cpx": 1.5,
    "clinopyroxene": 1.5,
    "ol": 1.5,
    "olivine": 1.5,
}


def get_aspect_ratio(phase_name):
    name = phase_name.lower()
    for key, ar in ASPECT_RATIOS.items():
        if key in name:
            return ar
    return 1.0


# ---------------------------------------------------------
# Streamlit UI
# ---------------------------------------------------------
st.title("MELTS Thermodynamic + Viscosity Calculator (rhyolite-MELTS 1.2.0)")

mode = st.radio(
    "Choose calculation mode:",
    ["Liquidus Temperature", "Sub-liquidus Equilibrium", "Liquid Viscosity (GRD)"]
)

st.header("Bulk Composition (wt%)")

SiO2 = st.number_input("SiO2", value=48.59)
TiO2 = st.number_input("TiO2", value=1.09)
Al2O3 = st.number_input("Al2O3", value=16.25)
FeOtot = st.number_input("FeO*", value=10.40)
MnO = st.number_input("MnO", value=0.16)
MgO = st.number_input("MgO", value=8.35)
CaO = st.number_input("CaO", value=12.91)
Na2O = st.number_input("Na2O", value=1.89)
K2O = st.number_input("K2O", value=0.26)
P2O5 = st.number_input("P2O5", value=0.12)
H2O = st.number_input("H2O", value=0.25)

Fe3FeT = st.slider("Fe³⁺/Feᵗ", 0.0, 0.5, 0.15)
pressure = st.number_input("Pressure (bars)", value=100)
buffer = st.selectbox("fO₂ buffer", ["QFM", "NNO", "IW"])

# Mode-specific inputs
if mode == "Sub-liquidus Equilibrium":
    T_sub = st.number_input("Equilibrium Temperature (°C)", value=1100)

    relvisc_on = st.checkbox("Compute crystal-bearing viscosity (Mueller et al. 2011)")

    if relvisc_on:
        liq_source = st.selectbox(
            "Liquid viscosity source:",
            ["GRD08 from bulk", "User VFT parameters"]
        )

        A_user = B_user = C_user = None
        if liq_source == "User VFT parameters":
            A_user = st.number_input("A (dimensionless)", value=-4.55)
            B_user = st.number_input("B (K)", value=10000.0)
            C_user = st.number_input("C (K)", value=500.0)

        use_override = st.checkbox("Override φm1 and b (Mueller Eq. 4)?")
        if use_override:
            phi_m1_user = st.number_input("φm1 (default 0.55)", value=0.55)
            b_user = st.number_input("b (default 1.0)", value=1.0)
        else:
            phi_m1_user = 0.55
            b_user = 1.0

if mode == "Liquid Viscosity (GRD)":
    temp_option = st.selectbox(
        "Temperature option for GRD viscosity:",
        ["Liquidus temperature", "User-defined temperature", "Temperature range", "Liquidus + user T"]
    )
    T_visc_user = st.number_input("User temperature (°C)", value=1200)
    T_range_min = st.number_input("Range min T (°C)", value=800)
    T_range_max = st.number_input("Range max T (°C)", value=1400)
    T_range_step = st.number_input("Range step (°C)", value=50)

run = st.button("Run Calculation")

# ---------------------------------------------------------
# RUN
# ---------------------------------------------------------
if run:

    # -----------------------------
    # Build bulk composition for MELTS
    # -----------------------------
    bulk = {
        'SiO2': SiO2,
        'TiO2': TiO2,
        'Al2O3': Al2O3,
        'FeO': FeOtot,
        'MnO': MnO,
        'MgO': MgO,
        'CaO': CaO,
        'Na2O': Na2O,
        'K2O': K2O,
        'P2O5': P2O5,
        'H2O': H2O
    }

    # Separate GRD bulk with FeO(T) (no Fe2O3 split)
    grd_bulk = {
        'SiO2': SiO2,
        'TiO2': TiO2,
        'Al2O3': Al2O3,
        'FeO': FeOtot,   # FeO(T)
        'MnO': MnO,
        'MgO': MgO,
        'CaO': CaO,
        'Na2O': Na2O,
        'K2O': K2O,
        'P2O5': P2O5,
        'H2O': H2O,
        'F2O-1': 0.0
    }

    # -----------------------------
    # Normalize Fe redox for MELTS
    # -----------------------------
    Fe_total = bulk['FeO']
    Fe2O3_wt = Fe_total * Fe3FeT * (159.69 / (2.0 * 71.844))
    FeO_wt   = Fe_total - Fe2O3_wt

    bulk['FeO']   = FeO_wt
    bulk['Fe2O3'] = Fe2O3_wt

    # -----------------------------
    # Initialize MELTS
    # -----------------------------
    melts = MELTSdynamic(4)

    for ox, val in bulk.items():
        melts.engine.setBulkComposition(ox.lower(), val)

    melts.engine.pressure = pressure
    melts.engine.setSystemProperties(
        "log fO2 path", buffer,
        "log fO2 offset", "0.0"
    )

    # =========================================================
    # LIQUIDUS SEARCH
    # =========================================================
    T = 1650
    melts.engine.temperature = T
    melts.engine.calcEquilibriumState()

    step = 25
    node = melts

    # Coarse search
    while True:
        if node.engine.status.failed:
            st.error(f"MELTS failed at T = {T} °C during coarse search.")
            st.stop()

        if len(node.engine.solidNames) > 0:
            break

        T -= step
        node.engine.temperature = T
        node.engine.calcEquilibriumState()

    # Fine search
    T = T + step
    step = 2
    node.engine.temperature = T
    node.engine.calcEquilibriumState()

    while True:
        if node.engine.status.failed:
            st.error(f"MELTS failed at T = {T} °C during fine search.")
            st.stop()

        if len(node.engine.solidNames) > 0:
            liquidus = T
            break

        T -= step
        node.engine.temperature = T
        node.engine.calcEquilibriumState()

    # =========================================================
    # MODE 1 — LIQUIDUS
    # =========================================================
    if mode == "Liquidus Temperature":
        st.success(f"Liquidus temperature (fO₂‑constrained) = {liquidus:.1f} °C")
        st.stop()

    # =========================================================
    # MODE 3 — LIQUID VISCOSITY (GRD)
    # =========================================================
    if mode == "Liquid Viscosity (GRD)":
        A, B, C, Tg, F = grd_params_from_wt(grd_bulk)

        st.subheader("GRD08 VFT Parameters")
        st.write(f"A = {A:.4f}")
        st.write(f"B = {B:.2f} K")
        st.write(f"C = {C:.2f} K")
        st.write(f"Tg = {Tg:.2f} K")
        st.write(f"F  = {F:.2f}")

        st.subheader("Viscosity Results (log₁₀ η, Pa·s)")

        # Liquidus temperature
        if temp_option in ["Liquidus temperature", "Liquidus + user T"]:
            eta_liq = viscosity_grd(grd_bulk, liquidus)
            st.write(f"At liquidus ({liquidus:.1f} °C): log₁₀ η = {math.log10(eta_liq):.3f}")

        # User-defined temperature
        if temp_option in ["User-defined temperature", "Liquidus + user T"]:
            eta_user = viscosity_grd(grd_bulk, T_visc_user)
            st.write(f"At {T_visc_user:.1f} °C: log₁₀ η = {math.log10(eta_user):.3f}")

        # Temperature range
        if temp_option == "Temperature range":
            if T_range_step <= 0 or T_range_max <= T_range_min:
                st.error("Check range and step: need T_max > T_min and step > 0.")
                st.stop()

            temps = list(range(int(T_range_min), int(T_range_max) + 1, int(T_range_step)))
            rows = []
            for T_C in temps:
                eta = viscosity_grd(grd_bulk, T_C)
                rows.append({
                    "T_C": T_C,
                    "log10_eta": math.log10(eta)
                })

            df = pd.DataFrame(rows)
            st.write("log₁₀(η) [Pa·s] vs T (°C):")
            st.dataframe(df)

            st.line_chart(df.set_index("T_C")["log10_eta"])
            st.caption("log₁₀(η) [Pa·s] as a function of temperature (°C).")

        st.stop()

    # =========================================================
    # MODE 2 — SUB-LIQUIDUS EQUILIBRIUM
    # =========================================================

    # 1. Prevent running above liquidus
    if T_sub >= liquidus:
        st.warning(
            f"Temperature {T_sub} °C is at or above the liquidus ({liquidus:.1f} °C). "
            "System is fully molten."
        )
        st.write("Melt fraction: 1.000")
        st.write("Solid fraction: 0.000")
        st.stop()

    # 2. Handle MELTS near-liquidus discontinuity
    deltaT = liquidus - T_sub
    if deltaT < 8:
        st.info(
            "Temperature is very close to the liquidus. "
            "Using smoothed melt fraction to avoid MELTS discontinuity."
        )
        melt_fraction = max(0.0, deltaT / 8.0)
        solid_fraction = 1.0 - melt_fraction

        st.subheader("Equilibrium Results (Smoothed)")
        st.write(f"Temperature: {T_sub} °C")
        st.write(f"Melt fraction: {melt_fraction:.3f}")
        st.write(f"Solid fraction: {solid_fraction:.3f}")
        st.write("Solids: Incipient crystallization (no stable phase assemblage).")
        st.stop()

    # 3. True MELTS equilibrium
    melts.engine.temperature = T_sub
    melts.engine.calcEquilibriumState()

    if melts.engine.status.failed:
        st.error("MELTS failed to converge at this temperature.")
        st.stop()

    liquids = list(melts.engine.liquidNames)
    solids = list(melts.engine.solidNames)

    phase_masses = {}
    for ph in liquids + solids:
        m = melts.engine.mass.get(ph, 0.0)
        if m > 1e-9:
            phase_masses[ph] = m

    total_mass = sum(phase_masses.values())
    melt_mass = sum(phase_masses[ph] for ph in liquids)
    melt_fraction = melt_mass / total_mass if total_mass > 0 else 0.0
    solid_fraction = 1.0 - melt_fraction

    solid_wtperc = {
        ph: 100.0 * phase_masses[ph] / total_mass
        for ph in solids if ph in phase_masses
    }
    solid_wtperc = dict(sorted(solid_wtperc.items(), key=lambda x: -x[1]))

    st.subheader("Equilibrium Results")
    st.write(f"Temperature: {T_sub} °C")
    st.write(f"Melt fraction: {melt_fraction:.3f}")
    st.write(f"Solid fraction: {solid_fraction:.3f}")

    st.subheader("Solid Phases (wt%)")
    if len(solid_wtperc) == 0:
        st.write("No solids present.")
    else:
        for ph, wt in solid_wtperc.items():
            st.write(f"- {ph}: {wt:.2f} wt%")

    # -----------------------------------------------------
    # Mueller et al. (2011) Relative + Suspension Viscosity
    # -----------------------------------------------------
    if 'relvisc_on' in locals() and relvisc_on:

        # 1. Solid mass fractions (system and normalized within solids)
        solid_mass_frac_sys = {ph: wt / 100.0 for ph, wt in solid_wtperc.items()}
        sum_solid = sum(solid_mass_frac_sys.values())
        if sum_solid > 0:
            solid_mass_frac = {ph: mf / sum_solid for ph, mf in solid_mass_frac_sys.items()}
        else:
            solid_mass_frac = {ph: 0.0 for ph in solid_mass_frac_sys}

        # 2. Weighted-average aspect ratio
        rp_avg = 0.0
        for ph, mf in solid_mass_frac.items():
            rp_avg += mf * get_aspect_ratio(ph)

        # 3. φm from averaged aspect ratio
        phi_m_mix = phi_m_from_rp_mueller(rp_avg, phi_m1=phi_m1_user, b=b_user)
        phi_c = solid_fraction

        # 4. Liquid viscosity at T_sub
        T_K = T_sub + 273.15
        if liq_source == "GRD08 from bulk":
            eta_liq = viscosity_grd(grd_bulk, T_sub)
        else:
            log_eta_liq = A_user + B_user / (T_K - C_user)
            eta_liq = 10.0 ** log_eta_liq

        # 5. Relative viscosity and suspension viscosity
        eta_rel = relative_viscosity_mueller(phi_c, phi_m_mix)

        st.subheader("Mueller et al. (2011) Crystal-bearing Viscosity")
        st.write(f"Crystal fraction φc (from MELTS): {phi_c:.3f}")
        st.write(f"Averaged aspect ratio rp_avg: {rp_avg:.3f}")
        st.write(f"Maximum packing φm (from rp_avg): {phi_m_mix:.3f}")
        st.write(f"log₁₀ η_rel: {math.log10(eta_rel) if math.isfinite(eta_rel) else '∞'}")

        if math.isfinite(eta_rel):
            eta_susp = eta_rel * eta_liq
            st.write(f"log₁₀ η_liq [Pa·s]: {math.log10(eta_liq):.3f}")
            st.write(f"log₁₀ η_susp [Pa·s]: {math.log10(eta_susp):.3f}")
        else:
            st.warning("φc ≥ φm → viscosity diverges (jammed suspension).")

        # 6. Phase table with AR (normalized within solids)
        rows_phase = []
        for ph, mf in solid_mass_frac.items():
            rp = get_aspect_ratio(ph)
            rows_phase.append({
                "Phase": ph,
                "Mass_frac_in_solids": mf,
                "Aspect_ratio": rp
            })

        st.subheader("Phase Aspect Ratios (normalized within solids)")
        df_phase = pd.DataFrame(rows_phase)
        st.dataframe(df_phase)

        # 7. Plot: log₁₀ η_rel vs φc
        st.subheader("Relative Viscosity Curve (Mueller Eq. 3)")
        if phi_m_mix is not None and phi_m_mix > 0:
            phi_vals = [i / 100.0 * phi_m_mix for i in range(0, 100)]
            eta_rel_vals = [relative_viscosity_mueller(phi, phi_m_mix) for phi in phi_vals]
            df_rel = pd.DataFrame({
                "phi_c": phi_vals,
                "log10_eta_rel": [
                    math.log10(v) if (v is not None and math.isfinite(v) and v > 0) else None
                    for v in eta_rel_vals
                ]
            })
            st.line_chart(df_rel.set_index("phi_c")["log10_eta_rel"])
            st.caption("log₁₀(η_rel) vs φc. φc is the crystal fraction; divergence as φc → φm.")

        # 8. Plot: φm vs aspect ratio
        st.subheader("φm vs Aspect Ratio (Mueller Eq. 4)")
        rp_vals = [0.5 + 0.1 * i for i in range(0, 51)]
        phi_m_vals = [phi_m_from_rp_mueller(rp, phi_m1=phi_m1_user, b=b_user) for rp in rp_vals]
        df_phi = pd.DataFrame({
            "rp": rp_vals,
            "phi_m": phi_m_vals
        })
        st.line_chart(df_phi.set_index("rp")["phi_m"])
        st.caption("φm decreases with increasing aspect ratio (more elongate crystals).")
