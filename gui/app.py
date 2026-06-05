"""SCOPE2 Sensitivity Explorer — Streamlit app.

Usage:
    cd SCOPE2/gui
    pip install -r requirements.txt
    streamlit run app.py
"""

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import streamlit as st

from utils import (
    PARAM_META, OUTPUT_META, PROFILES, PROFILE_COLORS,
    ANGLE_SCENARIOS, ANGLE_COLORS, ANGLE_DESCRIPTIONS,
    PROFILE_METRICS, PROFILE_SWEEP_VALUES,
    lut_ready, profiles_ready,
    load_sensitivity, load_sensitivity_all_profiles, load_reflectance,
    load_vertical_profiles,
)

# ── Page config ──────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="SCOPE2 Sensitivity Explorer",
    page_icon="🌿",
    layout="wide",
)

st.title("SCOPE2 Sensitivity Explorer")
st.caption(
    "Pre-computed mSCOPE simulations (10 canopy layers). "
    "All curves use default parameter values except the swept variable."
)

# ── Guard: check LUT exists ──────────────────────────────────────────────────
if not lut_ready():
    st.error(
        "**Lookup table not found.** "
        "Run `generate_lut.m` from the SCOPE2 root directory in MATLAB first, "
        "then restart this app."
    )
    st.stop()

# ── Sidebar ──────────────────────────────────────────────────────────────────
with st.sidebar:
    st.header("Layer Profile")
    profile = st.radio(
        "Canopy LAI distribution (top → bottom)",
        PROFILES,
        help=(
            "**Uniform**: equal LAI in each layer.\n\n"
            "**Top-heavy**: more foliage near the top (exponential decay downward).\n\n"
            "**Bottom-heavy**: more foliage near the bottom."
        ),
    )

    st.divider()
    st.header("Default Parameter Values")
    st.caption("Curves are computed at these values; only the selected sweep parameter varies.")
    for pname, meta in PARAM_META.items():
        st.markdown(f"**{pname}** = {meta['default']} {meta['units']}")

    st.divider()
    st.markdown(
        "**Profiles use 10 canopy layers.** "
        "Leaf biochemistry (Cab = 40 μg cm⁻², Vcmax₀ = 60 μmol m⁻² s⁻¹) "
        "is uniform across layers."
    )

# ── Tabs ─────────────────────────────────────────────────────────────────────
tab_sens, tab_refl, tab_vert = st.tabs([
    "📈 Sensitivity Curves", "🌈 Reflectance Spectra", "🌿 Vertical Profiles"
])

# ─────────────────────────────────────────────────────────────────────────────
# Tab 1 — Sensitivity curves
# ─────────────────────────────────────────────────────────────────────────────
with tab_sens:
    col_left, col_right = st.columns([1, 3])

    with col_left:
        x_param = st.selectbox(
            "X-axis (sweep parameter)",
            list(PARAM_META.keys()),
            format_func=lambda p: f"{p}  [{PARAM_META[p]['units']}]",
        )
        y_output = st.selectbox(
            "Y-axis (output)",
            list(OUTPUT_META.keys()),
            format_func=lambda o: f"{o}  [{OUTPUT_META[o]['units']}]",
        )
        overlay_all = st.checkbox("Overlay all 3 profiles", value=True)
        show_rd = (y_output == "Actot") and st.checkbox(
            "Show dark respiration (Rd) as dashed line", value=False
        )

    with col_right:
        fig, ax = plt.subplots(figsize=(8, 4.5))

        if overlay_all:
            for prof in PROFILES:
                df = load_sensitivity(x_param, prof)
                ax.plot(
                    df["value"], df[y_output],
                    color=PROFILE_COLORS[prof], linewidth=2, label=prof,
                )
                if show_rd:
                    ax.plot(
                        df["value"], -df["Rd"],
                        color=PROFILE_COLORS[prof], linewidth=1,
                        linestyle="--", alpha=0.6,
                    )
        else:
            df = load_sensitivity(x_param, profile)
            ax.plot(df["value"], df[y_output], color=PROFILE_COLORS[profile], linewidth=2)
            if show_rd:
                ax.plot(
                    df["value"], -df["Rd"],
                    color=PROFILE_COLORS[profile], linewidth=1,
                    linestyle="--", alpha=0.6, label="−Rd",
                )

        # Annotate default value
        default_val = PARAM_META[x_param]["default"]
        ax.axvline(default_val, color="grey", linestyle=":", linewidth=1, alpha=0.7,
                   label=f"Default ({default_val} {PARAM_META[x_param]['units']})")

        ax.set_xlabel(f"{x_param}  [{PARAM_META[x_param]['units']}]", fontsize=12)
        ax.set_ylabel(f"{y_output}  [{OUTPUT_META[y_output]['units']}]", fontsize=12)
        ax.set_title(
            f"{OUTPUT_META[y_output]['label']} vs {PARAM_META[x_param]['label']}",
            fontsize=13,
        )
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        st.pyplot(fig)
        plt.close(fig)

    # Current-value table
    st.subheader("Output values at all default parameters")
    rows = []
    for pname in PARAM_META:
        df_p = load_sensitivity(pname, profile)
        default_v = PARAM_META[pname]["default"]
        idx = (df_p["value"] - default_v).abs().idxmin()
        row_vals = {
            "Sweep param": f"{pname} = {default_v} {PARAM_META[pname]['units']}",
        }
        for out in OUTPUT_META:
            row_vals[f"{out} [{OUTPUT_META[out]['units']}]"] = f"{df_p.loc[idx, out]:.2f}"
        rows.append(row_vals)
    st.dataframe(rows, use_container_width=True, hide_index=True)

# ─────────────────────────────────────────────────────────────────────────────
# Tab 2 — Reflectance spectra
# ─────────────────────────────────────────────────────────────────────────────
with tab_refl:
    col_left2, col_right2 = st.columns([1, 3])

    with col_left2:
        refl_param = st.selectbox(
            "Vary parameter",
            list(PARAM_META.keys()),
            key="refl_param",
            format_func=lambda p: f"{p}  [{PARAM_META[p]['units']}]",
        )
        n_spectra = st.slider(
            "Number of spectra to show", min_value=5, max_value=40, value=10,
            help="Evenly spaced subset of the 40 simulated values.",
        )
        show_regions = st.checkbox("Shade spectral regions", value=True)

    with col_right2:
        step = max(1, 40 // n_spectra)
        wl, refl_matrix, param_values = load_reflectance(refl_param, profile, step=step)

        cmap   = cm.get_cmap("RdYlBu_r", len(param_values))
        colors = [cmap(k) for k in range(len(param_values))]

        fig2, ax2 = plt.subplots(figsize=(8, 4.5))

        # Shade spectral regions
        if show_regions:
            ax2.axvspan(400,  700, alpha=0.06, color="#5599ff", label="VIS (400–700 nm)")
            ax2.axvspan(700,  1300, alpha=0.05, color="#ff9900", label="NIR (700–1300 nm)")
            ax2.axvspan(1300, 2400, alpha=0.04, color="#aaaaaa", label="SWIR (1300–2400 nm)")

        for j, (refl_row, pval) in enumerate(zip(refl_matrix, param_values)):
            ax2.plot(wl, refl_row, color=colors[j], linewidth=0.9, alpha=0.85)

        # Colorbar
        sm = plt.cm.ScalarMappable(
            cmap="RdYlBu_r",
            norm=plt.Normalize(vmin=param_values.min(), vmax=param_values.max()),
        )
        sm.set_array([])
        cbar = fig2.colorbar(sm, ax=ax2, pad=0.02)
        cbar.set_label(f"{refl_param}  [{PARAM_META[refl_param]['units']}]", fontsize=10)

        ax2.set_xlim(400, 2400)
        ax2.set_ylim(0, None)
        ax2.set_xlabel("Wavelength  [nm]", fontsize=12)
        ax2.set_ylabel("Canopy reflectance  [—]", fontsize=12)
        ax2.set_title(
            f"Reflectance spectra as {PARAM_META[refl_param]['label']} varies  "
            f"({profile} profile)",
            fontsize=13,
        )
        if show_regions:
            ax2.legend(loc="upper right", fontsize=8)
        ax2.grid(True, alpha=0.25)
        fig2.tight_layout()
        st.pyplot(fig2)
        plt.close(fig2)

    st.caption(
        "Blue (cool) = low parameter value · Red (warm) = high parameter value. "
        "Spectra shown as bidirectional reflectance factor (BRF) computed from SCOPE's "
        "RTMo output: BRF = π · L↑ / (E_sun + E_sky)."
    )

# ─────────────────────────────────────────────────────────────────────────────
# Tab 3 — Vertical Profiles
# ─────────────────────────────────────────────────────────────────────────────
with tab_vert:
    if not profiles_ready():
        st.info(
            "**Vertical profile data not found.** "
            "Run `generate_profiles.m` from the SCOPE2 root directory in MATLAB, "
            "then refresh this page. "
            "The sensitivity and reflectance tabs above are unaffected."
        )
    else:
        col_left3, col_right3 = st.columns([1, 3])

        with col_left3:
            vert_metric = st.selectbox(
                "Metric",
                list(PROFILE_METRICS.keys()),
                format_func=lambda m: f"{m}  [{PROFILE_METRICS[m]['units']}]",
            )
            vert_sweep = st.radio(
                "Vary background parameter",
                ["Total LAI", "Incoming radiation (Rin)"],
            )
            sweep_key = "LAI" if vert_sweep.startswith("Total") else "Rin"
            sweep_vals = PROFILE_SWEEP_VALUES[sweep_key]
            sweep_labels = {
                "LAI": [f"{v} m² m⁻²" for v in sweep_vals],
                "Rin": [f"{int(v)} W m⁻²" for v in sweep_vals],
            }[sweep_key]

            # Slider index → actual value
            sweep_idx = st.select_slider(
                f"{sweep_key} value",
                options=list(range(len(sweep_vals))),
                format_func=lambda i: sweep_labels[i],
                value=sweep_vals.index(
                    min(sweep_vals, key=lambda v: abs(v - (3.0 if sweep_key == "LAI" else 600)))
                ),
            )
            selected_sweep_val = sweep_vals[sweep_idx]

            show_ps_axis = st.checkbox("Overlay sunlit fraction (Ps) as dashed lines", value=False)

            st.divider()
            st.markdown("**Leaf angle scenarios:**")
            for sc in ANGLE_SCENARIOS:
                color = ANGLE_COLORS[sc]
                st.markdown(
                    f'<span style="color:{color}">■</span> **{sc}**: '
                    f'<small>{ANGLE_DESCRIPTIONS[sc]}</small>',
                    unsafe_allow_html=True,
                )

        with col_right3:
            df_vert = load_vertical_profiles(sweep_key, selected_sweep_val)
            actual_val = df_vert["sweep_value"].iloc[0]

            metric_meta = PROFILE_METRICS[vert_metric]

            fig3, ax3 = plt.subplots(figsize=(7, 5))

            ax3_twin = ax3.twiny() if show_ps_axis else None

            for sc in ANGLE_SCENARIOS:
                sub = df_vert[df_vert["scenario"] == sc].sort_values("depth_rel")
                color = ANGLE_COLORS[sc]
                ax3.plot(
                    sub[vert_metric], sub["depth_rel"],
                    color=color, linewidth=2, label=sc,
                )
                if show_ps_axis and ax3_twin is not None:
                    ax3_twin.plot(
                        sub["Ps"], sub["depth_rel"],
                        color=color, linewidth=1, linestyle="--", alpha=0.45,
                    )

            # Y-axis: 0 at top, 1 at bottom
            ax3.set_ylim(1, 0)
            ax3.set_xlabel(f"{vert_metric}  [{metric_meta['units']}]", fontsize=12)
            ax3.set_ylabel("Relative depth in canopy  [0 = top, 1 = bottom]", fontsize=11)
            sweep_unit = "m² m⁻²" if sweep_key == "LAI" else "W m⁻²"
            ax3.set_title(
                f"{metric_meta['label']} profile  —  {sweep_key} = {actual_val:.1f} {sweep_unit}",
                fontsize=13,
            )
            ax3.legend(fontsize=9, loc="lower right")
            ax3.grid(True, alpha=0.3)

            if show_ps_axis and ax3_twin is not None:
                ax3_twin.set_xlabel("Sunlit fraction (Ps)  [—]", fontsize=10, color="grey")
                ax3_twin.set_xlim(0, 1)
                ax3_twin.tick_params(axis="x", colors="grey")

            fig3.tight_layout()
            st.pyplot(fig3)
            plt.close(fig3)

        st.caption(
            "Each line shows the canopy-depth profile of the selected metric under a "
            "different leaf angle distribution scenario (uniform LAI, 10 mSCOPE layers). "
            "Vertical-on-top (erectophile near sky) allows more direct light to penetrate "
            "deeper; Horizontal-on-top (planophile near sky) intercepts more light at the top."
        )
