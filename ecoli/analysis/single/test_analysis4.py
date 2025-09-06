import os
import numpy as np
from typing import Any

from duckdb import DuckDBPyConnection
import altair as alt
import pandas as pd

from scipy.constants import N_A
from vivarium.library.units import units


COLORS_256 = [  # From colorbrewer2.org, qualitative 8-class set 1
    [228, 26, 28],
    [55, 126, 184],
    [77, 175, 74],
    [152, 78, 163],
    [255, 127, 0],
    [255, 255, 51],
    [166, 86, 40],
    [247, 129, 191],
]

COLORS = ["#%02x%02x%02x" % (color[0], color[1], color[2]) for color in COLORS_256]


def plot(
    params: dict[str, Any],
    conn: DuckDBPyConnection,
    history_sql: str,
    config_sql: str,
    success_sql: str,
    sim_data_paths: dict[str, dict[int, str]],
    validation_data_paths: list[str],
    outdir: str,
    variant_metadata: dict[str, dict[int, Any]],
    variant_names: dict[str, str],
):
    query = f"""
    SELECT time,bulk FROM ({history_sql})
    ORDER BY time ASC;
    """

    output_df = conn.sql(query).df()
    time_minutes = output_df["time"].to_numpy() / 60
    bulk_matrix = np.stack(output_df["bulk"].values).astype(int)
    bulk_df = pd.DataFrame(bulk_matrix)
    bulk_df.insert(0, "Time (min)", time_minutes)

    # dntps
    dntp_cols = [6211, 6222, 6299, 12449]
    selected_dntp_cols = ["Time (min)"] + [bulk_df.columns[i] for i in dntp_cols]
    dntp_df = bulk_df[selected_dntp_cols].copy()

    # Rename for clarity
    dntp_df.columns = ["Time (min)", "DATP (C)", "DCTP (C)", "DGTP (C)", "TTP (C)"]

    # Add total dNTPs column
    dntp_df["Total dNTPs"] = dntp_df[
        ["DATP (C)", "DCTP (C)", "DGTP (C)", "TTP (C)"]
    ].sum(axis=1)

    query2 = f"""
                        SELECT time, listeners__mass__volume
                        FROM ({history_sql})
                        ORDER BY time ASC
                        """

    output_df2 = conn.sql(query2).df()

    # 2. Add volume to dNTP dataframe
    dntp_df["Volume (fL)"] = output_df2["listeners__mass__volume"].values

    # Computing concentration
    AVOGADRO = N_A / units.mol
    for col in ["DATP (C)", "DCTP (C)", "DGTP (C)", "TTP (C)", "Total dNTPs"]:
        mols = dntp_df[col].values / AVOGADRO
        volume_fl = dntp_df["Volume (fL)"].values * units.fL
        conc = (mols / volume_fl).to(units.micromolar)
        new_col_name = col + " (µM)"
        dntp_df[new_col_name] = conc.magnitude

    # Columns with concentration values ending with " (µM)"
    conc_cols = [col for col in dntp_df.columns if col.endswith(" (µM)")]

    plot_df = pd.melt(
        dntp_df,
        id_vars=["Time (min)"],
        value_vars=conc_cols,
        var_name="dNTP",
        value_name="Concentration (µM)",
    )

    chart = (
        alt.Chart(plot_df)
        .mark_line()
        .encode(
            x=alt.X("Time (min):Q", title="Time (min)"),
            y=alt.Y("Concentration (µM):Q", title="Concentration (µM)"),
            color=alt.Color("dNTP:N", title="dNTP"),
        )
        .properties(title="dNTP Concentrations Over Time", width=700, height=400)
    )
    dntp_df.to_csv(os.path.join(outdir, "dntp.csv"), index=False)
    chart.save(os.path.join(outdir, "dNTP_Conc.html"))

    # polymerized dntps
    polymerized_dntp_cols = [16008, 16017, 16026, 16188]
    selected_polymerized_dntp_cols = ["Time (min)"] + [
        bulk_df.columns[i] for i in polymerized_dntp_cols
    ]
    polymerized_dntp_df = bulk_df[selected_polymerized_dntp_cols].copy()

    # Rename for clarity
    polymerized_dntp_df.columns = [
        "Time (min)",
        "polymerized_DATP (C)",
        "polymerized_DCTP (C)",
        "polymerized_DGTP (C)",
        "polymerized_TTP (C)",
    ]

    # Add total dNTPs column
    polymerized_dntp_df["polymerized_Total dNTPs"] = polymerized_dntp_df[
        [
            "polymerized_DATP (C)",
            "polymerized_DCTP (C)",
            "polymerized_DGTP (C)",
            "polymerized_TTP (C)",
        ]
    ].sum(axis=1)

    # 2. Add volume to polymerized_dNTP dataframe
    polymerized_dntp_df["Volume (fL)"] = output_df2["listeners__mass__volume"].values

    # Computing concentration
    for col in [
        "polymerized_DATP (C)",
        "polymerized_DCTP (C)",
        "polymerized_DGTP (C)",
        "polymerized_TTP (C)",
        "polymerized_Total dNTPs",
    ]:
        mols = polymerized_dntp_df[col].values / AVOGADRO
        volume_fl = polymerized_dntp_df["Volume (fL)"].values * units.fL
        conc = (mols / volume_fl).to(units.micromolar)
        new_col_name = col + " (µM)"
        polymerized_dntp_df[new_col_name] = conc.magnitude

    # Columns with concentration values ending with " (µM)"
    polymerized_conc_cols = [
        col for col in polymerized_dntp_df.columns if col.endswith(" (µM)")
    ]

    plot_df2 = pd.melt(
        polymerized_dntp_df,
        id_vars=["Time (min)"],
        value_vars=polymerized_conc_cols,
        var_name="polymerized_dNTP",
        value_name="Concentration (µM)",
    )

    chart_polymerized_dntp = (
        alt.Chart(plot_df2)
        .mark_line()
        .encode(
            x=alt.X("Time (min):Q", title="Time (min)"),
            y=alt.Y("Concentration (µM):Q", title="Concentration (µM)"),
            color=alt.Color("polymerized_dNTP:N", title="polymerized_dNTP"),
        )
        .properties(
            title="Polymerized dNTP Concentrations Over Time", width=700, height=400
        )
    )
    polymerized_dntp_df.to_csv(
        os.path.join(outdir, "polymerized_dntp.csv"), index=False
    )
    chart_polymerized_dntp.save(os.path.join(outdir, "polymerized_dNTP_Conc.html"))

    # ntps
    ntp_cols = [797, 5979, 10001, 15592]
    selected_ntp_cols = ["Time (min)"] + [bulk_df.columns[i] for i in ntp_cols]
    ntp_df = bulk_df[selected_ntp_cols].copy()

    # Rename for clarity
    ntp_df.columns = ["Time (min)", "ATP (C)", "CTP (C)", "GTP (C)", "UTP (C)"]

    # Add total NTPs column
    ntp_df["Total NTPs"] = ntp_df[["ATP (C)", "CTP (C)", "GTP (C)", "UTP (C)"]].sum(
        axis=1
    )

    # 2. Add volume to NTP dataframe
    ntp_df["Volume (fL)"] = output_df2["listeners__mass__volume"].values

    # Computing concentration
    for col in ["ATP (C)", "CTP (C)", "GTP (C)", "UTP (C)", "Total NTPs"]:
        mols = ntp_df[col].values / AVOGADRO
        volume_fl = ntp_df["Volume (fL)"].values * units.fL
        conc = (mols / volume_fl).to(units.micromolar)
        new_col_name = col + " (µM)"
        ntp_df[new_col_name] = conc.magnitude

    # Columns with concentration values ending with " (µM)"
    ntp_conc_cols = [col for col in ntp_df.columns if col.endswith(" (µM)")]

    plot_df3 = pd.melt(
        ntp_df,
        id_vars=["Time (min)"],
        value_vars=ntp_conc_cols,
        var_name="NTP",
        value_name="Concentration (µM)",
    )

    chart_ntp = (
        alt.Chart(plot_df3)
        .mark_line()
        .encode(
            x=alt.X("Time (min):Q", title="Time (min)"),
            y=alt.Y("Concentration (µM):Q", title="Concentration (µM)"),
            color=alt.Color("NTP:N", title="NTP"),
        )
        .properties(title="NTP Concentrations Over Time", width=700, height=400)
    )
    ntp_df.to_csv(os.path.join(outdir, "ntp.csv"), index=False)
    chart_ntp.save(os.path.join(outdir, "NTP_Conc.html"))

    # polymerized_NTPs
    polymerized_ntp_cols = [15981, 15990, 16062, 16206]
    selected_polymerized_ntp_cols = ["Time (min)"] + [
        bulk_df.columns[i] for i in polymerized_ntp_cols
    ]
    polymerized_ntp_df = bulk_df[selected_polymerized_ntp_cols].copy()

    # Rename for clarity
    polymerized_ntp_df.columns = [
        "Time (min)",
        "polymerized_ATP (C)",
        "polymerized_CTP (C)",
        "polymerized_GTP (C)",
        "polymerized_UTP (C)",
    ]

    # Add total NTPs column
    polymerized_ntp_df["Total polymerized_NTPs"] = polymerized_ntp_df[
        [
            "polymerized_ATP (C)",
            "polymerized_CTP (C)",
            "polymerized_GTP (C)",
            "polymerized_UTP (C)",
        ]
    ].sum(axis=1)

    # 2. Add volume to NTP dataframe
    polymerized_ntp_df["Volume (fL)"] = output_df2["listeners__mass__volume"].values

    # Computing concentration
    for col in [
        "polymerized_ATP (C)",
        "polymerized_CTP (C)",
        "polymerized_GTP (C)",
        "polymerized_UTP (C)",
        "Total polymerized_NTPs",
    ]:
        mols = polymerized_ntp_df[col].values / AVOGADRO
        volume_fl = polymerized_ntp_df["Volume (fL)"].values * units.fL
        conc = (mols / volume_fl).to(units.micromolar)
        new_col_name = col + " (µM)"
        polymerized_ntp_df[new_col_name] = conc.magnitude

    # Columns with concentration values ending with " (µM)"
    polymerized_ntp_conc_cols = [
        col for col in polymerized_ntp_df.columns if col.endswith(" (µM)")
    ]

    plot_df4 = pd.melt(
        polymerized_ntp_df,
        id_vars=["Time (min)"],
        value_vars=polymerized_ntp_conc_cols,
        var_name="polymerized_NTP",
        value_name="Concentration (µM)",
    )

    chart_polymerized_ntp = (
        alt.Chart(plot_df4)
        .mark_line()
        .encode(
            x=alt.X("Time (min):Q", title="Time (min)"),
            y=alt.Y("Concentration (µM):Q", title="Concentration (µM)"),
            color=alt.Color("polymerized_NTP:N", title="polymerized_NTP"),
        )
        .properties(
            title="Polymerized NTP Concentrations Over Time", width=700, height=400
        )
    )

    chart_polymerized_ntp.save(os.path.join(outdir, "Polymerized_NTP_Conc.html"))

    #  RNAPs
    inactive_RNAP_cols = [691]
    selected_inactive_RNAP_cols = ["Time (min)"] + [
        bulk_df.columns[i] for i in inactive_RNAP_cols
    ]
    RNAP_df = bulk_df[selected_inactive_RNAP_cols].copy()

    # Rename for clarity
    RNAP_df.columns = ["Time (min)", "inactive_RNAP (C)"]

    query3 = f"""
               SELECT time, listeners__unique_molecule_counts__active_RNAP
               FROM ({history_sql})
               ORDER BY time ASC
               """

    output_df3 = conn.sql(query3).df()

    # 2. Add active to RNAP dataframe
    RNAP_df["active_RNAP (C)"] = output_df3[
        "listeners__unique_molecule_counts__active_RNAP"
    ].values

    # 2. Add volume to RNAP dataframe
    RNAP_df["Volume (fL)"] = output_df2["listeners__mass__volume"].values

    # Computing concentration
    for col in ["inactive_RNAP (C)", "active_RNAP (C)"]:
        mols = RNAP_df[col].values / AVOGADRO
        volume_fl = RNAP_df["Volume (fL)"].values * units.fL
        conc = (mols / volume_fl).to(units.micromolar)
        new_col_name = col + " (µM)"
        RNAP_df[new_col_name] = conc.magnitude

    # Columns with concentration values ending with " (µM)"
    RNAP_conc_cols = [col for col in RNAP_df.columns if col.endswith(" (µM)")]

    plot_df5 = pd.melt(
        RNAP_df,
        id_vars=["Time (min)"],
        value_vars=RNAP_conc_cols,
        var_name="RNAP",
        value_name="Concentration (µM)",
    )

    chart_RNAP = (
        alt.Chart(plot_df5)
        .mark_line()
        .encode(
            x=alt.X("Time (min):Q", title="Time (min)"),
            y=alt.Y("Concentration (µM):Q", title="Concentration (µM)"),
            color=alt.Color("RNAP:N", title="RNAP"),
        )
        .properties(title="RNAP Concentrations Over Time", width=700, height=400)
    )

    chart_RNAP.save(os.path.join(outdir, "RNAP_Conc.html"))

    # Computing fraction of active RNAPs
    RNAP_df["active_RNAP fraction"] = RNAP_df["active_RNAP (C)"] / (
        RNAP_df["inactive_RNAP (C)"] + RNAP_df["active_RNAP (C)"]
    )

    plot_df_frac = RNAP_df[["Time (min)", "active_RNAP fraction"]]

    chart_frac = (
        alt.Chart(plot_df_frac)
        .mark_line(color="green")
        .encode(
            x=alt.X("Time (min):Q", title="Time (min)"),
            y=alt.Y("active_RNAP fraction:Q", title="Fraction of Active RNAPs"),
        )
        .properties(title="Fraction of Active RNAPs Over Time", width=700, height=400)
    )

    # Save fraction plot
    chart_frac.save(os.path.join(outdir, "RNAP_ActiveFraction.html"))

    # Rnases H
