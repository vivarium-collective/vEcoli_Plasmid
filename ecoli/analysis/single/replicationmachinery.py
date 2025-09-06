import os
import numpy as np
from typing import Any

from duckdb import DuckDBPyConnection
import altair as alt
import pandas as pd


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

    # replisome monomer subunits
    replisome_monomer_cols = [5377, 6569, 7233, 7189]
    selected_replisome_monomer_cols = ["Time (min)"] + [
        bulk_df.columns[i] for i in replisome_monomer_cols
    ]
    replisome_monomer_df = bulk_df[selected_replisome_monomer_cols].copy()

    # Rename for clarity
    replisome_monomer_df.columns = [
        "Time (min)",
        "CPLX0-3621[c]",
        "EG10239-MONOMER[c]",
        "EG11500-MONOMER[c]",
        "EG11412-MONOMER[c]",
    ]

    plot_df = pd.melt(
        replisome_monomer_df,
        id_vars=["Time (min)"],
        value_vars=[
            "CPLX0-3621[c]",
            "EG10239-MONOMER[c]",
            "EG11500-MONOMER[c]",
            "EG11412-MONOMER[c]",
        ],
        var_name="Monomer",
        value_name="Count",
    )

    chart = (
        alt.Chart(plot_df)
        .mark_line()
        .encode(
            x=alt.X("Time (min):Q", title="Time (min)"),
            y=alt.Y("Count:Q", title="Counts"),
            color=alt.Color("Monomer:N", title="Replisome Subunit"),
        )
        .properties(title="Replisome Monomer counts Over Time", width=700, height=400)
    )

    chart.save(os.path.join(outdir, "replisome_monomers.html"))

    # replisome trimer subunits
    replisome_trimer_cols = [5305, 5383]
    selected_replisome_trimer_cols = ["Time (min)"] + [
        bulk_df.columns[i] for i in replisome_trimer_cols
    ]
    replisome_trimer_df = bulk_df[selected_replisome_trimer_cols].copy()

    # Rename for clarity
    replisome_trimer_df.columns = ["Time (min)", "CPLX0-2361[c]", "CPLX0-3761[c]"]

    plot_dftrimer = pd.melt(
        replisome_trimer_df,
        id_vars=["Time (min)"],
        value_vars=["CPLX0-2361[c]", "CPLX0-3761[c]"],
        var_name="Trimer",
        value_name="Count",
    )

    charttrimer = (
        alt.Chart(plot_dftrimer)
        .mark_line()
        .encode(
            x=alt.X("Time (min):Q", title="Time (min)"),
            y=alt.Y("Count:Q", title="Counts"),
            color=alt.Color("Trimer:N", title="Replisome Subunit"),
        )
        .properties(title="Replisome Trimer counts Over Time", width=700, height=400)
    )

    charttrimer.save(os.path.join(outdir, "replisome_trimers.html"))
