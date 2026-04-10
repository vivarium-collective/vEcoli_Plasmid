import os
from typing import Any
import numpy as np
import pandas as pd
import altair as alt

from duckdb import DuckDBPyConnection


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
    SELECT time, bulk FROM ({history_sql})
    ORDER BY time ASC;
    """

    output_df = conn.sql(query).df()
    bulk_matrix = np.stack(output_df["bulk"].values).astype(int)
    bulk_df = pd.DataFrame(bulk_matrix)
    time_minutes = output_df["time"].to_numpy() / 60
    bulk_df.insert(0, "Time (min)", time_minutes)
    bulk_df.to_csv(os.path.join(outdir, "bulk_matrix.csv"), index=False)

    dna_primase = bulk_df.columns[6614]
    dna_primase_df = bulk_df[["Time (min)", dna_primase]].copy()
    dna_primase_df.rename(columns={dna_primase: "DNA_primase"}, inplace=True)
    chart = (
        alt.Chart(dna_primase_df)
        .mark_line()
        .encode(
            x=alt.X("Time (min):Q", title="Time (min)"),
            y=alt.Y("DNA_primase:Q", title="Count"),
        )
        .properties(title="DNA Primase over Time")
    )
    chart.save(os.path.join(outdir, "dna_primase_counts.html"))
