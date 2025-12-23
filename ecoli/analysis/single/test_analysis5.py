import os
from typing import Any

from duckdb import DuckDBPyConnection
import altair as alt


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
    query3 = f"""
               SELECT time, listeners__unique_molecule_counts__full_plasmid
               FROM ({history_sql})
               ORDER BY time ASC
               """

    output_df = conn.sql(query3).df()
    # Convert time from seconds to minutes
    output_df["Time (min)"] = output_df["time"] / 60

    # Create Altair line chart
    chart = (
        alt.Chart(output_df)
        .mark_line()
        .encode(
            x=alt.X("Time (min):Q", title="Time (min)"),
            y=alt.Y(
                "listeners__unique_molecule_counts__full_plasmid:Q",
                title="Full plasmid counts",
            ),
        )
        .properties(
            title="Full plasmid counts",
            width=600,
            height=400,
        )
    )

    html_path = os.path.join(outdir, "full plasmid counts.html")
    chart.save(html_path)

    # first 20 seconds
    output_df_60s = output_df[output_df["time"] <= 60].copy()
    output_df_60s["Time (sec)"] = output_df_60s["time"]

    # Create Altair line chart
    chart_full_plasmids_60s = (
        alt.Chart(output_df_60s)
        .mark_line()
        .encode(
            x=alt.X("Time (sec):Q", title="Time (sec)"),
            y=alt.Y(
                "listeners__unique_molecule_counts__full_plasmid:Q",
                title="Full plasmid counts",
            ),
        )
        .properties(
            title="Full plasmid counts (first 60 seconds)",
            width=600,
            height=400,
        )
    )

    html_path2 = os.path.join(outdir, "full plasmid counts (60 seconds).html")
    chart_full_plasmids_60s.save(html_path2)

    # plasmid active replisome counts
    query4 = f"""
                   SELECT time, listeners__unique_molecule_counts__plasmid_active_replisome
                   FROM ({history_sql})
                   ORDER BY time ASC
                   """

    output_df2 = conn.sql(query4).df()
    # Convert time from seconds to minutes
    output_df2["Time (min)"] = output_df2["time"] / 60

    # Create Altair line chart
    chart_plasmid_active_replisome = (
        alt.Chart(output_df2)
        .mark_line()
        .encode(
            x=alt.X("Time (min):Q", title="Time (min)"),
            y=alt.Y(
                "listeners__unique_molecule_counts__plasmid_active_replisome:Q",
                title="Counts",
            ),
        )
        .properties(
            title="Plasmid active replisome counts",
            width=600,
            height=400,
        )
    )

    html_path3 = os.path.join(outdir, "plasmid_active_replisomes.html")
    chart_plasmid_active_replisome.save(html_path3)

    # first 20 seconds
    output_df2_60s = output_df2[output_df2["time"] <= 60].copy()
    output_df2_60s["Time (sec)"] = output_df2_60s["time"]

    # Create Altair line chart
    chart_plasmid_active_replisome_60s = (
        alt.Chart(output_df2_60s)
        .mark_line()
        .encode(
            x=alt.X("Time (sec):Q", title="Time (sec)"),
            y=alt.Y(
                "listeners__unique_molecule_counts__plasmid_active_replisome:Q",
                title="Counts",
            ),
        )
        .properties(
            title="Plasmid active replisome counts(first 60 seconds)",
            width=600,
            height=400,
        )
    )

    html_path4 = os.path.join(
        outdir, "plasmid active replisome counts (60 seconds).html"
    )
    chart_plasmid_active_replisome_60s.save(html_path4)
