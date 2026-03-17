import os
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

    # first 60 seconds
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

    # first 60 seconds
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

    # plasmid active replisomes around chromosome initiation
    # Chromosome initiation time
    chromosome_initiation_s = 1758
    chromosome_initiation_min = round(chromosome_initiation_s / 60, 2)  # 29.27 min

    # Window ±30 s
    window = 30
    start_time = chromosome_initiation_s - window
    end_time = chromosome_initiation_s + window

    # Convert to minutes for plotting
    start_min = start_time / 60
    end_min = end_time / 60

    # Filter the dataframe
    output_df2_comp = output_df2[
        (output_df2["time"] >= start_time) & (output_df2["time"] <= end_time)
    ].copy()
    output_df2_comp["Time (min)"] = output_df2_comp["time"] / 60

    # Define x-axis ticks: include start, end, and chromosome initiation
    x_ticks = sorted(
        {round(start_min, 2), round(end_min, 2), chromosome_initiation_min}
    )

    # Altair line chart
    chart_plasmid_active_replisome_comp = (
        alt.Chart(output_df2_comp)
        .mark_line()
        .encode(
            x=alt.X(
                "Time (min):Q",
                title="Time (min)",
                scale=alt.Scale(
                    domain=(start_min, end_min), zero=False
                ),  # force axis to start/end exactly at window
                axis=alt.Axis(
                    values=x_ticks,
                    labelExpr=f'datum.value == {chromosome_initiation_min} ? "Chromosome initiation at {chromosome_initiation_min} min" : datum.value',
                ),
            ),
            y=alt.Y(
                "listeners__unique_molecule_counts__plasmid_active_replisome:Q",
                title="Active plasmid replisomes",
            ),
        )
        .properties(
            title="Plasmid active replisome counts around chromosome initiation",
            width=600,
            height=400,
        )
    )

    # Add vertical line for chromosome initiation
    vline = (
        alt.Chart(pd.DataFrame({"x": [chromosome_initiation_min]}))
        .mark_rule(color="red", strokeDash=[5, 5])
        .encode(
            x="x:Q",
            tooltip=alt.Tooltip(["x"], format=".2f"),  # shows exact minutes
        )
    )

    # Combine line chart and vertical line
    final_chart = chart_plasmid_active_replisome_comp + vline

    # Save HTML
    html_path_comp2 = os.path.join(
        outdir, "plasmid_active_replisome_tight_window_minutes.html"
    )
    final_chart.save(html_path_comp2)

    ## Extracting 3 seconds upto initiation of chromosome and plasmid dynamics
    times_of_interest = [1754, 1755, 1756, 1757, 1758, 1759, 1760]
    df_plasmid = output_df2[output_df2["time"].isin(times_of_interest)].copy()

    query = f"""
        SELECT time, listeners__unique_molecule_counts__active_replisome
        FROM ({history_sql})
        ORDER BY time ASC
        """

    outputc_df = conn.sql(query).df()
    df_chr = outputc_df[outputc_df["time"].isin(times_of_interest)].copy()

    df_plasmid = df_plasmid.rename(
        columns={
            "listeners__unique_molecule_counts__plasmid_active_replisome": "plasmid_replisomes"
        }
    )

    df_chr = df_chr.rename(
        columns={
            "listeners__unique_molecule_counts__active_replisome": "chromosome_replisomes"
        }
    )

    df_replisome_merged = df_chr.merge(df_plasmid, on="time")
    df_replisome_merged["Time (min)"] = (df_replisome_merged["time"] / 60).round(2)

    df_replisome_merged["TimeLabel"] = (
        df_replisome_merged["Time (min)"].astype(str) + " min"
    )

    df_replisome_merged.loc[df_replisome_merged["time"] == 1758, "TimeLabel"] = (
        f"{round(1758 / 60, 2)} min – Chromosome initiated"
    )

    df_plot = df_replisome_merged.melt(
        id_vars=["time", "Time (min)", "TimeLabel"],
        value_vars=["chromosome_replisomes", "plasmid_replisomes"],
        var_name="Replication",
        value_name="Active replisomes",
    )

    df_plot["Replication"] = df_plot["Replication"].replace(
        {"chromosome_replisomes": "Chromosome", "plasmid_replisomes": "Plasmid"}
    )

    df_plot.to_csv(os.path.join(outdir, "replisome_dynamics_times.csv"), index=False)

    # Extracting 60s window from chromosome replisomes
    # plasmid active replisomes around chromosome initiation
    # Chromosome initiation time
    # Filter the dataframe
    output_dfc_60 = outputc_df[
        (outputc_df["time"] >= start_time) & (outputc_df["time"] <= end_time)
    ].copy()
    output_dfc_60["Time (min)"] = output_dfc_60["time"] / 60

    # Define x-axis ticks: include start, end, and chromosome initiation
    x_ticks = sorted(
        {round(start_min, 2), round(end_min, 2), chromosome_initiation_min}
    )

    # Altair line chart
    chart_active_replisome_comp = (
        alt.Chart(output_dfc_60)
        .mark_line()
        .encode(
            x=alt.X(
                "Time (min):Q",
                title="Time (min)",
                scale=alt.Scale(
                    domain=(start_min, end_min), zero=False
                ),  # force axis to start/end exactly at window
                axis=alt.Axis(
                    values=x_ticks
                    # labelExpr=f'datum.value == {chromosome_initiation_min} ? "Chromosome initiation at {chromosome_initiation_min} min" : datum.value'
                ),
            ),
            y=alt.Y(
                "listeners__unique_molecule_counts__active_replisome:Q",
                title="Active chromosome replisomes",
            ),
        )
        .properties(title="Replisomes active on chromosome", width=600, height=400)
    )
    html_path_comp3 = os.path.join(
        outdir, "chromosome_active_replisome_tight_window_minutes.html"
    )
    chart_active_replisome_comp.save(html_path_comp3)
