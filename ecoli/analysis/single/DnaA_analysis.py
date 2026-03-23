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
    #  RNAPs

    query = f"""
               SELECT time, listeners__unique_molecule_counts__DnaA_box
               FROM ({history_sql})
               ORDER BY time ASC
               """

    output_df = conn.sql(query).df()
    # Convert time from seconds to minutes
    output_df["Time (min)"] = output_df["time"] / 60

    # Create Altair line chart
    chart = (
        alt.Chart(output_df)
        .mark_line()
        .encode(
            x=alt.X("Time (min):Q", title="Time (min)"),
            y=alt.Y(
                "listeners__unique_molecule_counts__DnaA_box:Q",
                title="DnaA box counts",
            ),
        )
        .properties(
            title="Count of DnaA boxes Over Time",
            width=600,
            height=400,
        )
    )
    output_df.to_csv(os.path.join(outdir, "DnaAboxescounts.csv"), index=False)
    html_path = os.path.join(outdir, "DnaAboxescounts.html")
    chart.save(html_path)

    query2 = f"""
                   SELECT time, listeners__replication_data__free_DnaA_boxes
                   FROM ({history_sql})
                   ORDER BY time ASC
                   """

    output_df2 = conn.sql(query2).df()
    # Convert time from seconds to minutes
    output_df2["Time (min)"] = output_df2["time"] / 60

    # Create Altair line chart
    chart2 = (
        alt.Chart(output_df2)
        .mark_line()
        .encode(
            x=alt.X("Time (min):Q", title="Time (min)"),
            y=alt.Y(
                "listeners__replication_data__free_DnaA_boxes:Q",
                title="Free DnaA box counts",
            ),
        )
        .properties(
            title="Count of free DnaA boxes Over Time",
            width=600,
            height=400,
        )
    )
    output_df2.to_csv(os.path.join(outdir, "freeDnaAboxescounts.csv"), index=False)
    html_path2 = os.path.join(outdir, "freeDnaAboxescounts.html")
    chart2.save(html_path2)

    query3 = f"""
                       SELECT time, listeners__replication_data__total_DnaA_boxes
                       FROM ({history_sql})
                       ORDER BY time ASC
                       """

    output_df3 = conn.sql(query3).df()
    # Convert time from seconds to minutes
    output_df3["Time (min)"] = output_df3["time"] / 60

    # Create Altair line chart
    chart3 = (
        alt.Chart(output_df3)
        .mark_line()
        .encode(
            x=alt.X("Time (min):Q", title="Time (min)"),
            y=alt.Y(
                "listeners__replication_data__total_DnaA_boxes:Q",
                title="Total DnaA box counts",
            ),
        )
        .properties(
            title="Count of total DnaA boxes Over Time",
            width=600,
            height=400,
        )
    )
    output_df3.to_csv(os.path.join(outdir, "totalDnaAboxescounts.csv"), index=False)
    html_path3 = os.path.join(outdir, "totalDnaAboxescounts.html")
    chart3.save(html_path3)

    query4 = f"""
                           SELECT time, listeners__unique_molecule_counts__oriC
                           FROM ({history_sql})
                           ORDER BY time ASC
                           """

    output_df4 = conn.sql(query4).df()
    # Convert time from seconds to minutes
    output_df4["Time (min)"] = output_df4["time"] / 60

    # Create Altair line chart
    chart4 = (
        alt.Chart(output_df4)
        .mark_line()
        .encode(
            x=alt.X("Time (min):Q", title="Time (min)"),
            y=alt.Y(
                "listeners__unique_molecule_counts__oriC:Q",
                title="oriC counts",
            ),
        )
        .properties(
            title="Count of oriC Over Time",
            width=600,
            height=400,
        )
    )

    html_path4 = os.path.join(outdir, "oriCcounts.html")
    chart4.save(html_path4)

    # DnaA mrna cistron
    query5 = f"""
                               SELECT time, listeners__rna_counts__mRNA_cistron_counts
                               FROM ({history_sql})
                               ORDER BY time ASC
                               """

    output_df5 = conn.sql(query5).df()
    # Convert time from seconds to minutes
    time_minutes = output_df5["time"] / 60
    mrna_cistron_data = np.vstack(
        output_df5["listeners__rna_counts__mRNA_cistron_counts"].values
    ).astype(int)
    mrna_cistron_df = pd.DataFrame(mrna_cistron_data)
    mrna_cistron_df.insert(0, "Time (min)", time_minutes)
    mrna_cistron_df.to_csv(os.path.join(outdir, "mrna_cistron.csv"), index=False)

    dnaA_rna_id = mrna_cistron_df.columns[228]
    dnaA_rna_df = mrna_cistron_df[["Time (min)", dnaA_rna_id]].copy()
    dnaA_rna_df.rename(columns={dnaA_rna_id: "DnaA_cistron_mRNA"}, inplace=True)

    chart5 = (
        alt.Chart(dnaA_rna_df)
        .mark_line()
        .encode(
            x=alt.X("Time (min):Q", title="Time (min)"),
            y=alt.Y("DnaA_cistron_mRNA:Q", title="DnaA_cistron_mRNA (counts)"),
        )
        .properties(title="DnaA_cistron_mRNA counts over Time")
    )
    chart5.save(os.path.join(outdir, "DnaA_cistron_mRNA.html"))

    query6 = f"""
            SELECT time, bulk
            FROM ({history_sql})
            ORDER BY time ASC
            """
    output_df6 = conn.sql(query6).df()
    time_minutes = output_df6["time"].to_numpy() / 60
    bulk_matrix = np.stack(output_df6["bulk"].values).astype(int)
    bulk_df = pd.DataFrame(bulk_matrix)
    bulk_df.insert(0, "Time (min)", time_minutes)
    bulk_df.to_csv(os.path.join(outdir, "bulk_matrix.csv"), index=False)

    # DnaA proteins
    DnaA_cols = [11524, 10781]
    # Create a DataFrame with just Time and DnaAproteins
    selected_cols = ["Time (min)"] + [bulk_df.columns[i] for i in DnaA_cols]
    DnaAprotein_df = bulk_df[selected_cols].copy()
    DnaAprotein_df.columns = ["Time (min)", "DnaA", "DnaA-ATP"]

    # Melt for plotting
    melted_DnaA_df = pd.melt(
        DnaAprotein_df, id_vars=["Time (min)"], var_name="Form", value_name="Counts"
    )

    # Plot
    chart6 = (
        alt.Chart(melted_DnaA_df)
        .mark_line()
        .encode(
            x=alt.X("Time (min):Q", title="Time (min)"),
            y=alt.Y("Counts:Q", title="DnaA (counts)"),
            color=alt.Color("Form:N", scale=alt.Scale(range=COLORS)),
        )
        .properties(title="DnaA Counts Over Time")
    )

    # Save or display
    chart6.save(os.path.join(outdir, "DnaA_counts.html"))

    # Full chromosome counts
    query7 = f"""
                               SELECT time, listeners__unique_molecule_counts__full_chromosome
                               FROM ({history_sql})
                               ORDER BY time ASC
                               """

    output_df7 = conn.sql(query7).df()
    # Convert time from seconds to minutes
    output_df7["Time (min)"] = output_df7["time"] / 60

    # Create Altair line chart
    chart7 = (
        alt.Chart(output_df7)
        .mark_line()
        .encode(
            x=alt.X("Time (min):Q", title="Time (min)"),
            y=alt.Y(
                "listeners__unique_molecule_counts__full_chromosome:Q",
                title="Full chromosome counts",
            ),
        )
        .properties(
            title="Count of Full Chromosomes Over Time",
            width=600,
            height=400,
        )
    )

    html_path7 = os.path.join(outdir, "fullchromosomecounts.html")
    chart7.save(html_path7)

    # active replisomes
    query8 = f"""
                                   SELECT time, listeners__unique_molecule_counts__active_replisome
                                   FROM ({history_sql})
                                   ORDER BY time ASC
                                   """

    output_df8 = conn.sql(query8).df()
    # Convert time from seconds to minutes
    output_df8["Time (min)"] = output_df8["time"] / 60

    # Create Altair line chart
    chart8 = (
        alt.Chart(output_df8)
        .mark_line()
        .encode(
            x=alt.X("Time (min):Q", title="Time (min)"),
            y=alt.Y(
                "listeners__unique_molecule_counts__active_replisome:Q",
                title="Active replisome counts",
            ),
        )
        .properties(
            title="Count of Active Replisomes Over Time",
            width=600,
            height=400,
        )
    )

    html_path8 = os.path.join(outdir, "activereplisomecounts.html")
    chart8.save(html_path8)

    # Prepare free + total box dfs
    free_df = output_df2[
        ["Time (min)", "listeners__replication_data__free_DnaA_boxes"]
    ].copy()
    free_df.columns = ["Time (min)", "Counts"]
    free_df["Type"] = "Free DnaA boxes"

    total_df = output_df3[
        ["Time (min)", "listeners__replication_data__total_DnaA_boxes"]
    ].copy()
    total_df.columns = ["Time (min)", "Counts"]
    total_df["Type"] = "Total DnaA boxes"

    # DnaA protein already melted
    dnaA_df = melted_DnaA_df.copy()
    dnaA_df.columns = ["Time (min)", "Type", "Counts"]

    # Combine all
    top_df = pd.concat([dnaA_df, free_df, total_df], ignore_index=True)

    top_chart = (
        alt.Chart(top_df)
        .mark_line()
        .encode(
            x=alt.X("Time (min):Q"),
            y=alt.Y("Counts:Q", title="Counts"),
            color=alt.Color("Type:N", scale=alt.Scale(range=COLORS)),
        )
        .properties(
            title="DnaA Dynamics and DnaA Box Availability",
            width=700,
            height=300,
        )
    )

    ori_df = output_df4[
        ["Time (min)", "listeners__unique_molecule_counts__oriC"]
    ].copy()
    ori_df.columns = ["Time (min)", "Counts"]
    ori_df["Type"] = "oriC"

    chrom_df = output_df7[
        ["Time (min)", "listeners__unique_molecule_counts__full_chromosome"]
    ].copy()
    chrom_df.columns = ["Time (min)", "Counts"]
    chrom_df["Type"] = "Full chromosomes"

    repl_df = output_df8[
        ["Time (min)", "listeners__unique_molecule_counts__active_replisome"]
    ].copy()
    repl_df.columns = ["Time (min)", "Counts"]
    repl_df["Type"] = "Active replisomes"

    bottom_df = pd.concat([ori_df, chrom_df, repl_df], ignore_index=True)

    bottom_chart = (
        alt.Chart(bottom_df)
        .mark_line()
        .encode(
            x=alt.X("Time (min):Q"),
            y=alt.Y("Counts:Q", title="Counts"),
            color=alt.Color("Type:N", scale=alt.Scale(range=COLORS)),
            strokeDash=alt.StrokeDash(
                "Type:N",
                scale=alt.Scale(
                    domain=["Full chromosomes", "oriC", "Active replisomes"],
                    range=[
                        [1, 0],  # chromosomes → solid
                        [4, 2],  # oriC → dashed
                        [1, 2],  # replisomes → dotted
                    ],
                ),
                legend=None,
            ),
        )
        .properties(
            title="Replication State Dynamics",
            width=700,
            height=300,
        )
    )

    final_chart = (
        alt.vconcat(top_chart, bottom_chart)
        .configure_axis(labelFontSize=14, titleFontSize=16)
        .configure_legend(titleFontSize=14, labelFontSize=12)
        .configure_title(fontSize=18, anchor="start")
    )

    final_chart.save(os.path.join(outdir, "combined_replication_figure.html"))
