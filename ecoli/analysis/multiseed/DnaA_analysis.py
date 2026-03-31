import os
from typing import Any

from duckdb import DuckDBPyConnection
import altair as alt

from ecoli.library.parquet_emitter import (
    read_stacked_columns,
)


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

    DnaA_box_counts_sql = read_stacked_columns(
        history_sql,
        [
            "listeners__unique_molecule_counts__DnaA_box AS DnaA_box_counts",
        ],
        order_results=False,
    )

    query = f"""
               SELECT time,lineage_seed, listeners__unique_molecule_counts__DnaA_box
               FROM ({DnaA_box_counts_sql})
               ORDER BY lineage_seed,time 
               """

    output_df = conn.sql(query).df()
    # Convert time from seconds to minutes
    output_df["Time (min)"] = output_df["time"] / 60

    # Rename for convenience
    output_df = output_df.rename(
        columns={"listeners__unique_molecule_counts__DnaA_box": "counts"}
    )

    # --- STEP 1: Compute median trajectory ---
    median_df = (
        output_df.groupby("Time (min)")["counts"]
        .median()
        .reset_index()
        .rename(columns={"counts": "median_counts"})
    )

    # --- STEP 2: Find median seed ---
    merged = output_df.merge(median_df, on="Time (min)")
    merged["diff"] = (merged["counts"] - merged["median_counts"]) ** 2

    seed_scores = merged.groupby("lineage_seed")["diff"].mean().reset_index()

    median_seed = seed_scores.sort_values("diff").iloc[0]["lineage_seed"]
    print("Median seed:", median_seed)

    # --- STEP 3: Split data ---
    background_df = output_df[output_df["lineage_seed"] != median_seed]
    highlight_df = output_df[output_df["lineage_seed"] == median_seed]

    # --- STEP 4: Plot ---
    background_chart = (
        alt.Chart(background_df)
        .mark_line(opacity=0.1, strokeWidth=1)
        .encode(
            x=alt.X("Time (min):Q", title="Time (min)"),
            y=alt.Y("counts:Q", title="DnaA box counts"),
            detail="lineage_seed:N",  # IMPORTANT
        )
    )

    highlight_chart = (
        alt.Chart(highlight_df)
        .mark_line(strokeWidth=3, color="black")
        .encode(x="Time (min):Q", y="counts:Q")
    )

    chart = (background_chart + highlight_chart).properties(
        title="DnaA Box Counts Across Multiple Seeds",
        width=600,
        height=400,
    )

    # --- SAVE ---
    output_df.to_csv(os.path.join(outdir, "DnaAboxescounts_multiseed.csv"), index=False)

    html_path = os.path.join(outdir, "DnaAboxescounts_multiseed.html")
    chart.save(html_path)

    # query2 = f"""
    #                SELECT time, listeners__replication_data__free_DnaA_boxes
    #                FROM ({history_sql})
    #                ORDER BY time ASC
    #                """
    #
    # output_df2 = conn.sql(query2).df()
    # # Convert time from seconds to minutes
    # output_df2["Time (min)"] = output_df2["time"] / 60
    #
    # # Create Altair line chart
    # chart2 = (
    #     alt.Chart(output_df2)
    #     .mark_line()
    #     .encode(
    #         x=alt.X("Time (min):Q", title="Time (min)"),
    #         y=alt.Y(
    #             "listeners__replication_data__free_DnaA_boxes:Q",
    #             title="Free DnaA box counts",
    #         ),
    #     )
    #     .properties(
    #         title="Count of free DnaA boxes Over Time",
    #         width=600,
    #         height=400,
    #     )
    # )
    # output_df2.to_csv(os.path.join(outdir, "freeDnaAboxescounts.csv"), index=False)
    # html_path2 = os.path.join(outdir, "freeDnaAboxescounts.html")
    # chart2.save(html_path2)
    #
    # query3 = f"""
    #                    SELECT time, listeners__replication_data__total_DnaA_boxes
    #                    FROM ({history_sql})
    #                    ORDER BY time ASC
    #                    """
    #
    # output_df3 = conn.sql(query3).df()
    # # Convert time from seconds to minutes
    # output_df3["Time (min)"] = output_df3["time"] / 60
    #
    # # Create Altair line chart
    # chart3 = (
    #     alt.Chart(output_df3)
    #     .mark_line()
    #     .encode(
    #         x=alt.X("Time (min):Q", title="Time (min)"),
    #         y=alt.Y(
    #             "listeners__replication_data__total_DnaA_boxes:Q",
    #             title="Total DnaA box counts",
    #         ),
    #     )
    #     .properties(
    #         title="Count of total DnaA boxes Over Time",
    #         width=600,
    #         height=400,
    #     )
    # )
    # output_df3.to_csv(os.path.join(outdir, "totalDnaAboxescounts.csv"), index=False)
    # html_path3 = os.path.join(outdir, "totalDnaAboxescounts.html")
    # chart3.save(html_path3)
    #
    # query6 = f"""
    #         SELECT time, bulk
    #         FROM ({history_sql})
    #         ORDER BY time ASC
    #         """
    # output_df6 = conn.sql(query6).df()
    # time_minutes = output_df6["time"].to_numpy() / 60
    # bulk_matrix = np.stack(output_df6["bulk"].values).astype(int)
    # bulk_df = pd.DataFrame(bulk_matrix)
    # bulk_df.insert(0, "Time (min)", time_minutes)
    # bulk_df.to_csv(os.path.join(outdir, "bulk_matrix.csv"), index=False)
    #
    # # DnaA proteins
    # DnaA_cols = [11524, 10781]
    # # Create a DataFrame with just Time and DnaAproteins
    # selected_cols = ["Time (min)"] + [bulk_df.columns[i] for i in DnaA_cols]
    # DnaAprotein_df = bulk_df[selected_cols].copy()
    # DnaAprotein_df.columns = ["Time (min)", "DnaA", "DnaA-ATP"]
    #
    # # Melt for plotting
    # melted_DnaA_df = pd.melt(
    #     DnaAprotein_df, id_vars=["Time (min)"], var_name="Form", value_name="Counts"
    # )
    #
    # # Plot
    # chart6 = (
    #     alt.Chart(melted_DnaA_df)
    #     .mark_line()
    #     .encode(
    #         x=alt.X("Time (min):Q", title="Time (min)"),
    #         y=alt.Y("Counts:Q", title="DnaA (counts)"),
    #         color=alt.Color("Form:N", scale=alt.Scale(range=COLORS)),
    #     )
    #     .properties(title="DnaA Counts Over Time")
    # )
    #
    # # Save or display
    # chart6.save(os.path.join(outdir, "DnaA_counts.html"))
    #
    # # Full chromosome counts
    # query7 = f"""
    #                            SELECT time, listeners__unique_molecule_counts__full_chromosome
    #                            FROM ({history_sql})
    #                            ORDER BY time ASC
    #                            """
    #
    # output_df7 = conn.sql(query7).df()
    # # Convert time from seconds to minutes
    # output_df7["Time (min)"] = output_df7["time"] / 60
    #
    # # Create Altair line chart
    # chart7 = (
    #     alt.Chart(output_df7)
    #     .mark_line()
    #     .encode(
    #         x=alt.X("Time (min):Q", title="Time (min)"),
    #         y=alt.Y(
    #             "listeners__unique_molecule_counts__full_chromosome:Q",
    #             title="Full chromosome counts",
    #         ),
    #     )
    #     .properties(
    #         title="Count of Full Chromosomes Over Time",
    #         width=600,
    #         height=400,
    #     )
    # )
    #
    # html_path7 = os.path.join(outdir, "fullchromosomecounts.html")
    # chart7.save(html_path7)
    #
    # # active replisomes
    # query8 = f"""
    #                                SELECT time, listeners__unique_molecule_counts__active_replisome
    #                                FROM ({history_sql})
    #                                ORDER BY time ASC
    #                                """
    #
    # output_df8 = conn.sql(query8).df()
    # # Convert time from seconds to minutes
    # output_df8["Time (min)"] = output_df8["time"] / 60
    #
    # # Create Altair line chart
    # chart8 = (
    #     alt.Chart(output_df8)
    #     .mark_line()
    #     .encode(
    #         x=alt.X("Time (min):Q", title="Time (min)"),
    #         y=alt.Y(
    #             "listeners__unique_molecule_counts__active_replisome:Q",
    #             title="Active replisome counts",
    #         ),
    #     )
    #     .properties(
    #         title="Count of Active Replisomes Over Time",
    #         width=600,
    #         height=400,
    #     )
    # )
    #
    # html_path8 = os.path.join(outdir, "activereplisomecounts.html")
    # chart8.save(html_path8)
    #
    # # Prepare free + total box dfs
    # free_df = output_df2[
    #     ["Time (min)", "listeners__replication_data__free_DnaA_boxes"]
    # ].copy()
    # free_df.columns = ["Time (min)", "Counts"]
    # free_df["Type"] = "Free DnaA boxes"
    #
    # total_df = output_df3[
    #     ["Time (min)", "listeners__replication_data__total_DnaA_boxes"]
    # ].copy()
    # total_df.columns = ["Time (min)", "Counts"]
    # total_df["Type"] = "Total DnaA boxes"
    #
    # # DnaA protein already melted
    # dnaA_df = melted_DnaA_df.copy()
    # dnaA_df.columns = ["Time (min)", "Type", "Counts"]
    #
    # # Combine all
    # top_df = pd.concat([dnaA_df, free_df, total_df], ignore_index=True)
    #
    # top_chart = (
    #     alt.Chart(top_df)
    #     .mark_line()
    #     .encode(
    #         x=alt.X("Time (min):Q"),
    #         y=alt.Y("Counts:Q", title="Counts"),
    #         color=alt.Color("Type:N", scale=alt.Scale(range=COLORS)),
    #     )
    #     .properties(
    #         title="DnaA Dynamics and DnaA Box Availability",
    #         width=700,
    #         height=300,
    #     )
    # )
