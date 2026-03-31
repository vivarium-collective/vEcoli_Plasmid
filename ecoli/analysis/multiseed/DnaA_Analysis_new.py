# %%
import os
import numpy as np
import duckdb

import altair as alt
import pandas as pd

from ecoli.library.parquet_emitter import (
    read_stacked_columns,
)

conn = duckdb.connect()
query_dict = {
    "experiment_id": "Multiseed_DnaA",
    "variant": 0,
    "generation": 1,
}
history_sql = f"""
SELECT * FROM read_parquet('out/{query_dict["experiment_id"]}/mdna/history/*/*/*/*/*/*.pq', hive_partitioning=true)
"""
DnaA_box_counts_sql = read_stacked_columns(
    history_sql,
    [
        "listeners__unique_molecule_counts__DnaA_box AS DnaA_box_counts",
    ],
    order_results=False,
)

query = f"""
SELECT lineage_seed, time, DnaA_box_counts FROM ({DnaA_box_counts_sql})
ORDER BY lineage_seed, time
"""
output_df = conn.sql(query).df()

# Convert time from seconds to minutes
output_df["Time (min)"] = output_df["time"] / 60

# Rename for convenience
output_df = output_df.rename(columns={"DnaA_box_counts": "counts"})

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

log_scale = alt.Scale(type="log")
background_chart = (
    alt.Chart(background_df)
    .mark_line(opacity=0.1, strokeWidth=1)
    .encode(
        x=alt.X("Time (min):Q", title="Time (min)"),
        y=alt.Y("counts:Q", title="Log scale counts", scale=log_scale),
        detail="lineage_seed:N",  # IMPORTANT
    )
)

highlight_chart = (
    alt.Chart(highlight_df)
    .mark_line(strokeWidth=3, color="black")
    .encode(x="Time (min):Q", y=alt.Y("counts:Q", scale=log_scale))
)

chart = (background_chart + highlight_chart).properties(
    title="Counts Across Multiple Seeds (log scale)",
    width=600,
    height=400,
)
# %%
outdir = "out/manual_DnaA_analysis"
os.makedirs(outdir, exist_ok=True)

# %%
# --- SAVE ---
output_df.to_csv(os.path.join(outdir, "DnaAboxescounts_multiseed.csv"), index=False)

html_path = os.path.join(outdir, "DnaAboxescounts_multiseed.html")
chart.save(html_path)

# %%
bulk_sql = read_stacked_columns(
    history_sql,
    [
        "bulk",
    ],
    order_results=False,
)
query2 = f"""
            SELECT time, bulk, lineage_seed
            FROM ({bulk_sql})
            ORDER BY lineage_seed,time
            """
output_df2 = conn.sql(query2).df()
time_minutes = output_df2["time"].to_numpy() / 60
bulk_matrix = np.stack(output_df2["bulk"].values).astype(int)
bulk_df = pd.DataFrame(bulk_matrix)
bulk_df.insert(0, "Time (min)", time_minutes)
bulk_df.insert(1, "lineage_seed", output_df2["lineage_seed"].to_numpy())
bulk_df.to_csv(os.path.join(outdir, "bulk_matrix.csv"), index=False)

# DnaA proteins
DnaA_cols = [11525, 10782]
# Create a DataFrame with just Time and DnaAproteins
selected_cols = (
    ["Time (min)"] + [bulk_df.columns[i] for i in DnaA_cols] + ["lineage_seed"]
)
DnaAprotein_df = bulk_df[selected_cols].copy()
DnaAprotein_df.columns = ["Time (min)", "DnaA", "DnaA-ATP", "lineage_seed"]

# Nucleotides
nucleotide_cols = [550, 6212, 6223, 6300, 12450, 798, 5980, 10002, 15593, 3499]

nucleotide_selected_cols = (
    ["Time (min)"] + [bulk_df.columns[i] for i in nucleotide_cols] + ["lineage_seed"]
)
nucleotide_df = bulk_df[nucleotide_selected_cols].copy()
# Rename nucleotide columns for clarity
nucleotide_names = [
    "ADP",
    "DATP",
    "DCTP",
    "DGTP",
    "TTP",
    "ATP",
    "CTP",
    "GTP",
    "UTP",
    "Phosphatidylglycerol",
    "lineage_seed",
]
nucleotide_df.columns = ["Time (min)"] + nucleotide_names

# Merge DnaA and nucleotide columns into one wide table
combined_df = DnaAprotein_df.merge(nucleotide_df, on="Time (min)")

# Merge DnaA proteins + nucleotides into one wide table
bulk_combined_df = pd.merge(
    DnaAprotein_df, nucleotide_df, on=["Time (min)", "lineage_seed"]
)

# %%
output_df = output_df.rename(columns={"DnaA_box_counts": "DnaA-box"})
# Merge the DnaA box counts with bulk_combined_df

# %%
full_df = pd.merge(
    output_df[["Time (min)", "lineage_seed", "DnaA-box"]],
    bulk_combined_df,
    on=["Time (min)", "lineage_seed"],
)

# %%
log_scale = alt.Scale(type="log")

# Melt full_df if not already done
melted_full_df = full_df.melt(
    id_vars=["Time (min)", "lineage_seed"], var_name="Molecule", value_name="Counts"
)
melted_full_df["Counts"] = melted_full_df["Counts"] + 1
# Pick first seed
first_seed = melted_full_df["lineage_seed"].unique()[0]

# Split into background and highlight
background_df = melted_full_df[melted_full_df["lineage_seed"] != first_seed]
highlight_df = melted_full_df[melted_full_df["lineage_seed"] == first_seed]

# %%
# Background lines: all seeds except median
background_chart = (
    alt.Chart(background_df)
    .mark_line(opacity=0.25, strokeWidth=1)
    .encode(
        x=alt.X("Time (min):Q", title="Time (min)"),
        y=alt.Y(
            "Counts:Q", scale=alt.Scale(type="log"), title="Counts (+1, log scale)"
        ),
        color=alt.Color(
            "Molecule:N", scale=alt.Scale(scheme="tableau20"), legend=None
        ),  # hides legend
        detail="lineage_seed:N",
    )
)

# Highlight lines: fully opaque, controls legend
highlight_chart = (
    alt.Chart(highlight_df)
    .mark_line(strokeWidth=3)
    .encode(
        x="Time (min):Q",
        y=alt.Y("Counts:Q", scale=alt.Scale(type="log")),
        color=alt.Color(
            "Molecule:N",
            scale=alt.Scale(scheme="tableau20"),
            legend=alt.Legend(
                title="Molecule",
                labelFontSize=12,
                titleFontSize=14,
                symbolType="stroke",
            ),
        ),
    )
)


# --- Combine charts ---
chart = (
    alt.layer(background_chart, highlight_chart)
    .configure_axis(
        grid=False,
        labelFontSize=12,
        titleFontSize=14,
    )
    .properties(width=800, height=500)
)

# Save figure
html_path = os.path.join(outdir, "all_molecules_multiseed.html")
chart.save(html_path)

# %%
legend_chart = (
    alt.Chart(melted_full_df)
    .mark_line(strokeWidth=4)  # thick so legend is bold
    .encode(
        color=alt.Color(
            "Molecule:N",
            scale=alt.Scale(scheme="tableau20"),
            legend=alt.Legend(
                orient="right",
                labelFontSize=14,
                titleFontSize=16,
                symbolType="stroke",  # line legend
            ),
        )
    )
    .transform_aggregate(
        count="count()", groupby=["Molecule"]
    )  # one entry per molecule
    .encode(y=alt.value(0))  # collapse into a single row
    .properties(width=200, height=200)
    .configure_axis(grid=False, domain=False, ticks=False, labels=False)
)

legend_chart.save(os.path.join(outdir, "legend_only.html"))

# %%
