# %%
import os
import numpy as np
import duckdb

import altair as alt
import matplotlib.pyplot as plt
import pandas as pd

import seaborn as sns


# from ecoli.analysis.multiseed.DnaA_Analysis_new import highlight_df
# from ecoli.analysis.multiseed.DnaA_Analysis_new import nucleotide_cols
from ecoli.library.parquet_emitter import (
    read_stacked_columns,
)

# Loading Simdata to obtain the bulk id labels
from ecoli.library.sim_data import LoadSimData
from ecoli.library.schema import bulk_name_to_idx
from scipy.constants import N_A


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
# %%
sim_data_default = "out/multiseed_noDperiod/parca/kb/simData.cPickle"
sim_data = LoadSimData(sim_data_default).sim_data

bulk_molecule_ids = sim_data.internal_state.bulk_molecules.bulk_data["id"].tolist()
bulk_id_df = pd.DataFrame(bulk_molecule_ids, columns=["bulk_id"])

conn = duckdb.connect()
query_dict = {
    "experiment_id": "multiseed_noDperiod",
    "variant": 0,
    "generation": 1,
}
history_sql = f"""
SELECT * FROM read_parquet('out/{query_dict["experiment_id"]}/history/*/*/*/*/*/*.pq', hive_partitioning=true)
"""
outdir = "out/manual_DnaA_analysis"
os.makedirs(outdir, exist_ok=True)
# %%
cell_mass_sql = read_stacked_columns(
    history_sql,
    ["listeners__mass__cell_mass"],
    order_results=False,
)
dna_mass_sql = read_stacked_columns(
    history_sql,
    ["listeners__mass__dna_mass"],
    order_results=False,
)
mrna_mass_sql = read_stacked_columns(
    history_sql,
    ["listeners__mass__mRna_mass"],
    order_results=False,
)
protein_mass_sql = read_stacked_columns(
    history_sql,
    ["listeners__mass__protein_mass"],
    order_results=False,
)
smallMolecule_mass_sql = read_stacked_columns(
    history_sql,
    ["listeners__mass__smallMolecule_mass"],
    order_results=False,
)

query_2 = f"""
SELECT time, lineage_seed, listeners__mass__cell_mass FROM ({cell_mass_sql})
ORDER BY lineage_seed, time
"""
output_df2 = conn.sql(query_2).df()
# Convert time from seconds to minutes
output_df2["Time (min)"] = output_df2["time"] / 60

query_3 = f"""
SELECT time,lineage_seed, listeners__mass__dna_mass FROM ({dna_mass_sql})
ORDER BY lineage_seed, time
"""

output_df3 = conn.sql(query_3).df()

query_4 = f"""
SELECT time, lineage_seed, listeners__mass__mRna_mass FROM ({mrna_mass_sql})
ORDER BY lineage_seed, time 
"""
output_df4 = conn.sql(query_4).df()

query_5 = f"""
SELECT time, lineage_seed, listeners__mass__protein_mass FROM ({protein_mass_sql})
ORDER BY lineage_seed, time
"""
output_df5 = conn.sql(query_5).df()

query_6 = f"""
SELECT time, lineage_seed, listeners__mass__smallMolecule_mass FROM ({smallMolecule_mass_sql})
ORDER BY lineage_seed, time
"""
output_df6 = conn.sql(query_6).df()
# %%
bulk_sql = read_stacked_columns(
    history_sql,
    ["bulk"],
    order_results=False,
)
query_7 = f"""
            SELECT time, bulk, lineage_seed
            FROM ({bulk_sql})
            ORDER BY lineage_seed,time
            """
output_df7 = conn.sql(query_7).df()
time_minutes = output_df7["time"].to_numpy() / 60
bulk_matrix = np.stack(output_df7["bulk"].values).astype(int)
bulk_df = pd.DataFrame(bulk_matrix)
bulk_df.insert(0, "Time (min)", time_minutes)
bulk_df.insert(1, "lineage_seed", output_df7["lineage_seed"].to_numpy())

# %%
dntp_ids = bulk_name_to_idx(
    ["DATP[c]", "DCTP[c]", "DGTP[c]", "TTP[c]"], bulk_molecule_ids
)

ntp_ids = bulk_name_to_idx(["ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]"], bulk_molecule_ids)
adp_id = bulk_name_to_idx("ADP[c]", bulk_molecule_ids)
# %%
dnucleotide_selected_cols = (
    ["Time (min)"] + ["lineage_seed"] + [bulk_df.columns[i + 2] for i in dntp_ids]
)
dnucleotide_df = bulk_df[dnucleotide_selected_cols].copy()
# Rename nucleotide columns for clarity
dnucleotide_names = ["DATP", "DCTP", "DGTP", "TTP"]
dnucleotide_df.columns = ["Time (min)"] + ["lineage_seed"] + dnucleotide_names
dnucleotide_df["Free_dNTP"] = (
    dnucleotide_df["DATP"]
    + dnucleotide_df["DCTP"]
    + dnucleotide_df["DGTP"]
    + dnucleotide_df["TTP"]
)

ntp_selected_cols = (
    ["Time (min)"] + ["lineage_seed"] + [bulk_df.columns[i + 2] for i in ntp_ids]
)
ntp_df = bulk_df[ntp_selected_cols].copy()
# Rename nucleotide columns for clarity
ntp_names = ["ATP", "CTP", "GTP", "UTP"]
ntp_df.columns = ["Time (min)"] + ["lineage_seed"] + ntp_names
ntp_df["Free_NTP"] = ntp_df["ATP"] + ntp_df["CTP"] + ntp_df["GTP"] + ntp_df["UTP"]

adp_selected_cols = ["Time (min)"] + ["lineage_seed"] + [bulk_df.columns[adp_id + 2]]
adp_df = bulk_df[adp_selected_cols].copy()
adp_df.columns = ["Time (min)"] + ["lineage_seed"] + ["ADP"]
# %%
for df in [output_df2, output_df3, output_df4, output_df5, output_df6]:
    df["Time (min)"] = df["time"] / 60

# Step 2: Rename the mass columns for clarity
output_df2 = output_df2.rename(columns={"listeners__mass__cell_mass": "Cell_mass"})
output_df3 = output_df3.rename(columns={"listeners__mass__dna_mass": "DNA_mass"})
output_df4 = output_df4.rename(columns={"listeners__mass__mRna_mass": "mRNA_mass"})
output_df5 = output_df5.rename(
    columns={"listeners__mass__protein_mass": "Protein_mass"}
)
output_df6 = output_df6.rename(
    columns={"listeners__mass__smallMolecule_mass": "SmallMolecule_mass"}
)
# %%
# Step 3: Merge all dataframes on 'Time (min)'
merged_df = (
    output_df2[["Time (min)", "Cell_mass", "lineage_seed"]]
    .merge(
        output_df3[["Time (min)", "DNA_mass", "lineage_seed"]],
        on=["Time (min)", "lineage_seed"],
    )
    .merge(
        output_df4[["Time (min)", "mRNA_mass", "lineage_seed"]],
        on=["Time (min)", "lineage_seed"],
    )
    .merge(
        output_df5[["Time (min)", "Protein_mass", "lineage_seed"]],
        on=["Time (min)", "lineage_seed"],
    )
    .merge(
        output_df6[["Time (min)", "SmallMolecule_mass", "lineage_seed"]],
        on=["Time (min)", "lineage_seed"],
    )
    .merge(
        dnucleotide_df[["Time (min)", "Free_dNTP", "lineage_seed"]],
        on=["Time (min)", "lineage_seed"],
    )
    .merge(
        ntp_df[["Time (min)", "Free_NTP", "lineage_seed"]],
        on=["Time (min)", "lineage_seed"],
    )
    .merge(
        adp_df[["Time (min)", "ADP", "lineage_seed"]], on=["Time (min)", "lineage_seed"]
    )
)

# Step 4 (optional): Convert to long format for plotting multiple molecules
long_df_mass = merged_df.melt(
    id_vars=["Time (min)", "lineage_seed"], var_name="Molecule", value_name="Mass"
)
# %%
# Creating separate plot for nucleotide counts
merged_df = (
    output_df4[["Time (min)", "mRNA_mass", "lineage_seed"]]
    .merge(
        dnucleotide_df[["Time (min)", "Free_dNTP", "lineage_seed"]],
        on=["Time (min)", "lineage_seed"],
    )
    .merge(
        ntp_df[["Time (min)", "Free_NTP", "lineage_seed"]],
        on=["Time (min)", "lineage_seed"],
    )
    .merge(
        adp_df[["Time (min)", "ADP", "lineage_seed"]], on=["Time (min)", "lineage_seed"]
    )
)

long_df_mass = merged_df.melt(
    id_vars=["Time (min)", "lineage_seed"], var_name="Molecule", value_name="Mass"
)
# %%
df = long_df_mass.copy()
# Get t0 per (seed, molecule)
t0_df = (
    df.sort_values("Time (min)")
    .groupby(["lineage_seed", "Molecule"], as_index=False)
    .first()
)
t0_map = t0_df.set_index(["lineage_seed", "Molecule"])["Mass"]
# Map baseline to all rows
df["t0"] = df.set_index(["lineage_seed", "Molecule"]).index.map(t0_map)

# Fold change
pseudocount = 1e-6
df["FoldChange"] = (df["Mass"] + pseudocount) / (df["t0"] + pseudocount)
# %%
# Compute median per time and molecule
median_df = (
    df.groupby(["Time (min)", "Molecule"], as_index=False)["FoldChange"]
    .median()
    .rename(columns={"FoldChange": "median_FC"})
)

# Merge median back onto df
df_with_median = df.merge(median_df, on=["Time (min)", "Molecule"])

# Compute absolute difference from median
df_with_median["abs_diff"] = np.abs(
    df_with_median["FoldChange"] - df_with_median["median_FC"]
)

# Compute average absolute difference per seed and molecule
seed_distance = (
    df_with_median.groupby(["lineage_seed", "Molecule"], as_index=False)["abs_diff"]
    .mean()  # use mean to avoid bias toward longer trajectories
    .rename(columns={"abs_diff": "avg_diff"})
)

# Choose reference molecule (e.g., DNA_mass)
reference_molecule = "DNA_mass"
highlight_seed = (
    seed_distance.loc[seed_distance["Molecule"] == reference_molecule]
    .sort_values("avg_diff")
    .iloc[0]["lineage_seed"]
)
# %%
# Select seed for nucleotides
# highlight_seed = 4
print("Selected seed for all molecules:", highlight_seed)

# Filter df for the selected seed across all molecules
highlight_df = df[df["lineage_seed"] == highlight_seed]
# %%
all_molecules = highlight_df["Molecule"].unique()

color_encoding = alt.Color(
    "Molecule:N",
    scale=alt.Scale(
        domain=all_molecules,  # ensures every molecule is assigned a color
        range=COLORS,
    ),
    legend=alt.Legend(title="Molecule", labelColor="black", titleColor="black"),
)
# Set x-axis scale to terminate at max time
t_max = df["Time (min)"].max()
x_scale = alt.X(
    "Time (min):Q",
    scale=alt.Scale(domain=[0, t_max]),
    axis=alt.Axis(title="Time (min)"),
)
# %%
# deviation from representative
rep_df = highlight_df.rename(columns={"FoldChange": "Rep_FC"})

merged = df.merge(
    rep_df[["Time (min)", "Molecule", "Rep_FC"]], on=["Time (min)", "Molecule"]
)

merged["Deviation"] = merged["FoldChange"] - merged["Rep_FC"]

band_df = merged.groupby(["Time (min)", "Molecule"], as_index=False).agg(
    lower=("Deviation", lambda x: x.quantile(0.25)),
    upper=("Deviation", lambda x: x.quantile(0.75)),
    rep=("Rep_FC", "first"),
)

band_df["lower"] = band_df["rep"] + band_df["lower"]
band_df["upper"] = band_df["rep"] + band_df["upper"]

band_chart = (
    alt.Chart(band_df)
    .mark_area(opacity=0.15)
    .encode(x=x_scale, y="lower:Q", y2="upper:Q", color=color_encoding)
)
# %%
# Get max time
# t_max = df["Time (min)"].max()

# Background (blurred seeds)
background_chart = (
    alt.Chart(df)
    .mark_line(opacity=0.2, strokeWidth=1.5)
    .encode(
        x=x_scale,
        y=alt.Y(
            "FoldChange:Q",
            title="Fold change (t=0 = 1)",
            scale=alt.Scale(domain=[0.8, 3]),
            # axis=alt.Axis(values=[0.1,0.2,0.4,0.6, 1, 2,3]
            # ),
        ),
        color=alt.Color(color_encoding.shorthand, legend=None),
        detail=["lineage_seed:N", "Molecule:N"],
    )
)

# Highlighted seed
median_line = (
    alt.Chart(highlight_df)
    .mark_line(strokeWidth=4)
    .encode(
        x=x_scale,
        y=alt.Y(
            "FoldChange:Q",
            scale=alt.Scale(domain=[0.8, 3]),
            # axis=alt.Axis(values=[0.1,0.2,0.4,0.6, 1,2, 3]
            # ),
        ),
        color=color_encoding,
        # strokeDash=alt.StrokeDash("Molecule:N")
    )
)
# Times in minutes
# t_rif = 10
t_rep_done = 32.12
# new initiation
t_init = 38.25

# Darker gray: Rifampicin
# shade_rif = (
#     alt.Chart(pd.DataFrame({
#         "x_start": [t_rif],
#         "x_end": [110]
#     }))
#     .mark_rect(opacity=0.10, color="gray")  # darker
#     .encode(
#         x="x_start:Q",
#         x2="x_end:Q"
#     )
# )

# Lighter gray: post-replication
shade_rep = (
    alt.Chart(pd.DataFrame({"x_start": [t_rep_done], "x_end": [t_init]}))
    .mark_rect(opacity=0.04, color="gray")  # lighter
    .encode(x="x_start:Q", x2="x_end:Q")
)

shade_init = (
    alt.Chart(pd.DataFrame({"x_start": [t_init], "x_end": [t_max]}))
    .mark_rect(opacity=0.1, color="gray")  # lighter
    .encode(x="x_start:Q", x2="x_end:Q")
)

# Combine
chart = (
    alt.layer(
        # shade_rif,
        shade_rep,
        shade_init,
        background_chart,
        # band_chart,
        median_line,
    )
    .configure_axis(grid=False, labelFontSize=14, titleFontSize=16)
    .configure_legend(labelFontSize=14, titleFontSize=16)
    .configure_title(fontSize=18, fontWeight="bold")
    .properties(
        width=800,
        height=500,
        padding={"bottom": 70},
        title="Fold change of mass after rifampicin addition (10 min)",
    )
    .resolve_scale(x="shared")
)

# Save
chart.save(os.path.join(outdir, "foldchange_combined_noDperiod.svg"))

# %%

# --- 2. Create a small "dummy" chart for the legend ---
# Use just one row per molecule
legend_df = pd.DataFrame(
    {
        "Molecule": df["Molecule"].unique(),
        "y": [1] * len(df["Molecule"].unique()),  # dummy y value
    }
)

legend_chart = (
    alt.Chart(legend_df)
    .mark_line(strokeWidth=4)
    .encode(
        x=alt.X("Molecule:N", axis=None),
        y=alt.Y("y:Q", axis=None),
        color=color_encoding,
    )
    .properties(width=800, height=500, title="Molecule legend")
    .configure_title(fontSize=16, fontWeight="bold")
)

legend_chart.save(os.path.join(outdir, "10minrifamfoldchange_legend.html"))

# %%
# Generating heatmap for representative seed


# Use the representative seed
rep_seed = 4
rep_df = df[df["lineage_seed"] == rep_seed].copy()

# Pivot so rows are Molecules, columns are Time, values are FoldChange
heatmap_df = rep_df.pivot(index="Molecule", columns="Time (min)", values="FoldChange")

# Sort molecules if you want
heatmap_df = heatmap_df.loc[
    [
        "Cell_mass",
        "DNA_mass",
        "mRNA_mass",
        "Protein_mass",
        "SmallMolecule_mass",
        "Free_dNTP",
        "Free_NTP",
        "ADP",
    ]
]

# Set up figure
plt.figure(figsize=(30, 8))
# Draw heatmap
ax = sns.heatmap(
    heatmap_df,
    cmap="bwr",
    center=1,  # red-white-blue diverging
    cbar_kws={"label": "Fold change"},
    # linewidths=0.1,
    # linecolor='gray'
)

# Original time values
time_cols = heatmap_df.columns.values

# Keep only integer times
integer_times = [t for t in time_cols if t % 1 == 0]

# Select every 5 minutes
tick_times = [t for t in integer_times if t % 5 == 0]

# Get positions in the heatmap (indices of columns)
tick_positions = [i for i, t in enumerate(time_cols) if t in tick_times]

# Tick labels
tick_labels = [str(int(time_cols[i])) for i in tick_positions]

# Apply ticks
ax.set_xticks(tick_positions)
ax.set_xticklabels(tick_labels, rotation=0, fontsize=20)
ax.set_yticklabels(ax.get_yticklabels(), fontsize=20)
plt.xlabel("Time (min)", fontsize=16)

# Colorbar customization
cbar = ax.collections[0].colorbar
cbar.set_label("Fold change", fontsize=20)  # colorbar label
cbar.ax.tick_params(labelsize=20)

plt.savefig("heatmap_foldchange.png", dpi=300)
# %%
query_dict_2 = {
    "experiment_id": "DnaAboxes",
    "variant": 0,
    "lineage_seed": 0,
    "generation": 1,
}
history_sql_2 = f"""
SELECT * FROM read_parquet('out/{query_dict_2["experiment_id"]}/history/*/*/*/*/*/*.pq', hive_partitioning=true)
"""
query_7 = f"""
SELECT time, listeners__mass__cell_mass FROM ({history_sql_2})
ORDER BY time ASC
"""
output_df7 = conn.sql(query_7).df()

query_8 = f"""
SELECT time, listeners__mass__dna_mass FROM ({history_sql_2})
ORDER BY time ASC
"""
output_df8 = conn.sql(query_8).df()

query_9 = f"""
SELECT time, listeners__mass__rna_mass FROM ({history_sql_2})
ORDER BY time ASC"""

output_df9 = conn.sql(query_9).df()

query_10 = f"""
SELECT time, listeners__mass__protein_mass FROM ({history_sql_2})
ORDER BY time ASC"""

output_df10 = conn.sql(query_10).df()

query_11 = f"""
SELECT time, listeners__mass__smallMolecule_mass FROM ({history_sql_2})
"""
output_df11 = conn.sql(query_11).df()

# %%
for df in [output_df7, output_df8, output_df9, output_df10, output_df11]:
    df["Time (min)"] = df["time"] / 60

# Step 2: Rename the mass columns for clarity
output_df7 = output_df7.rename(
    columns={"listeners__mass__cell_mass": "Cell_mass_default"}
)
output_df8 = output_df8.rename(
    columns={"listeners__mass__dna_mass": "DNA_mass_default"}
)
output_df9 = output_df9.rename(
    columns={"listeners__mass__rna_mass": "RNA_mass_default"}
)
output_df10 = output_df10.rename(
    columns={"listeners__mass__protein_mass": "Protein_mass_default"}
)
output_df11 = output_df11.rename(
    columns={"listeners__mass__smallMolecule_mass": "SmallMolecule_mass_default"}
)

# Step 3: Merge all dataframes on 'Time (min)'
merged_df = (
    output_df7[["Time (min)", "Cell_mass_default"]]
    .merge(output_df8[["Time (min)", "DNA_mass_default"]], on="Time (min)")
    .merge(output_df9[["Time (min)", "RNA_mass_default"]], on="Time (min)")
    .merge(output_df10[["Time (min)", "Protein_mass_default"]], on="Time (min)")
    .merge(output_df11[["Time (min)", "SmallMolecule_mass_default"]], on="Time (min)")
)

# Step 4 (optional): Convert to long format for plotting multiple molecules
long_df_default_mass = merged_df.melt(
    id_vars=["Time (min)"], var_name="Molecule", value_name="Mass"
)

# %%
long_df_default_mass["Molecule"] = long_df_default_mass["Molecule"].str.replace(
    "_default", ""
)
first_seed = long_df_mass["lineage_seed"].unique()[0]

background_df = long_df_mass[long_df_mass["lineage_seed"] != first_seed]
highlight_df = long_df_mass[long_df_mass["lineage_seed"] == first_seed]

color_encoding = alt.Color(
    "Molecule:N",
    scale=alt.Scale(scheme="tableau20"),
    legend=alt.Legend(title="Molecule", labelColor="black", titleColor="black"),
)
# %%
background_chart = (
    alt.Chart(background_df)
    .mark_line(opacity=0.25, strokeWidth=1.5)
    .encode(
        x=alt.X("Time (min):Q", axis=alt.Axis(titlePadding=30)),
        y=alt.Y("Mass:Q", scale=alt.Scale(type="log")),
        color=alt.Color("Molecule:N", scale=alt.Scale(scheme="tableau20"), legend=None),
        detail="lineage_seed:N",
    )
)

highlight_chart = (
    alt.Chart(highlight_df)
    .mark_line(strokeWidth=4)
    .encode(
        x="Time (min):Q",
        y=alt.Y("Mass:Q", scale=alt.Scale(type="log")),
        color=color_encoding,
    )
)

default_chart = (
    alt.Chart(long_df_default_mass)
    .mark_line(strokeDash=[5, 5], strokeWidth=2)
    .encode(
        x="Time (min):Q",
        y=alt.Y("Mass:Q", scale=alt.Scale(type="log")),
        color=alt.Color("Molecule:N", scale=alt.Scale(scheme="tableau20"), legend=None),
    )
)
# %%
t_inhibit = 10 / 60
t_replication_done = 1327 / 60

events_df = pd.DataFrame(
    {
        "Time (min)": [t_inhibit, t_replication_done],
        "Event": ["Transcription Inhibition at 10s", "Ongoing Replication Completed"],
    }
)

vlines = (
    alt.Chart(events_df)
    .mark_rule(color="black", opacity=0.3, size=2)
    .encode(x="Time (min):Q")
)

labels = (
    alt.Chart(events_df)
    .mark_text(angle=0, dy=15, align="center", fontSize=14)
    .encode(
        x="Time (min):Q", y=alt.value(520), text="Event:N", color=alt.value("black")
    )
)

chart = (
    alt.layer(background_chart, default_chart, highlight_chart, vlines, labels)
    .configure_axis(grid=False, labelFontSize=14, titleFontSize=16)
    .configure_legend(labelFontSize=14, titleFontSize=16)
    .configure_title(fontSize=18, fontWeight="bold")
    .properties(
        width=800,
        height=500,
        padding={"bottom": 70},
        title="Mass trajectories: transcription inhibition vs baseline",
    )
)

# Save
html_path = os.path.join(outdir, "runout_multiseed_with_default.html")
chart.save(html_path)

# %%
# --- Separate Legend Chart ---
legend_chart = (
    alt.Chart(long_df_mass)
    .mark_line(strokeWidth=4)  # thick lines for visibility
    .encode(
        color=alt.Color(
            "Molecule:N",
            scale=alt.Scale(scheme="tableau20"),  # SAME as main chart
            legend=alt.Legend(
                title="Molecule",
                orient="right",
                labelFontSize=14,
                titleFontSize=16,
                symbolType="stroke",  # line legend (not circles)
                labelColor="black",
                titleColor="black",
            ),
        )
    )
    .transform_aggregate(
        count="count()", groupby=["Molecule"]
    )  # one entry per molecule
    .encode(y=alt.value(0))  # collapse into one row
    .properties(width=200, height=200)
    .configure_axis(grid=False, domain=False, ticks=False, labels=False)
)

legend_chart.save(os.path.join(outdir, "legend_only_runout_multiseed.html"))

# %%
long_df_default_1min = long_df_default_mass[
    long_df_default_mass["Time (min)"] <= 10
].copy()
long_df_default_1min["Molecule"] = long_df_default_1min["Molecule"].str.replace(
    "_default", ""
)
# %%
# 1 min figure
long_df_1min = long_df_mass[long_df_mass["Time (min)"] <= 10].copy()
events_df_1min = events_df[events_df["Time (min)"] <= 1].copy()
long_df_1min["Time (s)"] = long_df_1min["Time (min)"] * 60
long_df_default_1min["Time (s)"] = long_df_default_1min["Time (min)"] * 60
events_df_1min["Time (s)"] = events_df_1min["Time (min)"] * 60

# %%
# Highlight first lineage seed
first_seed = long_df_1min["lineage_seed"].unique()[0]
long_df_1min["is_highlight"] = long_df_1min["lineage_seed"] == first_seed
background_df_1min = long_df_1min[~long_df_1min["is_highlight"]]
highlight_df_1min = long_df_1min[long_df_1min["is_highlight"]]

color_encoding = alt.Color(
    "Molecule:N",
    scale=alt.Scale(scheme="tableau20"),
    legend=alt.Legend(title="Molecule", labelColor="black", titleColor="black"),
)

background_chart = (
    alt.Chart(background_df_1min)
    .mark_line(opacity=0.25, strokeWidth=1.5)
    .encode(
        x=alt.X("Time (s):Q", axis=alt.Axis(title="Time (s)", titlePadding=30)),
        y=alt.Y("Mass:Q", scale=alt.Scale(type="log")),
        color=alt.Color("Molecule:N", scale=alt.Scale(scheme="tableau20"), legend=None),
        detail="lineage_seed:N",
    )
)

highlight_chart = (
    alt.Chart(highlight_df_1min)
    .mark_line(strokeWidth=4)
    .encode(
        x=alt.X("Time (s):Q", axis=alt.Axis(title="Time (s)", titlePadding=30)),
        y=alt.Y("Mass:Q", scale=alt.Scale(type="log")),
        color=color_encoding,
    )
)

default_chart = (
    alt.Chart(long_df_default_1min)
    .mark_line(strokeDash=[5, 5], strokeWidth=2)
    .encode(
        x=alt.X("Time (s):Q", axis=alt.Axis(title="Time (s)", titlePadding=30)),
        y=alt.Y("Mass:Q", scale=alt.Scale(type="log")),
        color=alt.Color("Molecule:N", scale=alt.Scale(scheme="tableau20"), legend=None),
        detail="Molecule:N",
    )
)

# Events (if within first 60s)
vlines_1min = (
    alt.Chart(events_df_1min)
    .mark_rule(color="black", opacity=0.3, size=2)
    .encode(x="Time (s):Q")
)
labels_1min = (
    alt.Chart(events_df_1min)
    .mark_text(angle=0, dy=15, align="center", fontSize=14)
    .encode(x="Time (s):Q", y=alt.value(520), text="Event:N", color=alt.value("black"))
)

chart = (
    alt.layer(
        background_chart,
        default_chart,
        highlight_chart,
        vlines_1min,
        labels_1min,
    )
    .configure_axis(grid=False, labelFontSize=14, titleFontSize=16)
    .configure_legend(labelFontSize=14, titleFontSize=16)
    .configure_title(fontSize=18, fontWeight="bold")
    .properties(
        width=800,
        height=500,
        padding={"bottom": 70},
        title="Mass trajectories (first 60s): transcription inhibition vs baseline",
    )
)

chart.save(os.path.join(outdir, "runout_zoom_60s.html"))

# %%
DnaA_box_counts_sql = read_stacked_columns(
    history_sql,
    [
        "listeners__unique_molecule_counts__DnaA_box",
    ],
    order_results=False,
)

query_12 = f"""
SELECT lineage_seed, time, listeners__unique_molecule_counts__DnaA_box FROM ({DnaA_box_counts_sql})
ORDER BY lineage_seed, time
"""
output_df12 = conn.sql(query_12).df()
# %%
# Convert time from seconds to minutes
output_df12["Time (min)"] = output_df12["time"] / 60
# %%

bulk_sql = read_stacked_columns(
    history_sql,
    [
        "bulk",
    ],
    order_results=False,
)
query_13 = f"""
            SELECT time, bulk, lineage_seed
            FROM ({bulk_sql})
            ORDER BY lineage_seed,time
            """
output_df13 = conn.sql(query_13).df()
time_minutes = output_df13["time"].to_numpy() / 60
bulk_matrix = np.stack(output_df13["bulk"].values).astype(int)
bulk_df = pd.DataFrame(bulk_matrix)
bulk_df.insert(0, "Time (min)", time_minutes)
bulk_df.insert(1, "lineage_seed", output_df13["lineage_seed"].to_numpy())

# %%
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
output_df12 = output_df12.rename(
    columns={"listeners__unique_molecule_counts__DnaA_box": "DnaA-box"}
)
# Merge the DnaA box counts with bulk_combined_df
# %%
full_df = pd.merge(
    output_df12[["Time (min)", "lineage_seed", "DnaA-box"]],
    bulk_combined_df,
    on=["Time (min)", "lineage_seed"],
)

# %%
query_14 = f"""
SELECT time, listeners__unique_molecule_counts__DnaA_box FROM ({history_sql_2})
ORDER BY time ASC
"""
output_df14 = conn.sql(query_14).df()
output_df14["Time (min)"] = output_df14["time"] / 60

query_15 = f"""
SELECT time, bulk, FROM ({history_sql_2})
ORDER BY time ASC
"""
output_df15 = conn.sql(query_15).df()
time_minutes = output_df15["time"].to_numpy() / 60
bulk_matrix_default = np.stack(output_df15["bulk"].values).astype(int)
bulk_df_default = pd.DataFrame(bulk_matrix_default)
bulk_df_default.insert(0, "Time (min)", time_minutes)

# DnaA proteins
DnaA_cols = [11524, 10781]
# Create a DataFrame with just Time and DnaAproteins
selected_cols = ["Time (min)"] + [bulk_df_default.columns[i] for i in DnaA_cols]
DnaAprotein_df_default = bulk_df_default[selected_cols].copy()
DnaAprotein_df_default.columns = ["Time (min)", "DnaA", "DnaA-ATP"]

# Nucleotides
nucleotide_cols = [549, 6211, 6222, 6299, 12449, 797, 5979, 10001, 15592, 3498]

nucleotide_selected_cols = ["Time (min)"] + [
    bulk_df_default.columns[i] for i in nucleotide_cols
]
nucleotide_df_default = bulk_df_default[nucleotide_selected_cols].copy()
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
]
nucleotide_df_default.columns = ["Time (min)"] + nucleotide_names

# Merge DnaA and nucleotide columns into one wide table
combined_df_default = DnaAprotein_df_default.merge(
    nucleotide_df_default, on="Time (min)"
)

# Merge DnaA proteins + nucleotides into one wide table
bulk_combined_df_default = pd.merge(
    DnaAprotein_df_default, nucleotide_df_default, on="Time (min)"
)

# %%
long_df_counts = full_df.melt(
    id_vars=["Time (min)", "lineage_seed"], var_name="Molecule", value_name="Counts"
)
long_df_default_counts = bulk_combined_df_default.melt(
    id_vars=["Time (min)"], var_name="Molecule", value_name="Counts"
)
long_df_counts["Counts"] = long_df_counts["Counts"] + 1
long_df_default_counts["Counts"] = long_df_default_counts["Counts"] + 1

first_seed = long_df_counts["lineage_seed"].unique()[0]
long_df_counts["is_highlight"] = long_df_counts["lineage_seed"] == first_seed
background_df = long_df_counts[~long_df_counts["is_highlight"]]
highlight_df = long_df_counts[long_df_counts["is_highlight"]]

color_encoding = alt.Color(
    "Molecule:N",
    scale=alt.Scale(scheme="tableau20"),
    legend=alt.Legend(title="Molecule", labelColor="black", titleColor="black"),
)

background_chart = (
    alt.Chart(background_df)
    .mark_line(opacity=0.25, strokeWidth=1.5)
    .encode(
        x=alt.X("Time (min):Q", axis=alt.Axis(title="Time (min)", titlePadding=30)),
        y=alt.Y("Counts:Q", scale=alt.Scale(type="log"), title="Counts (log scale)"),
        color=alt.Color("Molecule:N", scale=alt.Scale(scheme="tableau20"), legend=None),
        detail="lineage_seed:N",
    )
)
highlight_chart = (
    alt.Chart(highlight_df)
    .mark_line(strokeWidth=4)
    .encode(
        x="Time (min):Q",
        y=alt.Y("Counts:Q", scale=alt.Scale(type="log")),
        color=color_encoding,
    )
)

default_chart = (
    alt.Chart(long_df_default_counts)
    .mark_line(strokeDash=[5, 5], strokeWidth=2)
    .encode(
        x="Time (min):Q",
        y=alt.Y("Counts:Q", scale=alt.Scale(type="log")),
        color=alt.Color("Molecule:N", scale=alt.Scale(scheme="tableau20"), legend=None),
        detail="Molecule:N",
    )
)

t_inhibit = 10 / 60  # 10 seconds
t_replication_done = 1327 / 60  # 1327 seconds

events_df = pd.DataFrame(
    {
        "Time (min)": [t_inhibit, t_replication_done],
        "Event": ["Transcription Inhibition at 10s", "Ongoing Replication Completed"],
    }
)

vlines = (
    alt.Chart(events_df)
    .mark_rule(color="gray", opacity=0.3, size=4)
    .encode(x="Time (min):Q")
)

labels = (
    alt.Chart(events_df)
    .mark_text(angle=0, dy=15, align="center", fontSize=14)
    .encode(
        x="Time (min):Q", y=alt.value(520), text="Event:N", color=alt.value("black")
    )
)

# --- 6. Combine all layers ---
chart = (
    alt.layer(background_chart, default_chart, highlight_chart, vlines, labels)
    .configure_axis(grid=False, labelFontSize=14, titleFontSize=16)
    .configure_legend(labelFontSize=14, titleFontSize=16)
    .configure_title(fontSize=18, fontWeight="bold")
    .properties(
        width=800,
        height=500,
        padding={"bottom": 70},
        title="Mass trajectories: transcription inhibition vs baseline",
    )
)

# --- 7. Save ---
html_path = os.path.join(outdir, "runout_counts.html")
chart.save(html_path)

# %%
legend_chart = (
    alt.Chart(long_df_counts)
    .mark_line()
    .encode(y=alt.Y("Molecule:N", axis=None), color=color_encoding)
    .properties(width=150, height=300, title="Legend")
)
legend_chart.save(os.path.join(outdir, "runout_counts_legend.html"))

# %%
long_df_default_mass["Molecule"] = long_df_default_mass["Molecule"].str.replace(
    "_default", ""
)

# %%
# combine inhibited data
long_df_combined = pd.concat([long_df_mass, long_df_counts], ignore_index=True)

long_df_default_combined = pd.concat(
    [long_df_default_mass, long_df_default_counts], ignore_index=True
)
# %%
# Combine Mass and Counts into a single column
long_df_combined["Value"] = long_df_combined["Mass"].fillna(long_df_combined["Counts"])

# Drop the old Mass/Counts columns to avoid confusion
long_df_combined = long_df_combined.drop(columns=["Mass", "Counts"])

long_df_default_combined["Value"] = long_df_default_combined["Mass"].fillna(
    long_df_default_combined["Counts"]
)
long_df_default_combined = long_df_default_combined.drop(columns=["Mass", "Counts"])
# %%
# Computing fold change
eps = 1  # to make t=0 = 1

# Multiseed
t0_df = long_df_combined[long_df_combined["Time (min)"] == 0][
    ["lineage_seed", "Molecule", "Value"]
].rename(columns={"Value": "Value_t0"})

long_df_combined = long_df_combined.merge(t0_df, on=["lineage_seed", "Molecule"])

long_df_combined["FoldChange"] = (long_df_combined["Value"] + eps) / (
    long_df_combined["Value_t0"] + eps
)

# Default (single seed)
t0_default = long_df_default_combined[long_df_default_combined["Time (min)"] == 0][
    ["Molecule", "Value"]
].rename(columns={"Value": "Value_t0"})

long_df_default_combined = long_df_default_combined.merge(t0_default, on="Molecule")

long_df_default_combined["FoldChange"] = (long_df_default_combined["Value"] + eps) / (
    long_df_default_combined["Value_t0"] + eps
)

# %%
# --- 1. Highlight first seed ---
first_seed = long_df_combined["lineage_seed"].unique()[0]
long_df_combined["is_highlight"] = long_df_combined["lineage_seed"] == first_seed

background_df = long_df_combined[~long_df_combined["is_highlight"]]
highlight_df = long_df_combined[long_df_combined["is_highlight"]]

# --- 2. Color encoding ---
molecules = long_df_combined["Molecule"].unique().tolist()
color_encoding = alt.Color(
    "Molecule:N",
    scale=alt.Scale(domain=molecules, scheme="tableau20"),
    legend=alt.Legend(title="Molecule", labelColor="black", titleColor="black"),
)

# --- 3. Background (blurred) lines ---
background_chart = (
    alt.Chart(background_df)
    .mark_line(opacity=0.25, strokeWidth=1.5)
    .encode(
        x=alt.X("Time (min):Q", axis=alt.Axis(title="Time (min)", titlePadding=30)),
        y=alt.Y(
            "FoldChange:Q", scale=alt.Scale(type="linear"), title="Fold change (t=0=1)"
        ),
        color=alt.Color("Molecule:N", scale=alt.Scale(scheme="tableau20"), legend=None),
        detail="lineage_seed:N",
    )
)

# --- 4. Highlighted line ---
highlight_chart = (
    alt.Chart(highlight_df)
    .mark_line(strokeWidth=4)
    .encode(x="Time (min):Q", y="FoldChange:Q", color=color_encoding)
)

# --- 5. Default (dotted) line ---
default_chart = (
    alt.Chart(long_df_default_combined)
    .mark_line(strokeDash=[5, 5], strokeWidth=2)
    .encode(
        x="Time (min):Q",
        y="FoldChange:Q",
        color=alt.Color(
            "Molecule:N",
            scale=alt.Scale(domain=molecules, scheme="tableau20"),
            legend=None,
        ),
        detail="Molecule:N",
    )
)

# --- 6. Events ---
t_inhibit = 10 / 60  # 10 s
t_replication_done = 1327 / 60

events_df = pd.DataFrame(
    {
        "Time (min)": [t_inhibit, t_replication_done],
        "Event": ["Transcription Inhibition at 10s", "Ongoing Replication Completed"],
    }
)

vlines = (
    alt.Chart(events_df)
    .mark_rule(color="gray", opacity=0.3, size=4)
    .encode(x="Time (min):Q")
)

labels = (
    alt.Chart(events_df)
    .mark_text(angle=0, dy=15, align="center", fontSize=14)
    .encode(
        x="Time (min):Q",
        y=alt.value(520),  # slightly above 1 on linear scale
        text="Event:N",
        color=alt.value("black"),
    )
)

# %%


# --- 7. Combine all layers ---
chart = (
    alt.layer(background_chart, default_chart, highlight_chart, vlines, labels)
    .configure_axis(grid=False, labelFontSize=14, titleFontSize=16)
    .configure_legend(labelFontSize=14, titleFontSize=16)
    .configure_title(fontSize=18, fontWeight="bold")
    .properties(
        width=800,
        height=500,
        padding={"bottom": 70},
        title="Fold change of mass and counts: transcription inhibition vs baseline",
    )
)

# --- 8. Save ---
html_path = os.path.join(outdir, "foldchange_combined.html")
chart.save(html_path)

# %%
# --- Separate legend chart: only one line per molecule ---
# Pick first lineage seed for each molecule to display in legend
legend_df = long_df_combined[long_df_combined["lineage_seed"] == first_seed].copy()

legend_chart = (
    alt.Chart(legend_df)
    .mark_line(strokeWidth=3)
    .encode(
        x=alt.X("Time (min):Q", axis=None),
        y=alt.Y("FoldChange:Q", axis=None),
        color=alt.Color(
            "Molecule:N",
            scale=alt.Scale(domain=molecules, scheme="tableau20"),
            legend=alt.Legend(title="Molecule"),
        ),
    )
    .properties(width=200, height=200, title="Legend")
)

# Save the legend separately
legend_chart.save(os.path.join(outdir, "foldchange_legend.html"))

# %%
# --- Compute log fold change ---
# Add a small pseudocount to avoid log(0)
pseudocount = 1e-6
long_df_combined["FoldChange_log"] = np.log10(
    long_df_combined["FoldChange"] + pseudocount
)
long_df_default_combined["FoldChange_log"] = np.log10(
    long_df_default_combined["FoldChange"] + pseudocount
)

# --- Highlight first seed ---
first_seed = long_df_combined["lineage_seed"].unique()[0]
long_df_combined["is_highlight"] = long_df_combined["lineage_seed"] == first_seed
background_df = long_df_combined[~long_df_combined["is_highlight"]]
highlight_df = long_df_combined[long_df_combined["is_highlight"]]

# Color encoding
molecules = long_df_combined["Molecule"].unique().tolist()
color_encoding = alt.Color(
    "Molecule:N",
    scale=alt.Scale(domain=molecules, scheme="tableau20"),
    legend=alt.Legend(title="Molecule", labelColor="black", titleColor="black"),
)

# --- Background lines ---
background_chart = (
    alt.Chart(background_df)
    .mark_line(opacity=0.25, strokeWidth=1.5)
    .encode(
        x=alt.X("Time (min):Q", axis=alt.Axis(title="Time (min)")),
        y=alt.Y("FoldChange_log:Q", title="Log10 Fold Change"),
        color=alt.Color("Molecule:N", scale=alt.Scale(scheme="tableau20"), legend=None),
        detail="lineage_seed:N",
    )
)

# --- Highlight line ---
highlight_chart = (
    alt.Chart(highlight_df)
    .mark_line(strokeWidth=4)
    .encode(x="Time (min):Q", y="FoldChange_log:Q", color=color_encoding)
)

# --- Default line (dotted) ---
default_chart = (
    alt.Chart(long_df_default_combined)
    .mark_line(strokeDash=[5, 5], strokeWidth=2)
    .encode(
        x="Time (min):Q",
        y="FoldChange_log:Q",
        color=alt.Color(
            "Molecule:N",
            scale=alt.Scale(domain=molecules, scheme="tableau20"),
            legend=None,
        ),
        detail="Molecule:N",
    )
)

# --- Events as gray bars ---
t_inhibit = 10 / 60
t_replication_done = 1327 / 60
events_df = pd.DataFrame(
    {
        "Time (min)": [t_inhibit, t_replication_done],
        "Event": ["Transcription Inhibition at 10s", "Ongoing Replication Completed"],
    }
)

vlines = (
    alt.Chart(events_df)
    .mark_rule(color="gray", opacity=0.3, size=4)
    .encode(x="Time (min):Q")
)

labels = (
    alt.Chart(events_df)
    .mark_text(
        angle=0, dy=-10, align="center", fontSize=14
    )  # put labels on top or bottom as needed
    .encode(
        x="Time (min):Q", y=alt.value(-10), text="Event:N", color=alt.value("black")
    )
)

# --- Combine ---
chart_log = (
    alt.layer(background_chart, default_chart, highlight_chart, vlines, labels)
    .configure_axis(grid=False, labelFontSize=14, titleFontSize=16)
    .configure_legend(labelFontSize=14, titleFontSize=16)
    .configure_title(fontSize=18, fontWeight="bold")
    .properties(
        width=800,
        height=500,
        title="Log10 Fold Change: Transcription Inhibition vs Baseline",
    )
)

# Save
chart_log.save(os.path.join(outdir, "runout_logfoldchange.html"))

# %%
long_df_combined["FoldChange_log2"] = np.log2(
    long_df_combined["FoldChange"] + pseudocount
)
long_df_default_combined["FoldChange_log2"] = np.log2(
    long_df_default_combined["FoldChange"] + pseudocount
)

# --- Highlight first seed ---
first_seed = long_df_combined["lineage_seed"].unique()[0]
long_df_combined["is_highlight"] = long_df_combined["lineage_seed"] == first_seed
background_df = long_df_combined[~long_df_combined["is_highlight"]]
highlight_df = long_df_combined[long_df_combined["is_highlight"]]

# Color encoding
molecules = long_df_combined["Molecule"].unique().tolist()
color_encoding = alt.Color(
    "Molecule:N",
    scale=alt.Scale(domain=molecules, scheme="tableau20"),
    legend=alt.Legend(title="Molecule", labelColor="black", titleColor="black"),
)

# --- Background lines ---
background_chart = (
    alt.Chart(background_df)
    .mark_line(opacity=0.25, strokeWidth=1.5)
    .encode(
        x=alt.X("Time (min):Q", axis=alt.Axis(title="Time (min)")),
        y=alt.Y("FoldChange_log2:Q", title="Log2 Fold Change"),
        color=alt.Color("Molecule:N", scale=alt.Scale(scheme="tableau20"), legend=None),
        detail="lineage_seed:N",
    )
)

# --- Highlight line ---
highlight_chart = (
    alt.Chart(highlight_df)
    .mark_line(strokeWidth=4)
    .encode(x="Time (min):Q", y="FoldChange_log2:Q", color=color_encoding)
)

# --- Default line (dotted) ---
default_chart = (
    alt.Chart(long_df_default_combined)
    .mark_line(strokeDash=[5, 5], strokeWidth=2)
    .encode(
        x="Time (min):Q",
        y="FoldChange_log2:Q",
        color=alt.Color(
            "Molecule:N",
            scale=alt.Scale(domain=molecules, scheme="tableau20"),
            legend=None,
        ),
        detail="Molecule:N",
    )
)

# --- Events as gray bars ---
t_inhibit = 10 / 60
t_replication_done = 1327 / 60
events_df = pd.DataFrame(
    {
        "Time (min)": [t_inhibit, t_replication_done],
        "Event": ["Transcription Inhibition at 10s", "Ongoing Replication Completed"],
    }
)

vlines = (
    alt.Chart(events_df)
    .mark_rule(color="gray", opacity=0.3, size=4)
    .encode(x="Time (min):Q")
)

labels = (
    alt.Chart(events_df)
    .mark_text(
        angle=0, dy=-10, align="center", fontSize=14
    )  # put labels on top or bottom as needed
    .encode(
        x="Time (min):Q", y=alt.value(-10), text="Event:N", color=alt.value("black")
    )
)

# --- Combine ---
chart_log = (
    alt.layer(background_chart, default_chart, highlight_chart, vlines, labels)
    .configure_axis(grid=False, labelFontSize=14, titleFontSize=16)
    .configure_legend(labelFontSize=14, titleFontSize=16)
    .configure_title(fontSize=18, fontWeight="bold")
    .properties(
        width=800,
        height=500,
        title="Log2 Fold Change: Transcription Inhibition vs Baseline",
    )
)

# Save
chart_log.save(os.path.join(outdir, "runout_log2foldchange.html"))

# %%
mrna_mass_sql = read_stacked_columns(
    history_sql,
    ["listeners__mass__mRna_mass"],
    order_results=False,
)
query_17 = f"""
SELECT lineage_seed, time, listeners__mass__mRna_mass FROM ({mrna_mass_sql})
ORDER BY lineage_seed, time
"""
output_df17 = conn.sql(query_17).df()
output_df17["Time (min)"] = output_df17["time"] / 60

first_seed = output_df17["lineage_seed"].unique()[0]
output_df17["is_highlight"] = output_df17["lineage_seed"] == first_seed

background_df = output_df17[~output_df17["is_highlight"]]
highlight_df = output_df17[output_df17["is_highlight"]]

# --- Background (blurred seeds) ---
background_chart = (
    alt.Chart(background_df)
    .mark_line(opacity=0.2, strokeWidth=1.5)
    .encode(
        x=alt.X("Time (min):Q", axis=alt.Axis(title="Time (min)", titlePadding=30)),
        y=alt.Y("listeners__mass__mRna_mass:Q", title="mRNA mass"),
        color=alt.value("gray"),
        detail="lineage_seed:N",
    )
)

# --- Highlighted seed ---
highlight_chart = (
    alt.Chart(highlight_df)
    .mark_line(strokeWidth=4, color="red")
    .encode(x="Time (min):Q", y="listeners__mass__mRna_mass:Q")
)

# --- Combine ---
chart = (
    (background_chart + highlight_chart)
    .configure_axis(grid=False, labelFontSize=14, titleFontSize=16)
    .properties(
        width=800, height=500, title="mRNA mass trajectories (first seed highlighted)"
    )
)

# Save
chart.save(os.path.join(outdir, "mrna_mass_multiseed.html"))

# %%
# Getting protein mass without ribosomes
ribosome_ids = bulk_name_to_idx(["CPLX0-3953[c]", "CPLX0-3962[c]"], bulk_molecule_ids)

bulk_data = sim_data.internal_state.bulk_molecules.bulk_data
# Access submass vector for specific rows
submass_5456 = bulk_data[5456]["mass"]
submass_5464 = bulk_data[5464]["mass"]

# Get 5th component (Python is 0-indexed → index 4)
val_5456 = submass_5456[5]
val_5464 = submass_5464[5]

print(val_5456)
print(val_5464)

# %%
AVOGADRO = N_A
mass_fg1 = (val_5456 / AVOGADRO) * 1e15
mass_fg2 = (val_5464 / AVOGADRO) * 1e15
ribosome_mass = mass_fg1 + mass_fg2
print(ribosome_mass)

# %%
ribosome_sql = read_stacked_columns(
    history_sql,
    ["listeners__unique_molecule_counts__active_ribosome"],
    order_results=False,
)
query_18 = f"""
SELECT time, lineage_seed, listeners__unique_molecule_counts__active_ribosome FROM ({ribosome_sql})
ORDER BY lineage_seed, time
"""
output_df18 = conn.sql(query_18).df()
# Convert time from seconds to minutes
output_df18["Time (min)"] = output_df18["time"] / 60

output_df18["Active Ribosome Mass (fg)"] = (
    output_df18["listeners__unique_molecule_counts__active_ribosome"] * ribosome_mass
)
# %%
# Free ribosome monomer masses
ribosome_monomer_ids = bulk_name_to_idx(
    [
        [
            [
                "EG10864-MONOMER[c]",
                "EG10865-MONOMER[c]",
                "EG10866-MONOMER[c]",
                "EG10867-MONOMER[c]",
                "EG10868-MONOMER[c]",
                "EG10869-MONOMER[c]",
                "EG10870-MONOMER[c]",
                "EG10871-MONOMER[c]",
                "EG10872-MONOMER[c]",
                "EG10873-MONOMER[c]",
                "EG10874-MONOMER[c]",
                "EG10875-MONOMER[c]",
                "EG10876-MONOMER[c]",
                "EG10877-MONOMER[c]",
                "EG10878-MONOMER[c]",
                "EG10879-MONOMER[c]",
                "EG10880-MONOMER[c]",
                "EG10881-MONOMER[c]",
                "EG10882-MONOMER[c]",
                "EG10883-MONOMER[c]",
                "EG10884-MONOMER[c]",
                "EG10885-MONOMER[c]",
                "EG10886-MONOMER[c]",
                "EG10887-MONOMER[c]",
                "EG10888-MONOMER[c]",
                "EG10889-MONOMER[c]",
                "EG10890-MONOMER[c]",
                "EG10891-MONOMER[c]",
                "EG10892-MONOMER[c]",
                "EG10900-MONOMER[c]",
                "EG10901-MONOMER[c]",
                "EG10902-MONOMER[c]",
                "EG10903-MONOMER[c]",
                "EG10904-MONOMER[c]",
                "EG10905-MONOMER[c]",
                "EG10906-MONOMER[c]",
                "EG10907-MONOMER[c]",
                "EG10908-MONOMER[c]",
                "EG10909-MONOMER[c]",
                "EG10910-MONOMER[c]",
                "EG10911-MONOMER[c]",
                "EG10912-MONOMER[c]",
                "EG10913-MONOMER[c]",
                "EG10914-MONOMER[c]",
                "EG10915-MONOMER[c]",
                "EG10916-MONOMER[c]",
                "EG10917-MONOMER[c]",
                "EG10918-MONOMER[c]",
                "EG10919-MONOMER[c]",
                "EG10920-MONOMER[c]",
                "EG11231-MONOMER[c]",
                "EG11232-MONOMER[c]",
                "EG50001-MONOMER[c]",
                "EG50002-MONOMER[c]",
            ]
        ]
    ],
    bulk_molecule_ids,
)
ribosome_monomer_masses = [bulk_data[i]["mass"] for i in ribosome_monomer_ids]
mass_array = np.array(ribosome_monomer_masses)
mass_array = mass_array.squeeze()
protein_mass_fg = np.array([m[5].asNumber() for m in mass_array]) / AVOGADRO * 1e15
# %%
ribosome_monomer_ids = ribosome_monomer_ids.squeeze()
inactive_ribosome_ids = ribosome_monomer_ids.tolist() + [5456, 5464]
ribosome_selected_cols = (
    ["Time (min)"]
    + ["lineage_seed"]
    + [bulk_df.columns[i + 2] for i in inactive_ribosome_ids]
)
ribosome_df = bulk_df[ribosome_selected_cols].copy()

# %%
count_cols = ribosome_df.columns[2:]  # skip Time and lineage_seed
counts_array = ribosome_df[count_cols].values  # shape (T, N)
protein_mass_fg = np.concatenate(
    [protein_mass_fg, [0.0005813422517875521, 0.0007908001340041741]]
)
counts_array = counts_array.astype(np.float64)
protein_mass_fg = np.array(protein_mass_fg, dtype=np.float64)
# %%
total_free_ribosome_mass = np.array(
    [np.dot(row, protein_mass_fg) for row in counts_array], dtype=np.float64
)
# %%
ribosome_df["Free Ribosome Mass (fg)"] = total_free_ribosome_mass
# Active ribosome mass
active_mass_df = output_df18[
    ["Time (min)", "lineage_seed", "Active Ribosome Mass (fg)"]
]

# Free ribosome mass
free_mass_df = ribosome_df[["Time (min)", "lineage_seed", "Free Ribosome Mass (fg)"]]
ribosome_mass_df = pd.merge(
    active_mass_df, free_mass_df, on=["Time (min)", "lineage_seed"]
)
ribosome_mass_df["Total Ribosome Mass (fg)"] = (
    ribosome_mass_df["Active Ribosome Mass (fg)"]
    + ribosome_mass_df["Free Ribosome Mass (fg)"]
)

# %%
protein_merged_df = pd.merge(
    output_df5[["Time (min)", "lineage_seed", "Protein_mass"]],
    ribosome_mass_df[["Time (min)", "lineage_seed", "Total Ribosome Mass (fg)"]],
    left_on=["Time (min)", "lineage_seed"],
    right_on=[ribosome_mass_df["Time (min)"], "lineage_seed"],
)
# %%
protein_merged_df["Protein Mass Excl. Ribosomes (fg)"] = (
    protein_merged_df["Protein_mass"] - protein_merged_df["Total Ribosome Mass (fg)"]
)
protein_mass_no_ribosomes_df = protein_merged_df[
    ["Time (min)", "lineage_seed", "Protein Mass Excl. Ribosomes (fg)"]
]
baseline = protein_mass_no_ribosomes_df[
    protein_mass_no_ribosomes_df["Time (min)"] == 0
][["lineage_seed", "Protein Mass Excl. Ribosomes (fg)"]].rename(
    columns={"Protein Mass Excl. Ribosomes (fg)": "baseline_mass"}
)
protein_mass_no_ribosomes_df = protein_mass_no_ribosomes_df.merge(
    baseline, on="lineage_seed", how="left"
)
protein_mass_no_ribosomes_df["Fold Change"] = (
    protein_mass_no_ribosomes_df["Protein Mass Excl. Ribosomes (fg)"]
    / protein_mass_no_ribosomes_df["baseline_mass"]
)
# %%

plt.figure(figsize=(8, 5))

for seed in protein_mass_no_ribosomes_df["lineage_seed"].unique():
    df_seed = protein_mass_no_ribosomes_df[
        protein_mass_no_ribosomes_df["lineage_seed"] == seed
    ]
    plt.plot(df_seed["Time (min)"], df_seed["Fold Change"], label=f"Lineage {seed}")

plt.xlabel("Time (min)")
plt.ylabel("Fold Change (Protein Mass excl. Ribosomes)")
plt.title("Fold Change of Protein Mass (excluding Ribosomes) over Time")
plt.legend()
plt.tight_layout()

# Save the figure
plt.savefig("protein_mass_fold_change.png", dpi=300)
