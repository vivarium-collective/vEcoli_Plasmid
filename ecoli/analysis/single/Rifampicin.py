# %%
import os
import duckdb

import altair as alt
import pandas as pd


conn = duckdb.connect()
query_dict = {
    "experiment_id": "rifampicin-exp",
    "variant": 0,
    "lineage_seed": 0,
    "generation": 1,
}
history_sql = f"""
SELECT * FROM read_parquet('out/{query_dict["experiment_id"]}/history/*/*/*/*/*/*.pq', hive_partitioning=true)
"""
query_1 = f"""
SELECT time, listeners__unique_molecule_counts__active_replisome FROM ({history_sql})
ORDER BY time ASC
"""
output_df8 = conn.sql(query_1).df()
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

outdir = "out/manual_rifampicin_analysis"
os.makedirs(outdir, exist_ok=True)

html_path8 = os.path.join(outdir, "activereplisomecounts.html")
chart8.save(html_path8)
# %%
query_2 = f"""
SELECT time, listeners__mass__cell_mass FROM ({history_sql})
ORDER BY time ASC
"""
output_df2 = conn.sql(query_2).df()
# Convert time from seconds to minutes
output_df2["Time (min)"] = output_df2["time"] / 60

query_3 = f"""
SELECT time, listeners__mass__dna_mass FROM ({history_sql})
ORDER BY time ASC
"""

output_df3 = conn.sql(query_3).df()

query_4 = f"""
SELECT time, listeners__mass__rna_mass FROM ({history_sql})
ORDER BY time ASC
"""
output_df4 = conn.sql(query_4).df()

query_5 = f"""
SELECT time, listeners__mass__protein_mass FROM ({history_sql})
ORDER BY time ASC
"""
output_df5 = conn.sql(query_5).df()

query_6 = f"""
SELECT time, listeners__mass__smallMolecule_mass FROM ({history_sql})
ORDER BY time ASC
"""
output_df6 = conn.sql(query_6).df()
# %%
for df in [output_df2, output_df3, output_df4, output_df5, output_df6]:
    df["Time (min)"] = df["time"] / 60

# Step 2: Rename the mass columns for clarity
output_df2 = output_df2.rename(columns={"listeners__mass__cell_mass": "Cell_mass"})
output_df3 = output_df3.rename(columns={"listeners__mass__dna_mass": "DNA_mass"})
output_df4 = output_df4.rename(columns={"listeners__mass__rna_mass": "RNA_mass"})
output_df5 = output_df5.rename(
    columns={"listeners__mass__protein_mass": "Protein_mass"}
)
output_df6 = output_df6.rename(
    columns={"listeners__mass__smallMolecule_mass": "SmallMolecule_mass"}
)

# Step 3: Merge all dataframes on 'Time (min)'
merged_df = (
    output_df2[["Time (min)", "Cell_mass"]]
    .merge(output_df3[["Time (min)", "DNA_mass"]], on="Time (min)")
    .merge(output_df4[["Time (min)", "RNA_mass"]], on="Time (min)")
    .merge(output_df5[["Time (min)", "Protein_mass"]], on="Time (min)")
    .merge(output_df6[["Time (min)", "SmallMolecule_mass"]], on="Time (min)")
)

# Step 4 (optional): Convert to long format for plotting multiple molecules
long_df = merged_df.melt(id_vars=["Time (min)"], var_name="Molecule", value_name="Mass")

t_inhibit = 10 / 60
t_replication_done = 1327 / 60

events_df = pd.DataFrame(
    {
        "Time (min)": [t_inhibit, t_replication_done],
        "Event": ["Transcription inhibition at 10s", "Ongoing Replication complete"],
    }
)

vlines = (
    alt.Chart(events_df)
    .mark_rule(strokeDash=[5, 5], size=2, color="black")
    .encode(
        x="Time (min):Q",
    )
)

labels = (
    alt.Chart(events_df)
    .mark_text(angle=0, dy=15, align="center", fontSize=14)
    .encode(
        x="Time (min):Q", y=alt.value(520), text="Event:N", color=alt.value("black")
    )
)

base_chart = (
    alt.Chart(long_df)
    .mark_line()
    .encode(
        x=alt.X("Time (min):Q", title="Time (min)", axis=alt.Axis(titlePadding=30)),
        y=alt.Y("Mass:Q", title="Mass (fg)", scale=alt.Scale(type="log")),
        color=alt.Color("Molecule:N", legend=alt.Legend(title="Molecule")),
    )
    .properties(
        width=800,
        height=500,
        title="Mass Trajectories with Transcription Inhibition at 10 s",
    )
)


chart = (
    (base_chart + vlines + labels)
    .configure_axis(labelFontSize=12, titleFontSize=14, grid=False)
    .configure_legend(labelFontSize=14, titleFontSize=16)
    .configure_title(fontSize=18, fontWeight="bold")
)

# Save
html_path = os.path.join(outdir, "mass_components_transcription_inhibition.html")
chart.save(html_path)

# %%
output_df3["Time (min)"] = output_df3["time"] / 60

# Event times (in minutes)
t_inhibit = 10 / 60
t_replication_done = 1327 / 60

base_chart = (
    alt.Chart(output_df3)
    .mark_line()
    .encode(
        x=alt.X("Time (min):Q", title="Time (min)"),
        y=alt.Y(
            "DNA_mass:Q",
            title="Mass (fg)",
            scale=alt.Scale(type="log"),
        ),
    )
)

# Create vertical lines dataframe
events_df = pd.DataFrame(
    {
        "Time (min)": [t_inhibit, t_replication_done],
        "Event": ["Transcription inhibition", "Replication complete"],
    }
)

# Vertical lines
vlines = (
    alt.Chart(events_df)
    .mark_rule(strokeDash=[5, 5], size=2)
    .encode(
        x="Time (min):Q", color=alt.Color("Event:N", legend=alt.Legend(title="Event"))
    )
)

# Optional: add text labels near the lines
labels = (
    alt.Chart(events_df)
    .mark_text(angle=90, dx=5, dy=-5)
    .encode(x="Time (min):Q", text="Event:N", color="Event:N")
)

# Combine everything
chart = (base_chart + vlines + labels).properties(
    title="DNA mass trajectory upon transcription inhibition at 10s",
    width=600,
    height=400,
)
# %%
outdir = "out/manual_rifampicin_analysis"
os.makedirs(outdir, exist_ok=True)
# %%
html_path8 = os.path.join(outdir, "dnamass.html")
chart.save(html_path8)

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
long_df_default = merged_df.melt(
    id_vars=["Time (min)"], var_name="Molecule", value_name="Mass"
)

# %%
color_encoding = alt.Color(
    "Molecule:N",
    scale=alt.Scale(scheme="tableau20"),
    legend=alt.Legend(title="Molecule", labelColor="black", titleColor="black"),
)

background_chart = (
    alt.Chart(long_df_default)
    .mark_line(strokeDash=[5, 5], strokeWidth=2)
    .encode(
        x=alt.X("Time (min):Q", title="Time (min)", axis=alt.Axis(titlePadding=30)),
        y=alt.Y("Mass:Q", scale=alt.Scale(type="log")),
        color=alt.Color("Molecule:N", scale=alt.Scale(scheme="tableau20"), legend=None),
    )
)

highlight_chart = (
    alt.Chart(long_df)
    .mark_line(strokeWidth=3)
    .encode(
        x="Time (min):Q",
        y=alt.Y("Mass:Q", scale=alt.Scale(type="log")),
        color=color_encoding,
    )
)

chart = (
    (background_chart + highlight_chart + vlines + labels)
    .configure_axis(grid=False, labelFontSize=14, titleFontSize=16)
    .configure_legend(labelFontSize=14, titleFontSize=16)
    .configure_title(fontSize=18, fontWeight="bold")
    .properties(
        width=800,
        height=500,
        padding={"bottom": 60},
        title="Mass trajectories: transcription inhibition vs baseline",
    )
)

html_path8 = os.path.join(outdir, "runoutsynthesis.html")
chart.save(html_path8)

# %%
# --- Separate Legend Chart with tableau20 ---
legend_chart = (
    alt.Chart(long_df)
    .mark_line(strokeWidth=4)  # thick lines for legend
    .encode(
        color=color_encoding  # REUSE the same color_encoding from highlight_chart
    )
    .transform_aggregate(
        count="count()", groupby=["Molecule"]
    )  # one entry per molecule
    .encode(y=alt.value(0))  # collapse all lines into a single row
    .properties(width=200, height=200)
    .configure_axis(grid=False, domain=False, ticks=False, labels=False)
)

# Save the legend separately
legend_chart.save(os.path.join(outdir, "legend_only_runout.html"))
