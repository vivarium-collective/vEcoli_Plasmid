# %%
import json
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import altair as alt
import os

# from ecoli.analysis.single.test_analysis3 import sim_data
merged = []

with open(
    "/Users/rashmidissasekara/Documents/code/vEcoli/allocator_debug_out/plasmid_allocator3.jsonl"
) as f:
    lines = [json.loads(line) for line in f]

i = 0
while i < len(lines) - 1:
    a = lines[i]
    b = lines[i + 1]

    # Safety check
    assert a["source"] == "plasmid_replication"
    assert b["source"] == "allocator"

    merged.append(
        {
            "time": a["time"],
            "active replisomes": a["active replisomes"],
            "active replisome id": a["active replisome id"],
            "ready plasmids": a["ready plasmids"],
            "idle plasmids": a["idle plasmids"],
            "total plasmids": a["total plasmids"],
            "requested": b["requested"],
            "allocated": b["allocated"],
        }
    )

    i += 2

# %%


df = pd.DataFrame(merged)
# %%
rows = []

for _, row in df.iterrows():
    time = row["time"]

    for (bulk_id, req), (_, alloc) in zip(row["requested"], row["allocated"]):
        rows.append(
            {
                "time": time,
                "bulk_id": int(bulk_id),
                "requested": int(req),
                "allocated": int(alloc),
            }
        )

df_long = pd.DataFrame(rows)

# %%
replisome_ids = {5304, 5382, 5376, 6568, 7232, 7188}
dntp_ids = {6210, 6221, 6298, 12448}

df_replisome = df_long[df_long["bulk_id"].isin(replisome_ids)].copy()
df_dntp = df_long[df_long["bulk_id"].isin(dntp_ids)].copy()

# %%
common_names = {
    5376: "replicative DNA helicase (M)",
    6568: "DNA primase (M)",
    7232: "DNA polymerase III subunit δ' (M)",
    7188: "DNA polymerase III subunit δ (M)",
    5304: "DNA polymerase III, core enzyme (T)",
    5382: "β sliding clamp (T)",
}
# %%


plt.figure(figsize=(13, 10))
ax = plt.gca()
ax.yaxis.set_major_locator(mticker.MaxNLocator(integer=True, nbins=20))

colors = plt.cm.tab10.colors
id_to_color = dict(zip(replisome_ids, colors))

bulk_handles = []

for bulk_id in replisome_ids:
    sub_df = df_replisome[df_replisome["bulk_id"] == bulk_id]
    color = id_to_color[bulk_id]

    # vertical connectors
    for _, row in sub_df.iterrows():
        plt.plot(
            [row["time"], row["time"]],
            [row["allocated"], row["requested"]],
            color=color,
            alpha=0.3,
            linewidth=1,
        )

    # requested
    plt.scatter(sub_df["time"], sub_df["requested"], color=color, marker="o", s=30)

    # allocated
    plt.scatter(sub_df["time"], sub_df["allocated"], color=color, marker="x", s=30)

    # bulk legend handle (color only)
    name = common_names.get(bulk_id, f"Bulk {bulk_id}")

    bulk_handles.append(plt.Line2D([0], [0], color=color, lw=3, label=name))

marker_handles = [
    plt.Line2D([0], [0], color="black", marker="o", linestyle="", label="Requested"),
    plt.Line2D([0], [0], color="black", marker="x", linestyle="", label="Allocated"),
]

plt.legend(
    handles=bulk_handles + marker_handles,
    title="Replisome subunits",
    ncol=2,
    frameon=False,
    fontsize=9,
    title_fontsize=10,
)

plt.xlabel("Time (s)")
plt.ylabel("Subunit count")
plt.title("Requested vs Allocated Replisome Subunits Over Time")

plt.tight_layout()
plt.savefig("replisome_requested_allocated.png", dpi=300)
plt.close()


# %%
# Extracting plasmid time series
total_plasmids = df["total plasmids"].to_numpy()
df["idle_plasmid_count"] = df["idle plasmids"].apply(len)
idle_plasmids = df["idle_plasmid_count"].to_numpy()
df["ready_plasmid_count"] = df["ready plasmids"].apply(len)
ready_plasmids = df["ready_plasmid_count"].to_numpy()
active_replisomes = df["active replisomes"].to_numpy()
times = df["time"].to_numpy()

# two vertically stacked plots (shared x-axis)
fig, (ax1, ax2) = plt.subplots(
    nrows=2,
    ncols=1,
    figsize=(13, 12),
    sharex=True,
    gridspec_kw={"height_ratios": [2, 1]},
)

for bulk_id in replisome_ids:
    sub_df = df_replisome[df_replisome["bulk_id"] == bulk_id]
    color = id_to_color[bulk_id]
    # vertical connectors
    for _, row in sub_df.iterrows():
        ax1.plot(
            [row["time"], row["time"]],
            [row["allocated"], row["requested"]],
            color=color,
            alpha=0.3,
            linewidth=1,
        )

    ax1.scatter(sub_df["time"], sub_df["requested"], color=color, marker="o", s=30)
    ax1.scatter(sub_df["time"], sub_df["allocated"], color=color, marker="x", s=30)

ax1.set_ylabel("Molecules")
ax1.set_title("Allocator: Requested vs Allocated")
ax1.legend(
    handles=bulk_handles + marker_handles,
    title="Replisome subunits",
    ncol=2,
    frameon=False,
    fontsize=9,
    title_fontsize=10,
)

ax2.plot(times, total_plasmids, label="Total plasmids", linewidth=2)
ax2.plot(times, idle_plasmids, label="Idle plasmids", linewidth=2)
ax2.plot(times, ready_plasmids, label="Ready plasmids", linewidth=2)
ax2.plot(times, active_replisomes, label="Active replisomes", linewidth=2)

ax2.set_ylabel("Count")
ax2.set_xlabel("Time (s)")
ax2.set_title("Plasmid replication state")
ax2.legend(ncol=3, frameon=False)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("allocator_vs_plasmid_states.png", dpi=300)
plt.close()

# %%
print(type(sub_df.iloc[0]["allocated"]), sub_df.iloc[0]["allocated"])

# %%
debug_outdir = "allocator_debug_out"
json_file = os.path.join(debug_outdir, "chromosome_plasmid_allocator.jsonl")

common_names = {
    5376: "replicative DNA helicase (M)",
    6568: "DNA primase (M)",
    7232: "DNA polymerase III δ' (M)",
    7188: "DNA polymerase III δ (M)",
    5304: "DNA polymerase III core (T)",
    5382: "β sliding clamp (T)",
    6210: "DATP[c]",
    6221: "DCTP[c]",
    6298: "DGTP[c]",
    12448: "TTP[c]",
}


def get_common_name(idx):
    return common_names.get(idx, f"id_{idx}")


records = []
with open(json_file) as f:
    for line in f:
        records.append(json.loads(line))

top3_records = records[:3]
times = [1756, 1757, 1758]

data = []

for record, t in zip(top3_records, times):
    for proc in ["ecoli-chromosome-replication", "ecoli-plasmid-replication"]:
        process_name = "Chromosome" if "chromosome" in proc else "Plasmid"

        req = {idx: val for idx, val in record[proc]["requested"]}
        alloc = {idx: val for idx, val in record[proc]["allocated"]}

        all_ids = set(req) | set(alloc)

        for idx in all_ids:
            data.append(
                {
                    "Time": t,
                    "Process": process_name,
                    "Molecule": get_common_name(idx),
                    "Type": "Requested",
                    "Count": req.get(idx, 0),
                }
            )

            data.append(
                {
                    "Time": t,
                    "Process": process_name,
                    "Molecule": get_common_name(idx),
                    "Type": "Allocated",
                    "Count": alloc.get(idx, 0),
                }
            )

df = pd.DataFrame(data)

df["Time_min"] = df["Time"] / 60

df["TimeLabel"] = df["Time_min"].round(2).astype(str) + " min"

# label the initiation time
df.loc[df["Time"] == 1758, "TimeLabel"] = (
    f"{round(1758 / 60, 2)} min – Chromosome initiated"
)

ratio_df = df.pivot_table(
    index=["Process", "TimeLabel", "Molecule"],
    columns="Type",
    values="Count",
    aggfunc="sum",
).reset_index()

ratio_df["ratio"] = ratio_df["Allocated"] / ratio_df["Requested"]

chart = (
    alt.Chart(ratio_df)
    .mark_bar()
    .encode(
        x=alt.X("ratio:Q", title=None, scale=alt.Scale(domain=[0, 1.1])),
        y=alt.Y("Molecule:N", sort="-x", title=None),
        color=alt.Color(
            "ratio:Q",
            scale=alt.Scale(scheme="redyellowgreen"),
            legend=alt.Legend(title=["Allocation ratio", "(Allocated / Requested)"]),
        ),
        column=alt.Column("TimeLabel:N", title="Time (min)"),
        row=alt.Row("Process:N", title=None),
        tooltip=["Process", "TimeLabel", "Molecule", "Requested", "Allocated", "ratio"],
    )
    .properties(width=230, height=180)
)
os.makedirs(debug_outdir, exist_ok=True)

html_path = os.path.join(debug_outdir, "replication_resource_allocation_ratio.html")

chart.save(html_path)

print("Figure saved to:", html_path)

# %%
# Stacked figure
# Load dataframe
df_plot = pd.read_csv(
    "out/mechanistic_chromosome_plasmid_2/analyses/variant=0/lineage_seed=0/generation=1/agent_id=0/plots/replisome_dynamics_times.csv"
)

# Keep full window for line continuity
df_plot["Time (min)"] = (df_plot["time"] / 60).round(2)

# Clean x-axis labels
df_plot["XLabel"] = df_plot["Time (min)"].round(2).astype(str)
# The three main points as strings matching XLabel
x_ticks = df_plot.loc[df_plot["time"].isin([1756, 1757, 1758]), "XLabel"].tolist()

# Tooltip keeps extra info
df_plot["TooltipLabel"] = df_plot["Time (min)"].round(2).astype(str)
df_plot.loc[df_plot["time"] == 1758, "TooltipLabel"] = (
    f"{round(1758 / 60, 2)} min – Chromosome initiated"
)

# Define categorical x-axis with all times in the window
df_plot["TimeLabel"] = df_plot["Time (min)"].astype(str)

# Top panel: existing chart
chart_top = chart

# Bottom panel: line chart
chart_bottom = (
    alt.Chart(df_plot)
    .mark_line(point=True)
    .encode(
        x=alt.X(
            "XLabel:N",
            title="Time (min)",
            axis=alt.Axis(
                values=x_ticks,  # only show ticks for main points
                labelAngle=0,
            ),
        ),
        y=alt.Y("Active replisomes:Q", title="Active replisomes"),
        color=alt.Color("Replication:N", title="Replication"),
        tooltip=["Replication", "Active replisomes", "TooltipLabel"],
    )
    .properties(
        width=600, height=250, title="Replisome dynamics around chromosome initiation"
    )
)

# Combine
final_chart = alt.vconcat(chart_top, chart_bottom)

# Save
html_path2 = os.path.join(debug_outdir, "combined_top_bottom_panel.html")
final_chart.save(html_path2)

# %%
debug_outdir = "allocator_debug_out"
json_file = os.path.join(debug_outdir, "chromosome_plasmid_allocator_no_priority.jsonl")

# Short molecule names
common_names = {
    5376: "replicative DNA helicase (M)",
    6568: "DNA primase (M)",
    7232: "DNA polymerase III δ' (M)",
    7188: "DNA polymerase III δ (M)",
    5304: "DNA polymerase III core (T)",
    5382: "β sliding clamp (T)",
    6210: "dATP",
    6221: "dCTP",
    6298: "dGTP",
    12448: "dTTP",
}


def get_name(idx):
    return common_names.get(idx, f"id_{idx}")


# Load JSONL records
records = [json.loads(line) for line in open(json_file)]
top3_records = records[:3]
times = [1756, 1757, 1758]

data = []
for record, t in zip(top3_records, times):
    for proc in ["ecoli-chromosome-replication", "ecoli-plasmid-replication"]:
        process_name = "Chromosome" if "chromosome" in proc else "Plasmid"
        req = {idx: val for idx, val in record[proc]["requested"]}
        alloc = {idx: val for idx, val in record[proc]["allocated"]}
        for idx in set(req) | set(alloc):
            requested = req.get(idx, 0)
            allocated = alloc.get(idx, 0)
            ratio = allocated / requested if requested > 0 else 0
            data.append(
                {
                    "Time": t,
                    "Process": process_name,
                    "Molecule": get_name(idx),
                    "Requested": requested,
                    "Allocated": allocated,
                    "ratio": ratio,
                }
            )

df = pd.DataFrame(data)

# Convert time to minutes
df["Time_min"] = df["Time"] / 60
df["TimeLabel"] = df["Time_min"].round(2).astype(str) + " min"
df.loc[df["Time"] == 1758, "TimeLabel"] = (
    f"{round(1758 / 60, 2)} min – Chromosome initiated"
)

subunit_order = [
    "replicative DNA helicase (M)",
    "DNA primase (M)",
    "DNA polymerase III δ' (M)",
    "DNA polymerase III δ (M)",
    "DNA polymerase III core (T)",
    "β sliding clamp (T)",
    "dATP",
    "dCTP",
    "dGTP",
    "dTTP",
]

# --------------------------
# Layered bars: requested vs allocated
# --------------------------
layered_bars = alt.layer(
    alt.Chart(df)
    .mark_bar(color="lightgray")
    .encode(
        x=alt.X("Requested:Q", title=None),
        y=alt.Y("Molecule:N", sort=subunit_order, title=None),
    ),
    alt.Chart(df)
    .mark_bar()
    .encode(
        x=alt.X("Allocated:Q", title=None),
        y=alt.Y("Molecule:N", sort=subunit_order),
        color=alt.Color(
            "ratio:Q",
            scale=alt.Scale(scheme="redyellowgreen", domain=[0, 1]),
            legend=alt.Legend(title=["Allocation ratio", "(Allocated / Requested)"]),
        ),
    ),
)

# Facet by Time and Process
chart = layered_bars.facet(
    row=alt.Row(
        "Process:N", title=None, header=alt.Header(labelFontSize=14, titleFontSize=14)
    ),
    column=alt.Column(
        "TimeLabel:N",
        title="Time (min)",
        header=alt.Header(labelFontSize=12, titleFontSize=14),
    ),
).resolve_scale(x="independent")

# Configure y-axis (subunit names) font
chart = chart.configure_axis(
    labelFontSize=12,  # y-axis labels (molecule names)
    titleFontSize=14,
)

# Add bottom title
chart = chart.configure_title(anchor="middle").properties(
    title=alt.TitleParams(
        "Requested molecule counts",
        orient="bottom",  # place title at bottom
        anchor="middle",  # center horizontally
        fontSize=14,
        offset=20,
    )
)

os.makedirs(debug_outdir, exist_ok=True)
html_path = os.path.join(
    debug_outdir, "replication_resource_allocation_ratio_no_priority.html"
)
chart.save(html_path)

print("Figure saved to:", html_path)

# %%
debug_outdir = "allocator_debug_out"

# --------------------------
# Top panel: Allocation ratios
# --------------------------
json_file = os.path.join(debug_outdir, "chromosome_plasmid_allocator_no_priority.jsonl")

# Short molecule names
common_names = {
    5376: "Helicase",
    6568: "Primase",
    7232: "Pol III δ′",
    7188: "Pol III δ",
    5304: "Pol III core",
    5382: "β clamp",
    6210: "dATP",
    6221: "dCTP",
    6298: "dGTP",
    12448: "dTTP",
}


def get_name(idx):
    return common_names.get(idx, f"id_{idx}")


# Load JSONL records
records = [json.loads(line) for line in open(json_file)]
top3_records = records[:3]
times = [1756, 1757, 1758]

data = []
for record, t in zip(top3_records, times):
    for proc in ["ecoli-chromosome-replication", "ecoli-plasmid-replication"]:
        process_name = "Chromosome" if "chromosome" in proc else "Plasmid"
        req = {idx: val for idx, val in record[proc]["requested"]}
        alloc = {idx: val for idx, val in record[proc]["allocated"]}
        for idx in set(req) | set(alloc):
            requested = req.get(idx, 0)
            allocated = alloc.get(idx, 0)
            ratio = allocated / requested if requested > 0 else 0
            data.append(
                {
                    "Time": t,
                    "Process": process_name,
                    "Molecule": get_name(idx),
                    "Requested": requested,
                    "Allocated": allocated,
                    "ratio": ratio,
                }
            )

df = pd.DataFrame(data)

# Convert time to minutes
df["Time_min"] = df["Time"] / 60
df["TimeLabel"] = df["Time_min"].round(2).astype(str) + " min"
df.loc[df["Time"] == 1758, "TimeLabel"] = (
    f"{round(1758 / 60, 2)} min – Chromosome initiated"
)

subunit_order = [
    "Helicase",
    "Primase",
    "Pol III core",
    "Pol III δ",
    "Pol III δ′",
    "β clamp",
    "dATP",
    "dCTP",
    "dGTP",
    "dTTP",
]

# Layered bars: requested vs allocated
top_layer = alt.layer(
    # full requested bar
    alt.Chart(df)
    .mark_bar(color="lightgray")
    .encode(
        x=alt.X("Requested:Q", title=None),
        y=alt.Y("Molecule:N", sort=subunit_order, title=None),
    ),
    # allocated overlay
    alt.Chart(df)
    .mark_bar()
    .encode(
        x=alt.X("Allocated:Q", title=None),
        y=alt.Y("Molecule:N", sort=subunit_order),
        color=alt.Color(
            "ratio:Q",
            scale=alt.Scale(scheme="redyellowgreen", domain=[0, 1]),
            legend=alt.Legend(title=["Allocation ratio", "(Allocated / Requested)"]),
        ),
    ),
)

# Facet the top layered chart
chart_top = top_layer.facet(
    row=alt.Row("Process:N", title=None),
    column=alt.Column("TimeLabel:N", title="Time (min)"),
).resolve_scale(x="independent")

# --------------------------
# Bottom panel: replisomes
# --------------------------
df_plot = pd.read_csv(
    "out/mechanistic_chromosome_plasmid_2/analyses/variant=0/lineage_seed=0/generation=1/agent_id=0/plots/replisome_dynamics_times.csv"
)

df_plot["Time (min)"] = (df_plot["time"] / 60).round(2)
df_plot["XLabel"] = df_plot["Time (min)"].round(2).astype(str)

# The three main points
x_ticks = df_plot.loc[df_plot["time"].isin([1756, 1757, 1758]), "XLabel"].tolist()

# Tooltip for extra info
df_plot["TooltipLabel"] = df_plot["Time (min)"].round(2).astype(str)
df_plot.loc[df_plot["time"] == 1758, "TooltipLabel"] = (
    f"{round(1758 / 60, 2)} min – Chromosome initiated"
)

# Bottom line chart
chart_bottom = (
    alt.Chart(df_plot)
    .mark_line(point=True)
    .encode(
        x=alt.X(
            "XLabel:N",
            title="Time (min)",
            axis=alt.Axis(
                values=x_ticks,  # only show main points
                labelAngle=0,
            ),
        ),
        y=alt.Y("Active replisomes:Q", title="Active replisomes"),
        color=alt.Color("Replication:N", title="Replication"),
        tooltip=["Replication", "Active replisomes", "TooltipLabel"],
    )
    .properties(
        width=850, height=250, title="Replisome dynamics around chromosome initiation"
    )
)
bottom_label = (
    alt.Chart(pd.DataFrame({"text": ["Requested molecule counts"]}))
    .mark_text(align="center", baseline="bottom", fontSize=14)
    .encode(
        text="text:N",
        x=alt.value(300),
        y=alt.value(0),  # center horizontally, roughly half of top chart width
    )
    .properties(
        width=850,
        height=0.5,
    )
)

# --------------------------
# Combine everything vertically
# --------------------------
final_chart = alt.vconcat(chart_top, bottom_label, chart_bottom).resolve_scale(
    x="independent"
)

# Save
os.makedirs(debug_outdir, exist_ok=True)
html_path = os.path.join(debug_outdir, "combined_top_bottom_panel_no_priority.html")
final_chart.save(html_path)
print("Figure saved to:", html_path)

# %%

html_path_bottom = os.path.join(debug_outdir, "replication_dynamics_no_priority.html")
chart_bottom.save(html_path_bottom)
print("Bottom panel saved to:", html_path_bottom)
