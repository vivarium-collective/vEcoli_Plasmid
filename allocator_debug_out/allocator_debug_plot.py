# %%
import json
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

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
