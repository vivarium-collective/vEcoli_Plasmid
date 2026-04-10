"""
Comparison of plasmid replication dynamics under two initial conditions:
  - "old": simulation starts with 1 pre-seeded replisome (b323d458)
  - "new": simulation starts with 0 replisomes (transformation experiment)

Run directly:
    python ecoli/analysis/multiseed/compare_plasmid_initial_conditions.py
"""

import os
import pandas as pd
import altair as alt

from ecoli.library.parquet_emitter import create_duckdb_conn, dataset_sql

OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "..", "out")
OLD_EXPERIMENT = "plasmid_multiseed_old"
NEW_EXPERIMENT = "plasmid_multiseed"
OLD_SINGLE_EXPERIMENT = "plasmid_multiseed_old_single"
NEW_SINGLE_EXPERIMENT = "mechanistic_chromosome_plasmid_priority_new"

COL_CHR_REPLISOME = "listeners__unique_molecule_counts__active_replisome"
COL_DOMAIN = "listeners__unique_molecule_counts__plasmid_domain"
OUTDIR = os.path.join(OUT_DIR, "plasmid_comparison")

COL_REPLISOME = "listeners__unique_molecule_counts__plasmid_active_replisome"
COL_PLASMID = "listeners__unique_molecule_counts__full_plasmid"


def load_data(conn, experiment_id, col):
    history_sql, _, _ = dataset_sql(OUT_DIR, [experiment_id])
    df = conn.sql(f"""
        SELECT time, lineage_seed, {col}
        FROM ({history_sql})
        ORDER BY lineage_seed, time ASC
    """).df()
    df["Time (min)"] = df["time"] / 60
    df["condition"] = experiment_id
    return df


def first_event_times(df, col, threshold=1):
    """Return time (min) when col first reaches >= threshold for each seed."""
    rows = []
    for seed, grp in df.groupby("lineage_seed"):
        grp = grp.sort_values("time")
        hit = grp[grp[col] >= threshold]
        t = float(hit["Time (min)"].iloc[0]) if len(hit) else float("nan")
        rows.append(
            {
                "lineage_seed": seed,
                "first_event_min": t,
                "condition": grp["condition"].iloc[0],
            }
        )
    return pd.DataFrame(rows)


def _closest_seed_to_mean(df, x_col, y_col):
    """Return the lineage_seed whose trajectory is closest to the mean trajectory."""
    mean_traj = df.groupby(x_col)[y_col].mean().rename("_mean")
    merged = df.merge(mean_traj, on=x_col)
    merged["_sq_err"] = (merged[y_col] - merged["_mean"]) ** 2
    return merged.groupby("lineage_seed")["_sq_err"].mean().idxmin()


def _single_trajectory_chart(df, col, y_title, condition_label):
    """Per-seed faded lines with the seed closest to the mean highlighted."""
    highlighted = _closest_seed_to_mean(df, "Time (min)", col)

    seed_lines = (
        alt.Chart(df[df["lineage_seed"] != highlighted])
        .mark_line(opacity=0.2, strokeWidth=1)
        .encode(
            x=alt.X("Time (min):Q"),
            y=alt.Y(f"{col}:Q", title=y_title),
            color=alt.Color("lineage_seed:N", legend=None),
        )
    )

    highlight_line = (
        alt.Chart(df[df["lineage_seed"] == highlighted])
        .mark_line(color="black", strokeWidth=3)
        .encode(
            x=alt.X("Time (min):Q"),
            y=alt.Y(f"{col}:Q"),
        )
    )

    return (seed_lines + highlight_line).properties(
        title=f"{condition_label} (highlighted seed: {highlighted})",
        width=500,
        height=350,
    )


def trajectory_chart(df_old, df_new, col, y_title, title):
    """Side-by-side trajectory charts for old and new conditions."""
    old_chart = _single_trajectory_chart(df_old, col, y_title, f"Old: {OLD_EXPERIMENT}")
    new_chart = _single_trajectory_chart(df_new, col, y_title, f"New: {NEW_EXPERIMENT}")
    return (
        alt.hconcat(old_chart, new_chart)
        .resolve_scale(y="shared")
        .properties(title=title)
    )


def histogram_chart(df_first, title, x_title):
    return (
        alt.Chart(df_first.dropna(subset=["first_event_min"]))
        .mark_bar(opacity=0.7)
        .encode(
            x=alt.X("first_event_min:Q", bin=alt.Bin(maxbins=20), title=x_title),
            y=alt.Y("count():Q", title="Seeds"),
            color=alt.Color(
                "condition:N",
                scale=alt.Scale(
                    domain=[OLD_EXPERIMENT, NEW_EXPERIMENT],
                    range=["#d62728", "#1f77b4"],
                ),
            ),
            column=alt.Column("condition:N", title="Initial condition"),
        )
        .properties(title=title, width=300, height=300)
    )


def variance_chart(df_old, df_new, col, y_title, title):
    """Standard deviation across seeds over time for both conditions."""
    rows = []
    for df in (df_old, df_new):
        std_df = df.groupby("Time (min)")[col].std().reset_index()
        std_df.columns = ["Time (min)", "std"]
        std_df["condition"] = df["condition"].iloc[0]
        rows.append(std_df)
    combined = pd.concat(rows)

    return (
        alt.Chart(combined)
        .mark_line()
        .encode(
            x=alt.X("Time (min):Q"),
            y=alt.Y("std:Q", title=y_title),
            color=alt.Color(
                "condition:N",
                scale=alt.Scale(
                    domain=[OLD_EXPERIMENT, NEW_EXPERIMENT],
                    range=["#d62728", "#1f77b4"],
                ),
            ),
        )
        .properties(title=title, width=600, height=350)
    )


def load_single_data(conn, experiment_id):
    """Load all relevant columns for a single simulation."""
    history_sql, _, _ = dataset_sql(OUT_DIR, [experiment_id])
    df = conn.sql(f"""
        SELECT time,
               {COL_REPLISOME},
               {COL_CHR_REPLISOME},
               {COL_PLASMID},
               {COL_DOMAIN}
        FROM ({history_sql})
        ORDER BY time ASC
    """).df()
    df["Time (min)"] = df["time"] / 60
    df["condition"] = experiment_id
    return df


def _single_line(df, col, y_title, color):
    return (
        alt.Chart(df)
        .mark_line(color=color, strokeWidth=2)
        .encode(
            x=alt.X("Time (min):Q"),
            y=alt.Y(f"{col}:Q", title=y_title),
        )
    )


def single_comparison(conn):
    """Side-by-side single-sim plots for old vs new initial conditions."""
    df_old = load_single_data(conn, OLD_SINGLE_EXPERIMENT)
    df_new = load_single_data(conn, NEW_SINGLE_EXPERIMENT)

    panels = [
        (COL_REPLISOME, "Plasmid active replisomes", "plasmid_replisome"),
        (COL_CHR_REPLISOME, "Chromosome active replisomes", "chr_replisome"),
        (COL_PLASMID, "Full plasmid count", "copy_number"),
        (COL_DOMAIN, "Plasmid domain count", "domain_count"),
    ]

    for col, y_title, fname in panels:
        old_chart = _single_line(df_old, col, y_title, "#d62728").properties(
            title=f"Old: {OLD_SINGLE_EXPERIMENT}", width=500, height=300
        )
        new_chart = _single_line(df_new, col, y_title, "#1f77b4").properties(
            title=f"New: {NEW_SINGLE_EXPERIMENT}", width=500, height=300
        )
        (
            alt.hconcat(old_chart, new_chart)
            .resolve_scale(y="shared")
            .properties(title=f"{y_title}: old vs new initial conditions")
        ).save(os.path.join(OUTDIR, f"single_{fname}_comparison.html"))


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    conn = create_duckdb_conn()

    # Load data
    df_rep_old = load_data(conn, OLD_EXPERIMENT, COL_REPLISOME)
    df_rep_new = load_data(conn, NEW_EXPERIMENT, COL_REPLISOME)
    df_plas_old = load_data(conn, OLD_EXPERIMENT, COL_PLASMID)
    df_plas_new = load_data(conn, NEW_EXPERIMENT, COL_PLASMID)

    # 1. Replisome trajectory comparison
    trajectory_chart(
        df_rep_old,
        df_rep_new,
        COL_REPLISOME,
        "Active plasmid replisomes",
        "Plasmid active replisomes: old vs new initial conditions",
    ).save(os.path.join(OUTDIR, "replisome_trajectory_comparison.html"))

    # 2. Plasmid copy number trajectory comparison
    trajectory_chart(
        df_plas_old,
        df_plas_new,
        COL_PLASMID,
        "Full plasmid count",
        "Plasmid copy number: old vs new initial conditions",
    ).save(os.path.join(OUTDIR, "copy_number_trajectory_comparison.html"))

    # 3. First initiation time histogram
    first_rep_old = first_event_times(df_rep_old, COL_REPLISOME, threshold=1)
    first_rep_new = first_event_times(df_rep_new, COL_REPLISOME, threshold=1)
    first_rep = pd.concat([first_rep_old, first_rep_new])

    histogram_chart(
        first_rep,
        "Distribution of first replication initiation time across seeds",
        "Time of first replisome assembly (min)",
    ).save(os.path.join(OUTDIR, "first_initiation_histogram.html"))

    # 4. First plasmid doubling time histogram
    first_plas_old = first_event_times(df_plas_old, COL_PLASMID, threshold=2)
    first_plas_new = first_event_times(df_plas_new, COL_PLASMID, threshold=2)
    first_plas = pd.concat([first_plas_old, first_plas_new])

    histogram_chart(
        first_plas,
        "Distribution of first plasmid doubling time across seeds",
        "Time of first copy number doubling (min)",
    ).save(os.path.join(OUTDIR, "first_doubling_histogram.html"))

    # 5. Cross-seed variance over time
    variance_chart(
        df_rep_old,
        df_rep_new,
        COL_REPLISOME,
        "Std dev of active plasmid replisomes",
        "Replisome count variability across seeds: old vs new",
    ).save(os.path.join(OUTDIR, "replisome_variance_comparison.html"))

    variance_chart(
        df_plas_old,
        df_plas_new,
        COL_PLASMID,
        "Std dev of plasmid copy number",
        "Copy number variability across seeds: old vs new",
    ).save(os.path.join(OUTDIR, "copy_number_variance_comparison.html"))

    # 6. Single simulation comparisons
    single_comparison(conn)

    print(f"Plots saved to {OUTDIR}/")


if __name__ == "__main__":
    main()
