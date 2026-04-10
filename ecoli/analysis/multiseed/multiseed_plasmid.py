import os
from typing import Any
import pandas as pd
import altair as alt

from duckdb import DuckDBPyConnection


def _closest_seed_to_mean(df, x_col, y_col):
    """Return the lineage_seed whose trajectory is closest to the mean trajectory."""
    mean_traj = df.groupby(x_col)[y_col].mean().rename("_mean")
    merged = df.merge(mean_traj, on=x_col)
    merged["_sq_err"] = (merged[y_col] - merged["_mean"]) ** 2
    return merged.groupby("lineage_seed")["_sq_err"].mean().idxmin()


def _make_multiseed_chart(
    df, x_col, y_col, title, x_title, y_title, highlighted_seed, width=600, height=400
):
    """Layer faded per-seed lines with highlighted_seed drawn on top."""
    seed_lines = (
        alt.Chart(df[df["lineage_seed"] != highlighted_seed])
        .mark_line(opacity=0.2, strokeWidth=1)
        .encode(
            x=alt.X(f"{x_col}:Q", title=x_title),
            y=alt.Y(f"{y_col}:Q", title=y_title),
            color=alt.Color("lineage_seed:N", legend=None),
        )
    )

    highlight_line = (
        alt.Chart(df[df["lineage_seed"] == highlighted_seed])
        .mark_line(color="black", strokeWidth=3)
        .encode(
            x=alt.X(f"{x_col}:Q", title=x_title),
            y=alt.Y(f"{y_col}:Q", title=y_title),
        )
    )

    full_title = f"{title} (highlighted seed: {highlighted_seed})"
    return (seed_lines + highlight_line).properties(
        title=full_title, width=width, height=height
    )


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
    # === Full plasmid counts ===
    df_plasmid = conn.sql(f"""
        SELECT time, lineage_seed, listeners__unique_molecule_counts__full_plasmid
        FROM ({history_sql})
        ORDER BY lineage_seed, time ASC
    """).df()
    df_plasmid["Time (min)"] = df_plasmid["time"] / 60

    # Pick the single highlighted seed from the full plasmid trajectory — used in all plots
    highlighted_seed = _closest_seed_to_mean(
        df_plasmid, "Time (min)", "listeners__unique_molecule_counts__full_plasmid"
    )

    _make_multiseed_chart(
        df_plasmid,
        x_col="Time (min)",
        y_col="listeners__unique_molecule_counts__full_plasmid",
        title="Full plasmid counts",
        x_title="Time (min)",
        y_title="Full plasmid counts",
        highlighted_seed=highlighted_seed,
    ).save(os.path.join(outdir, "full plasmid counts.html"))

    # First 60 seconds
    df_plasmid_60s = df_plasmid[df_plasmid["time"] <= 60].copy()
    df_plasmid_60s["Time (sec)"] = df_plasmid_60s["time"]

    _make_multiseed_chart(
        df_plasmid_60s,
        x_col="Time (sec)",
        y_col="listeners__unique_molecule_counts__full_plasmid",
        title="Full plasmid counts (first 60 seconds)",
        x_title="Time (sec)",
        y_title="Full plasmid counts",
        highlighted_seed=highlighted_seed,
    ).save(os.path.join(outdir, "full plasmid counts (60 seconds).html"))

    # === Plasmid active replisome counts ===
    df_replisome = conn.sql(f"""
        SELECT time, lineage_seed, listeners__unique_molecule_counts__plasmid_active_replisome
        FROM ({history_sql})
        ORDER BY lineage_seed, time ASC
    """).df()
    df_replisome["Time (min)"] = df_replisome["time"] / 60

    _make_multiseed_chart(
        df_replisome,
        x_col="Time (min)",
        y_col="listeners__unique_molecule_counts__plasmid_active_replisome",
        title="Plasmid active replisome counts",
        x_title="Time (min)",
        y_title="Counts",
        highlighted_seed=highlighted_seed,
    ).save(os.path.join(outdir, "plasmid_active_replisomes.html"))

    # First 60 seconds
    df_replisome_60s = df_replisome[df_replisome["time"] <= 60].copy()
    df_replisome_60s["Time (sec)"] = df_replisome_60s["time"]

    _make_multiseed_chart(
        df_replisome_60s,
        x_col="Time (sec)",
        y_col="listeners__unique_molecule_counts__plasmid_active_replisome",
        title="Plasmid active replisome counts (first 60 seconds)",
        x_title="Time (sec)",
        y_title="Counts",
        highlighted_seed=highlighted_seed,
    ).save(os.path.join(outdir, "plasmid active replisome counts (60 seconds).html"))

    # === Tight windows around chromosome initiation and completion ===
    df_chr = conn.sql(f"""
        SELECT time, lineage_seed, listeners__unique_molecule_counts__active_replisome
        FROM ({history_sql})
        ORDER BY lineage_seed, time ASC
    """).df()
    y_col_c = "listeners__unique_molecule_counts__active_replisome"
    y_col_p = "listeners__unique_molecule_counts__plasmid_active_replisome"

    # Derive event times from the highlighted seed's chromosome trajectory
    seed_chr = df_chr[df_chr["lineage_seed"] == highlighted_seed].sort_values("time")
    chromosome_completion_s = float(seed_chr[seed_chr[y_col_c] == 0]["time"].min())
    chromosome_initiation_s = float(
        seed_chr[
            (seed_chr["time"] > chromosome_completion_s) & (seed_chr[y_col_c] == 4)
        ]["time"].min()
    )

    def _make_tight_window_charts(event_time_s, event_label, outfile_prefix):
        """Build plasmid + chromosome tight-window plots around a single event time."""
        event_min = round(event_time_s / 60, 2)
        window = 30
        t_start = event_time_s - window
        t_end = event_time_s + window
        m_start = t_start / 60
        m_end = t_end / 60
        ticks = sorted({round(m_start, 2), round(m_end, 2), event_min})

        vline = (
            alt.Chart(pd.DataFrame({"x": [event_min]}))
            .mark_rule(color="red", strokeDash=[5, 5])
            .encode(x="x:Q", tooltip=alt.Tooltip(["x"], format=".2f"))
        )

        # --- Plasmid panel ---
        df_p_win = df_replisome[
            (df_replisome["time"] >= t_start) & (df_replisome["time"] <= t_end)
        ].copy()
        df_p_win["Time (min)"] = df_p_win["time"] / 60

        x_labeled = alt.X(
            "Time (min):Q",
            title="Time (min)",
            scale=alt.Scale(domain=(m_start, m_end), zero=False),
            axis=alt.Axis(
                values=ticks,
                labelExpr=f'datum.value == {event_min} ? "{event_label} at {event_min} min" : datum.value',
            ),
        )
        (
            alt.Chart(df_p_win[df_p_win["lineage_seed"] != highlighted_seed])
            .mark_line(opacity=0.2, strokeWidth=1)
            .encode(
                x=x_labeled,
                y=alt.Y(f"{y_col_p}:Q", title="Active plasmid replisomes"),
                color=alt.Color("lineage_seed:N", legend=None),
            )
            + alt.Chart(df_p_win[df_p_win["lineage_seed"] == highlighted_seed])
            .mark_line(color="black", strokeWidth=3)
            .encode(
                x=x_labeled,
                y=alt.Y(f"{y_col_p}:Q", title="Active plasmid replisomes"),
            )
            + vline
        ).properties(
            title=f"Plasmid active replisomes around {event_label.lower()} (highlighted seed: {highlighted_seed})",
            width=600,
            height=400,
        ).save(os.path.join(outdir, f"{outfile_prefix}_plasmid.html"))

        # --- Chromosome panel ---
        df_c_win = df_chr[
            (df_chr["time"] >= t_start) & (df_chr["time"] <= t_end)
        ].copy()
        df_c_win["Time (min)"] = df_c_win["time"] / 60

        x_plain = alt.X(
            "Time (min):Q",
            title="Time (min)",
            scale=alt.Scale(domain=(m_start, m_end), zero=False),
            axis=alt.Axis(values=ticks),
        )
        (
            alt.Chart(df_c_win[df_c_win["lineage_seed"] != highlighted_seed])
            .mark_line(opacity=0.2, strokeWidth=1)
            .encode(
                x=x_plain,
                y=alt.Y(f"{y_col_c}:Q", title="Active chromosome replisomes"),
                color=alt.Color("lineage_seed:N", legend=None),
            )
            + alt.Chart(df_c_win[df_c_win["lineage_seed"] == highlighted_seed])
            .mark_line(color="black", strokeWidth=3)
            .encode(
                x=x_plain,
                y=alt.Y(f"{y_col_c}:Q", title="Active chromosome replisomes"),
            )
            + vline
        ).properties(
            title=f"Chromosome active replisomes around {event_label.lower()} (highlighted seed: {highlighted_seed})",
            width=600,
            height=400,
        ).save(os.path.join(outdir, f"{outfile_prefix}_chromosome.html"))

    _make_tight_window_charts(
        chromosome_completion_s,
        "Chromosome completion",
        "tight_window_completion",
    )
    _make_tight_window_charts(
        chromosome_initiation_s,
        "Chromosome initiation",
        "tight_window_initiation",
    )
