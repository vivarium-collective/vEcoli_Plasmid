"""
Multigeneration analysis of plasmid copy number dynamics.

Plots:
1. Plasmid copy number over time across all generations
2. Copy number per cell at end of each generation (convergence to steady state)
3. RNA I and RNA II levels over time (copy number control mechanism)
"""

import os
from typing import Any

import altair as alt
import polars as pl
from duckdb import DuckDBPyConnection

from ecoli.library.parquet_emitter import read_stacked_columns


def plot(
    params: dict[str, Any],
    conn: DuckDBPyConnection,
    history_sql: str,
    config_sql: str,
    success_sql: str,
    sim_data_dict: dict[str, dict[int, str]],
    validation_data_paths: list[str],
    outdir: str,
    variant_metadata: dict[str, dict[int, Any]],
    variant_names: dict[str, str],
):
    # ------------------------------------------------------------------
    # Load plasmid count, RNA control state, and active replisome count
    # across all generations. `generation` column lets us track convergence.
    # ------------------------------------------------------------------
    data_columns = [
        "generation",
        "time",
        "listeners__unique_molecule_counts__full_plasmid AS plasmid_count",
        "listeners__unique_molecule_counts__plasmid_active_replisome AS active_replisomes",
        "process_state__plasmid_rna_control__rna_I AS rna_I",
        "process_state__plasmid_rna_control__rna_II AS rna_II",
    ]

    df = pl.DataFrame(
        read_stacked_columns(history_sql, data_columns, conn=conn)
    ).with_columns(
        # time in minutes within each generation (time resets at division)
        time_min=pl.col("time") / 60,
        generation=pl.col("generation").cast(pl.Int32),
    )

    n_generations = df["generation"].n_unique()

    # ------------------------------------------------------------------
    # Plot 1: Plasmid count vs time, one line per generation
    # ------------------------------------------------------------------
    chart_count = (
        alt.Chart(df.to_pandas())
        .mark_line()
        .encode(
            x=alt.X("time_min:Q", title="Time within generation (min)"),
            y=alt.Y("plasmid_count:Q", title="Plasmid copy number"),
            color=alt.Color(
                "generation:O",
                title="Generation",
                scale=alt.Scale(scheme="viridis"),
            ),
            tooltip=["generation:O", "time_min:Q", "plasmid_count:Q"],
        )
        .properties(
            title=f"Plasmid copy number across {n_generations} generations",
            width=700,
            height=400,
        )
    )
    chart_count.save(os.path.join(outdir, "plasmid_copy_number_by_generation.html"))

    # ------------------------------------------------------------------
    # Plot 2: End-of-generation copy number — shows convergence to ~23
    # ------------------------------------------------------------------
    end_of_gen = (
        df.group_by("generation")
        .agg(pl.col("plasmid_count").last().alias("final_count"))
        .sort("generation")
    )

    chart_convergence = (
        alt.Chart(end_of_gen.to_pandas())
        .mark_line(point=True)
        .encode(
            x=alt.X("generation:O", title="Generation"),
            y=alt.Y(
                "final_count:Q",
                title="Plasmid count at end of generation",
                scale=alt.Scale(zero=True),
            ),
            tooltip=["generation:O", "final_count:Q"],
        )
        .properties(
            title="Plasmid copy number convergence to steady state",
            width=600,
            height=350,
        )
    )

    # Add reference line at 23 (Ataai & Shuler 1986 prediction)
    ref_line = (
        alt.Chart(pl.DataFrame({"y": [23]}).to_pandas())
        .mark_rule(color="red", strokeDash=[6, 3])
        .encode(y="y:Q")
    )

    (chart_convergence + ref_line).save(
        os.path.join(outdir, "plasmid_copy_number_convergence.html")
    )

    # ------------------------------------------------------------------
    # Plot 3: RNA I and RNA II levels across generations
    # ------------------------------------------------------------------
    rna_df = (
        df.select(["generation", "time_min", "rna_I", "rna_II"])
        .unpivot(
            on=["rna_I", "rna_II"],
            index=["generation", "time_min"],
            variable_name="species",
            value_name="count",
        )
        .with_columns(
            pl.col("species").replace(
                {"rna_I": "RNA I (inhibitor)", "rna_II": "RNA II (primer)"}
            )
        )
    )

    chart_rna = (
        alt.Chart(rna_df.to_pandas())
        .mark_line()
        .encode(
            x=alt.X("time_min:Q", title="Time within generation (min)"),
            y=alt.Y("count:Q", title="Molecule count"),
            color=alt.Color("species:N", title="Species"),
            facet=alt.Facet("generation:O", title="Generation", columns=5),
        )
        .properties(
            title="RNA I / RNA II dynamics (Ataai-Shuler copy number control)",
            width=180,
            height=120,
        )
    )
    chart_rna.save(os.path.join(outdir, "rna_control_dynamics.html"))

    # ------------------------------------------------------------------
    # Plot 4: Active replisomes per generation
    # ------------------------------------------------------------------
    chart_replisome = (
        alt.Chart(df.to_pandas())
        .mark_line()
        .encode(
            x=alt.X("time_min:Q", title="Time within generation (min)"),
            y=alt.Y("active_replisomes:Q", title="Active plasmid replisomes"),
            color=alt.Color(
                "generation:O",
                title="Generation",
                scale=alt.Scale(scheme="viridis"),
            ),
        )
        .properties(
            title="Active plasmid replisomes across generations",
            width=700,
            height=350,
        )
    )
    chart_replisome.save(
        os.path.join(outdir, "plasmid_active_replisomes_multigen.html")
    )
