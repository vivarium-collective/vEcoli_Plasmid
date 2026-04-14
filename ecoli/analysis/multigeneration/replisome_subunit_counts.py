"""
Multigeneration analysis of replisome subunit counts.

Tracks the bulk counts of each chromosome replisome subunit across generations
to diagnose stochastic depletion that can cause TimeLimitError (cell never
divides because it cannot assemble a replisome).

Monomer subunits (2 required per oriC at initiation):
  - EG10239-MONOMER[c]  DnaG  (DNA primase)
  - EG11500-MONOMER[c]  DnaN  (beta clamp)
  - EG11412-MONOMER[c]  SSB   (single-strand binding protein)
  - CPLX0-3621[c]       Pol III clamp loader (gamma/delta complex)

Trimer subunits (6 required per oriC at initiation):
  - CPLX0-2361[c]       DnaB helicase complex
  - CPLX0-3761[c]       Tau complex

Plots:
1. Count of each monomer subunit at end of each generation
2. Count of each trimer subunit at end of each generation
3. Minimum count (across all subunit types) per generation — the actual
   bottleneck for replisome assembly
4. Full time-series of the most-depleted subunit (DnaG) across all generations
"""

import os
import pickle
from typing import Any

import altair as alt
import polars as pl
from duckdb import DuckDBPyConnection

from ecoli.library.parquet_emitter import open_arbitrary_sim_data, read_stacked_columns

# Bulk array indices for each replisome subunit (from simData).
# These are 0-based indices into the `bulk` array stored in each parquet row.
# Determined from sim_data.internal_state.bulk_molecules.bulk_data['id'].
MONOMER_SUBUNITS = {
    "DnaG (primase)": 6613,  # EG10239-MONOMER[c]
    "DnaN (beta clamp)": 7277,  # EG11500-MONOMER[c]
    "SSB": 7233,  # EG11412-MONOMER[c]
    "Clamp loader (gamma/d)": 5418,  # CPLX0-3621[c]
}

TRIMER_SUBUNITS = {
    "DnaB helicase": 5346,  # CPLX0-2361[c]
    "Tau complex": 5424,  # CPLX0-3761[c]
}

# Required counts per oriC at initiation
MONOMERS_PER_ORIC = 2
TRIMERS_PER_ORIC = 6


def _bulk_index_expr(idx: int) -> str:
    """DuckDB expression to extract one element from the bulk array column."""
    # DuckDB arrays are 1-indexed
    return f"bulk[{idx + 1}]"


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
    # Look up dnaG cistron indices from simData
    # ------------------------------------------------------------------
    with open_arbitrary_sim_data(sim_data_dict) as f:
        sim_data = pickle.load(f)

    transcription = sim_data.process.transcription
    all_cistron_ids = list(transcription.cistron_data["id"])
    cistron_is_mrna = transcription.cistron_data["is_mRNA"]

    dnag_cistron_idx = next(
        i for i, c in enumerate(all_cistron_ids) if c is not None and "EG10239" in c
    )

    mrna_cistron_ids = [
        c for c, is_mrna in zip(all_cistron_ids, cistron_is_mrna) if is_mrna
    ]
    dnag_mrna_idx = next(
        i for i, c in enumerate(mrna_cistron_ids) if c is not None and "EG10239" in c
    )

    # ------------------------------------------------------------------
    # Load bulk counts for all replisome subunits + generation label
    # ------------------------------------------------------------------
    all_subunits = {**MONOMER_SUBUNITS, **TRIMER_SUBUNITS}

    # Build SELECT expressions for each subunit
    subunit_selects = [
        f"{_bulk_index_expr(idx)} AS {name.replace(' ', '_').replace('/', '_').replace('(', '').replace(')', '')}"
        for name, idx in all_subunits.items()
    ]
    col_map = {
        name.replace(" ", "_").replace("/", "_").replace("(", "").replace(")", ""): name
        for name in all_subunits
    }

    data_columns = (
        ["generation", "time"]
        + subunit_selects
        + [
            "listeners__unique_molecule_counts__full_chromosome AS full_chromosome",
            "listeners__unique_molecule_counts__active_replisome AS active_replisome",
            f"listeners__rnap_data__rna_init_event_per_cistron[{dnag_cistron_idx + 1}]"
            " AS dnag_init_events",
            f"listeners__rna_counts__mRNA_cistron_counts[{dnag_mrna_idx + 1}]"
            " AS dnag_mrna_count",
        ]
    )

    df = pl.DataFrame(
        read_stacked_columns(history_sql, data_columns, conn=conn)
    ).with_columns(
        time_min=pl.col("time") / 60,
        generation=pl.col("generation").cast(pl.Int32),
    )

    n_generations = df["generation"].n_unique()

    # ------------------------------------------------------------------
    # End-of-generation counts: last recorded value for each generation
    # ------------------------------------------------------------------
    agg_exprs = [pl.col(col).last().alias(col) for col in col_map.keys()]
    end_of_gen = df.group_by("generation").agg(agg_exprs).sort("generation")

    # Reshape to long format for faceted plot
    end_long = end_of_gen.unpivot(
        on=list(col_map.keys()),
        index=["generation"],
        variable_name="col_name",
        value_name="count",
    ).with_columns(
        pl.col("col_name").replace(col_map).alias("subunit"),
        pl.col("col_name")
        .map_elements(
            lambda x: "Monomer (need ≥2/oriC)"
            if col_map[x] in MONOMER_SUBUNITS
            else "Trimer (need ≥6/oriC)",
            return_dtype=pl.Utf8,
        )
        .alias("subunit_type"),
    )

    # ------------------------------------------------------------------
    # Plot 1: End-of-generation count per subunit (monomer subunits)
    # ------------------------------------------------------------------
    monomer_end = end_long.filter(pl.col("subunit_type") == "Monomer (need ≥2/oriC)")

    chart_monomers = (
        alt.Chart(monomer_end.to_pandas())
        .mark_line(point=True)
        .encode(
            x=alt.X("generation:O", title="Generation"),
            y=alt.Y(
                "count:Q",
                title="Count at end of generation",
                scale=alt.Scale(zero=True),
            ),
            color=alt.Color("subunit:N", title="Subunit"),
            tooltip=["generation:O", "subunit:N", "count:Q"],
        )
        .properties(
            title=f"Replisome monomer subunit counts across {n_generations} generations",
            width=650,
            height=350,
        )
    )

    # Reference line: minimum needed per oriC
    ref_monomer = (
        alt.Chart(pl.DataFrame({"y": [MONOMERS_PER_ORIC]}).to_pandas())
        .mark_rule(color="red", strokeDash=[6, 3])
        .encode(y="y:Q")
    )

    (chart_monomers + ref_monomer).save(
        os.path.join(outdir, "replisome_monomer_subunit_counts.html")
    )

    # ------------------------------------------------------------------
    # Plot 2: End-of-generation count per subunit (trimer subunits)
    # ------------------------------------------------------------------
    trimer_end = end_long.filter(pl.col("subunit_type") == "Trimer (need ≥6/oriC)")

    chart_trimers = (
        alt.Chart(trimer_end.to_pandas())
        .mark_line(point=True)
        .encode(
            x=alt.X("generation:O", title="Generation"),
            y=alt.Y(
                "count:Q",
                title="Count at end of generation",
                scale=alt.Scale(zero=True),
            ),
            color=alt.Color("subunit:N", title="Subunit"),
            tooltip=["generation:O", "subunit:N", "count:Q"],
        )
        .properties(
            title=f"Replisome trimer subunit counts across {n_generations} generations",
            width=650,
            height=350,
        )
    )

    ref_trimer = (
        alt.Chart(pl.DataFrame({"y": [TRIMERS_PER_ORIC]}).to_pandas())
        .mark_rule(color="red", strokeDash=[6, 3])
        .encode(y="y:Q")
    )

    (chart_trimers + ref_trimer).save(
        os.path.join(outdir, "replisome_trimer_subunit_counts.html")
    )

    # ------------------------------------------------------------------
    # Plot 3: Effective bottleneck — minimum normalised availability
    # (count / required) per generation, for each subunit
    # ------------------------------------------------------------------
    # Compute "availability ratio" = count / required at end of each gen
    bottleneck_rows = []
    for row in end_of_gen.iter_rows(named=True):
        gen = row["generation"]
        for col, name in col_map.items():
            required = (
                MONOMERS_PER_ORIC if name in MONOMER_SUBUNITS else TRIMERS_PER_ORIC
            )
            bottleneck_rows.append(
                {
                    "generation": gen,
                    "subunit": name,
                    "count": row[col],
                    "required": required,
                    "ratio": row[col] / required,
                }
            )

    bottleneck_df = pl.DataFrame(bottleneck_rows)

    chart_bottleneck = (
        alt.Chart(bottleneck_df.to_pandas())
        .mark_line(point=True)
        .encode(
            x=alt.X("generation:O", title="Generation"),
            y=alt.Y(
                "ratio:Q",
                title="Count / required per oriC",
                scale=alt.Scale(zero=True),
            ),
            color=alt.Color("subunit:N", title="Subunit"),
            tooltip=["generation:O", "subunit:N", "count:Q", "ratio:Q"],
        )
        .properties(
            title="Replisome subunit availability (count ÷ required) — values <1 block initiation",
            width=700,
            height=400,
        )
    )

    ref_one = (
        alt.Chart(pl.DataFrame({"y": [1.0]}).to_pandas())
        .mark_rule(color="red", strokeDash=[6, 3])
        .encode(y="y:Q")
    )

    (chart_bottleneck + ref_one).save(
        os.path.join(outdir, "replisome_subunit_bottleneck.html")
    )

    # ------------------------------------------------------------------
    # Plot 4: Full time-series of DnaG (the most-depleted subunit)
    # ------------------------------------------------------------------
    dnag_col = next(col for col, name in col_map.items() if name == "DnaG (primase)")

    chart_dnag = (
        alt.Chart(df.select(["generation", "time_min", dnag_col]).to_pandas())
        .mark_line()
        .encode(
            x=alt.X("time_min:Q", title="Time within generation (min)"),
            y=alt.Y(f"{dnag_col}:Q", title="DnaG (DNA primase) count"),
            color=alt.Color(
                "generation:O",
                title="Generation",
                scale=alt.Scale(scheme="viridis"),
            ),
            tooltip=["generation:O", "time_min:Q", f"{dnag_col}:Q"],
        )
        .properties(
            title="DnaG (DNA primase) count across generations — stochastic depletion",
            width=700,
            height=400,
        )
    )

    ref_two = (
        alt.Chart(pl.DataFrame({"y": [MONOMERS_PER_ORIC]}).to_pandas())
        .mark_rule(color="red", strokeDash=[6, 3])
        .encode(y=alt.Y("y:Q", title="DnaG count"))
        .properties(title="Minimum needed for initiation (2/oriC)")
    )

    (chart_dnag + ref_two).save(os.path.join(outdir, "dnag_count_multigen.html"))

    # ------------------------------------------------------------------
    # Plot 5: Full chromosome count across generations
    # ------------------------------------------------------------------
    chart_chromosomes = (
        alt.Chart(df.select(["generation", "time_min", "full_chromosome"]).to_pandas())
        .mark_line()
        .encode(
            x=alt.X("time_min:Q", title="Time within generation (min)"),
            y=alt.Y("full_chromosome:Q", title="Full chromosome count"),
            color=alt.Color(
                "generation:O",
                title="Generation",
                scale=alt.Scale(scheme="viridis"),
            ),
            tooltip=["generation:O", "time_min:Q", "full_chromosome:Q"],
        )
        .properties(
            title="Full chromosome count across generations",
            width=700,
            height=350,
        )
    )
    chart_chromosomes.save(os.path.join(outdir, "full_chromosome_count_multigen.html"))

    # ------------------------------------------------------------------
    # Plot 6: Active chromosome replisomes across generations
    # ------------------------------------------------------------------
    chart_chr_replisome = (
        alt.Chart(df.select(["generation", "time_min", "active_replisome"]).to_pandas())
        .mark_line()
        .encode(
            x=alt.X("time_min:Q", title="Time within generation (min)"),
            y=alt.Y("active_replisome:Q", title="Active chromosome replisomes"),
            color=alt.Color(
                "generation:O",
                title="Generation",
                scale=alt.Scale(scheme="viridis"),
            ),
            tooltip=["generation:O", "time_min:Q", "active_replisome:Q"],
        )
        .properties(
            title="Active chromosome replisomes across generations",
            width=700,
            height=350,
        )
    )
    chart_chr_replisome.save(
        os.path.join(outdir, "active_chromosome_replisomes_multigen.html")
    )

    # ------------------------------------------------------------------
    # Plot 7: dnaG transcription initiation events per timestep
    # ------------------------------------------------------------------
    chart_dnag_init = (
        alt.Chart(df.select(["generation", "time_min", "dnag_init_events"]).to_pandas())
        .mark_line()
        .encode(
            x=alt.X("time_min:Q", title="Time within generation (min)"),
            y=alt.Y("dnag_init_events:Q", title="dnaG initiation events per timestep"),
            color=alt.Color(
                "generation:O",
                title="Generation",
                scale=alt.Scale(scheme="viridis"),
            ),
            tooltip=["generation:O", "time_min:Q", "dnag_init_events:Q"],
        )
        .properties(
            title="dnaG transcription initiation events across generations",
            width=700,
            height=350,
        )
    )
    chart_dnag_init.save(os.path.join(outdir, "dnag_initiation_events_multigen.html"))

    # ------------------------------------------------------------------
    # Plot 8: dnaG mRNA count over time
    # ------------------------------------------------------------------
    chart_dnag_mrna = (
        alt.Chart(df.select(["generation", "time_min", "dnag_mrna_count"]).to_pandas())
        .mark_line()
        .encode(
            x=alt.X("time_min:Q", title="Time within generation (min)"),
            y=alt.Y("dnag_mrna_count:Q", title="dnaG mRNA count"),
            color=alt.Color(
                "generation:O",
                title="Generation",
                scale=alt.Scale(scheme="viridis"),
            ),
            tooltip=["generation:O", "time_min:Q", "dnag_mrna_count:Q"],
        )
        .properties(
            title="dnaG mRNA count across generations",
            width=700,
            height=350,
        )
    )
    chart_dnag_mrna.save(os.path.join(outdir, "dnag_mrna_count_multigen.html"))
