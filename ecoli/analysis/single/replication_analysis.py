import os
from typing import Any

from duckdb import DuckDBPyConnection
import altair as alt


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
    query = f"""
    SELECT time, listeners__unique_molecule_counts__active_replisome
    FROM ({history_sql})
    ORDER BY time ASC
    """

    output_df = conn.sql(query).df()
    # Convert time from seconds to minutes
    output_df["Time (min)"] = output_df["time"] / 60

    # Create Altair line chart
    chart = (
        alt.Chart(output_df)
        .mark_line()
        .encode(
            x=alt.X("Time (min):Q", title="Time (min)"),
            y=alt.Y(
                "listeners__unique_molecule_counts__active_replisome:Q",
                title="active replisome counts",
            ),
        )
        .properties(
            title="Active replisome counts on chromosome",
            width=600,
            height=400,
        )
    )

    html_path = os.path.join(outdir, "active_replisome_counts_plasmid.html")
    chart.save(html_path)

    query2 = f"""
        SELECT time, listeners__mass__cell_mass
        FROM ({history_sql})
        ORDER BY time ASC
        """

    output_df2 = conn.sql(query2).df()
    # Convert time from seconds to minutes
    output_df2["Time (min)"] = output_df2["time"] / 60

    # Create Altair line chart
    chart_cellmass = (
        alt.Chart(output_df2)
        .mark_line()
        .encode(
            x=alt.X("Time (min):Q", title="Time (min)"),
            y=alt.Y(
                "listeners__mass__cell_mass:Q",
                title="cell mass (fg)",
            ),
        )
        .properties(
            title="cell mass (with plasmid)",
            width=600,
            height=400,
        )
    )

    html_path2 = os.path.join(outdir, "cell mass (plasmid).html")
    chart_cellmass.save(html_path2)

    query3 = f"""
            SELECT time, listeners__mass__dna_mass
            FROM ({history_sql})
            ORDER BY time ASC
            """

    output_df3 = conn.sql(query3).df()
    # Convert time from seconds to minutes
    output_df3["Time (min)"] = output_df3["time"] / 60

    # Create Altair line chart
    chart_dnamass = (
        alt.Chart(output_df3)
        .mark_line()
        .encode(
            x=alt.X("Time (min):Q", title="Time (min)"),
            y=alt.Y(
                "listeners__mass__dna_mass:Q",
                title="dna mass (fg)",
            ),
        )
        .properties(
            title="dna mass (with plasmid)",
            width=600,
            height=400,
        )
    )

    html_path3 = os.path.join(outdir, "dna mass (plasmid).html")
    chart_dnamass.save(html_path3)

    query4 = f"""
                SELECT time, listeners__mass__growth
                FROM ({history_sql})
                ORDER BY time ASC
                """

    output_df4 = conn.sql(query4).df()
    # Convert time from seconds to minutes
    output_df4["Time (min)"] = output_df4["time"] / 60

    # Create Altair line chart
    chart_growth = (
        alt.Chart(output_df4)
        .mark_line()
        .encode(
            x=alt.X("Time (min):Q", title="Time (min)"),
            y=alt.Y(
                "listeners__mass__growth:Q",
                title="growth",
            ),
        )
        .properties(
            title="growth (with plasmid)",
            width=600,
            height=400,
        )
    )

    html_path4 = os.path.join(outdir, "growth (plasmid).html")
    chart_growth.save(html_path4)

    query5 = f"""
                    SELECT time, listeners__mass__instantaneous_growth_rate
                    FROM ({history_sql})
                    ORDER BY time ASC
                    """

    output_df5 = conn.sql(query5).df()
    # Convert time from seconds to minutes
    output_df5["Time (min)"] = output_df5["time"] / 60

    # Create Altair line chart
    chart_instantaneous_growth = (
        alt.Chart(output_df5)
        .mark_line()
        .encode(
            x=alt.X("Time (min):Q", title="Time (min)"),
            y=alt.Y(
                "listeners__mass__instantaneous_growth_rate:Q",
                title="instantaneous growth",
            ),
        )
        .properties(
            title="instantaneous growth (with plasmid)",
            width=600,
            height=400,
        )
    )

    html_path5 = os.path.join(outdir, "instantaneous growth (plasmid).html")
    chart_instantaneous_growth.save(html_path5)

    # query6 = f"""
    #                     SELECT time, listeners__fba_results__reaction_fluxes
    #                     FROM ({history_sql})
    #                     ORDER BY time ASC
    #                     """
    #
    # output_df6 = conn.sql(query6).df()
    # # Convert time from seconds to minutes
    # time_minutes = output_df6["time"] / 60
    # flux_matrix = np.stack(
    #     output_df6["listeners__fba_results__reaction_fluxes"].values
    # ).astype(float)
    # flux_df = pd.DataFrame(flux_matrix)
    # flux_df.insert(0, "Time (min)", time_minutes)
    # flux_df.to_csv(os.path.join(outdir, "flux_matrix.csv"), index=False)
    #
    # #'adenosine deoxyribonucleotides de novo biosynthesis II'
    # adenosine_cols = [flux_df.columns[i] for i in [602, 814, 6771, 6779]]
    # adenosine_df = flux_df[["Time (min)"] + adenosine_cols].copy()
    #
    # # Rename for clarity
    # adenosine_df.columns = [
    #     "Time (min)",
    #     "ADPREDUCT-RXN",
    #     "DADPKIN-RXN",
    #     "RXN0-745",
    #     "RXN0-747",
    # ]
    #
    # # Melt for plotting
    # melted_adenosine_df = pd.melt(
    #     adenosine_df, id_vars=["Time (min)"], var_name="Reaction", value_name="Flux"
    # )
    #
    # # Extract flux values (exclude Time column)
    # flux_values = adenosine_df.iloc[:, 1:].values  # shape: (timepoints, reactions)
    #
    # # Find non-zero values
    # nonzero_fluxes = flux_values[flux_values != 0]
    #
    # # Determine y-axis domain
    # if nonzero_fluxes.size > 0:
    #     # Use max non-zero flux, padded slightly
    #     max_flux = nonzero_fluxes.max() * 1.05
    #     y_scale = alt.Scale(domain=[0, max_flux])
    # else:
    #     # All fluxes are zero → let Altair auto-scale (will show flat line at 0)
    #     y_scale = alt.Undefined
    #
    # # Plot
    # chart_adenosine = (
    #     alt.Chart(melted_adenosine_df)
    #     .mark_line()
    #     .encode(
    #         x=alt.X("Time (min):Q", title="Time (min)"),
    #         y=alt.Y("Flux:Q", title="Flux", scale=y_scale),
    #         color=alt.Color("Reaction:N", scale=alt.Scale(range=COLORS)),
    #     )
    #     .properties(
    #         title="Adenosine Deoxyribonucleotide Biosynthesis reaction fluxes (With plasmid) "
    #     )
    # )
    #
    # # Save or display
    # chart_adenosine.save(os.path.join(outdir, "adenosine_counts.html"))
    #
    # #'guanosine deoxyribonucleotides de novo biosynthesis II'
    # guanosine_cols = [flux_df.columns[i] for i in [862, 1034, 6778, 6780]]
    # guanosine_df = flux_df[["Time (min)"] + guanosine_cols].copy()
    #
    # # Rename for clarity
    # guanosine_df.columns = [
    #     "Time (min)",
    #     "DGDPKIN-RXN",
    #     "GDPREDUCT-RXN",
    #     "RXN0-746",
    #     "RXN0-748",
    # ]
    #
    # # Melt for plotting
    # melted_guanosine_df = pd.melt(
    #     guanosine_df, id_vars=["Time (min)"], var_name="Reaction", value_name="Flux"
    # )
    #
    # # Extract flux values (exclude Time column)
    # flux_values = guanosine_df.iloc[:, 1:].values  # shape: (timepoints, reactions)
    #
    # # Find non-zero values
    # nonzero_fluxes = flux_values[flux_values != 0]
    #
    # # Determine y-axis domain
    # if nonzero_fluxes.size > 0:
    #     # Use max non-zero flux, padded slightly
    #     max_flux = nonzero_fluxes.max() * 1.05
    #     y_scale = alt.Scale(domain=[0, max_flux])
    # else:
    #     # All fluxes are zero → let Altair auto-scale (will show flat line at 0)
    #     y_scale = alt.Undefined
    #
    # # Plot
    # chart_guanosine = (
    #     alt.Chart(melted_guanosine_df)
    #     .mark_line()
    #     .encode(
    #         x=alt.X("Time (min):Q", title="Time (min)"),
    #         y=alt.Y("Flux:Q", title="Flux", scale=y_scale),
    #         color=alt.Color("Reaction:N", scale=alt.Scale(range=COLORS)),
    #     )
    #     .properties(
    #         title="Guanosine Deoxyribonucleotide Biosynthesis reaction fluxes (With plasmid) "
    #     )
    # )
    #
    # # Save or display
    # chart_guanosine.save(os.path.join(outdir, "guanosine_counts.html"))
    #
    # # 'pyrimidine deoxyribonucleotides de novo biosynthesis II'
    # pyrimidine2_cols = [flux_df.columns[i] for i in [845, 918, 919, 922, 6054, 6062]]
    # pyrimidine2_df = flux_df[["Time (min)"] + pyrimidine2_cols].copy()
    #
    # # Rename for clarity
    # pyrimidine2_df.columns = [
    #     "Time (min)",
    #     "DCTP-DEAM-RXN",
    #     "DTDPKIN-RXN",
    #     "DTMPKI-RXN",
    #     "DUTP-PYROP-RXN",
    #     "RXN0-723",
    #     "RXN0-724",
    # ]
    #
    # # Melt for plotting
    # melted_pyrimidine2_df = pd.melt(
    #     pyrimidine2_df, id_vars=["Time (min)"], var_name="Reaction", value_name="Flux"
    # )
    #
    # # Extract flux values (exclude Time column)
    # flux_values = pyrimidine2_df.iloc[:, 1:].values  # shape: (timepoints, reactions)
    #
    # # Find non-zero values
    # nonzero_fluxes = flux_values[flux_values != 0]
    #
    # # Determine y-axis domain
    # if nonzero_fluxes.size > 0:
    #     # Use max non-zero flux, padded slightly
    #     max_flux = nonzero_fluxes.max() * 1.05
    #     y_scale = alt.Scale(domain=[0, max_flux])
    # else:
    #     # All fluxes are zero → let Altair auto-scale (will show flat line at 0)
    #     y_scale = alt.Undefined
    #
    # # Plot
    # chart_pyrimidine2 = (
    #     alt.Chart(melted_pyrimidine2_df)
    #     .mark_line()
    #     .encode(
    #         x=alt.X("Time (min):Q", title="Time (min)"),
    #         y=alt.Y("Flux:Q", title="Flux", scale=y_scale),
    #         color=alt.Color("Reaction:N", scale=alt.Scale(range=COLORS)),
    #     )
    #     .properties(
    #         title="Pyrimidine Deoxyribonucleotide Biosynthesis II reaction fluxes (With plasmid) "
    #     )
    # )
    #
    # # Save or display
    # chart_pyrimidine2.save(os.path.join(outdir, "pyrimidine2_counts.html"))
    #
    # #'pyrimidine deoxyribonucleotides de novo biosynthesis I'
    # pyrimidine1_cols = [
    #     flux_df.columns[i] for i in [762, 844, 918, 919, 920, 922, 1755, 9192]
    # ]
    # pyrimidine1_df = flux_df[["Time (min)"] + pyrimidine1_cols].copy()
    #
    # # Rename for clarity
    # pyrimidine1_df.columns = [
    #     "Time (min)",
    #     "CDPREDUCT-RXN",
    #     "DCDPKIN-RXN",
    #     "DTDPKIN-RXN",
    #     "DTMPKI-RXN",
    #     "DUDPKIN-RXN",
    #     "DUTP-PYROP-RXN",
    #     "RXN-12195",
    #     "UDPREDUCT-RXN",
    # ]
    #
    # # Melt for plotting
    # melted_pyrimidine1_df = pd.melt(
    #     pyrimidine1_df, id_vars=["Time (min)"], var_name="Reaction", value_name="Flux"
    # )
    #
    # # Extract flux values (exclude Time column)
    # flux_values = pyrimidine1_df.iloc[:, 1:].values  # shape: (timepoints, reactions)
    #
    # # Find non-zero values
    # nonzero_fluxes = flux_values[flux_values != 0]
    #
    # # Determine y-axis domain
    # if nonzero_fluxes.size > 0:
    #     # Use max non-zero flux, padded slightly
    #     max_flux = nonzero_fluxes.max() * 1.05
    #     y_scale = alt.Scale(domain=[0, max_flux])
    # else:
    #     # All fluxes are zero → let Altair auto-scale (will show flat line at 0)
    #     y_scale = alt.Undefined
    #
    # # Plot
    # chart_pyrimidine1 = (
    #     alt.Chart(melted_pyrimidine1_df)
    #     .mark_line()
    #     .encode(
    #         x=alt.X("Time (min):Q", title="Time (min)"),
    #         y=alt.Y("Flux:Q", title="Flux", scale=y_scale),
    #         color=alt.Color("Reaction:N", scale=alt.Scale(range=COLORS)),
    #     )
    #     .properties(
    #         title="Pyrimidine Deoxyribonucleotide Biosynthesis I reaction fluxes (With plasmid) "
    #     )
    # )
    #
    # # Save or display
    # chart_pyrimidine1.save(os.path.join(outdir, "pyrimidine1_counts.html"))
    #
    # # Replisome subunits
    #
    # query7 = f"""
    #         SELECT time, bulk
    #         FROM ({history_sql})
    #         ORDER BY time ASC
    #         """
    # output_df_replisome = conn.sql(query7).df()
    # time_minutes = output_df_replisome["time"].to_numpy() / 60
    # bulk_matrix = np.stack(output_df_replisome["bulk"].values).astype(int)
    # bulk_df = pd.DataFrame(bulk_matrix)
    # bulk_df.insert(0, "Time (min)", time_minutes)
    #
    # replisome_subunit_cols = [
    #     bulk_df.columns[i] for i in [5377, 6568, 7233, 7189, 5305, 5383]
    # ]
    # replisome_subunit_df = bulk_df[["Time (min)"] + replisome_subunit_cols].copy()
    #
    # # Rename for clarity
    # replisome_subunit_df.columns = [
    #     "Time (min)",
    #     "CPLX0-3621[c]",
    #     "EG10239-MONOMER[c]",
    #     "EG11500-MONOMER[c]",
    #     "EG11412-MONOMER[c]",
    #     "CPLX0-2361[c]",
    #     "CPLX0-3761[c]",
    # ]
    #
    # # Melt for plotting
    # melted_replisome_subunit_df = pd.melt(
    #     replisome_subunit_df,
    #     id_vars=["Time (min)"],
    #     var_name="Subunit",
    #     value_name="Counts",
    # )
    #
    # chart_replisome_subunits = (
    #     alt.Chart(melted_replisome_subunit_df)
    #     .mark_line()
    #     .encode(
    #         x=alt.X("Time (min):Q", title="Time (min)"),
    #         y=alt.Y("Counts:Q", title="Counts"),
    #         color=alt.Color("Subunit:N", scale=alt.Scale(range=COLORS)),
    #     )
    #     .properties(title="Replisome subunit counts over time (with plasmids) ")
    # )
    #
    # # Save or display
    # chart_replisome_subunits.save(os.path.join(outdir, "replisome_subunit_counts.html"))
    #
    # # dntps
    # dntp_cols = [bulk_df.columns[i] for i in [6211, 6222, 6299, 12449]]
    # dntp_df = bulk_df[["Time (min)"] + dntp_cols].copy()
    #
    # # Rename for clarity
    # dntp_df.columns = ["Time (min)", "DATP (C)", "DCTP (C)", "DGTP (C)", "TTP (C)"]
    #
    # # Melt for plotting
    # melted_dntp_df = pd.melt(
    #     dntp_df, id_vars=["Time (min)"], var_name="Molecule", value_name="Counts"
    # )
    #
    # chart_dntp = (
    #     alt.Chart(melted_dntp_df)
    #     .mark_line()
    #     .encode(
    #         x=alt.X("Time (min):Q", title="Time (min)"),
    #         y=alt.Y("Counts:Q", title="Counts"),
    #         color=alt.Color("Molecule:N", scale=alt.Scale(range=COLORS)),
    #     )
    #     .properties(title="dntp counts over time (with plasmids)")
    # )
    #
    # # Save or display
    # chart_dntp.save(os.path.join(outdir, "dntp_counts.html"))
