"""
Variant: tf_tu_exemption

Corrects cases where a TF fold change has been applied to transcription units
whose promoter TSS lies upstream of the TF's operator/binding site, meaning
the TF cannot actually block those promoters.

For each exempt (TF, TU) pair, the variant:
  1. Absorbs the delta_prob contribution back into basal_prob, preserving the
     net transcription probability at runtime.
  2. Removes the delta_prob entry, making the TU constitutive with respect to
     that TF.

This is motivated by the dnaG/rpsU case: the rpsU-dnaG-rpoD operon has three
transcription units (TU00352, TU00434, TU00435). The lexA operator sits at
genomic position 3,210,729–3,210,748. Only TU00435 (TSS 3,210,735) is actually
blocked by lexA. TU00352 (TSS 3,210,646) and TU00434 (TSS 3,210,712) have TSSs
upstream of the operator and should be constitutively expressed.

Default params apply this fix for lexA on TU00352 and TU00434. The variant is
general and can be configured for any (TF active ID, TU ID) pairs.
"""

from typing import Any, TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from reconstruction.ecoli.simulation_data import SimulationDataEcoli

# Default exemptions: TF active ID -> list of TU IDs whose promoters are
# upstream of the TF operator and should NOT carry the TF's fold change.
DEFAULT_EXEMPTIONS = {
    # lexA (active form PC00010): TU00352 and TU00434 have TSSs upstream of
    # the lexA operator (3,210,729–3,210,748). Only TU00435 is actually blocked.
    "PC00010": ["TU00352", "TU00434"],
}


def apply_variant(
    sim_data: "SimulationDataEcoli", params: dict[str, Any]
) -> "SimulationDataEcoli":
    """
    Remove TF delta_prob contributions from TUs whose promoters are upstream
    of the TF operator, absorbing the delta back into basal_prob.

    Args:
        sim_data: Simulation data (post-parca, basal_prob and delta_prob set).
        params: Parameter dictionary of the following format::

            {
                # Mapping of TF active ID to list of TU IDs to exempt.
                # If omitted, DEFAULT_EXEMPTIONS is used.
                "exemptions": {
                    "PC00010": ["TU00352", "TU00434"],
                    ...
                }
            }

    Returns:
        Simulation data with the following attributes modified::

            sim_data.process.transcription_regulation.basal_prob
            sim_data.process.transcription_regulation.delta_prob
    """
    exemptions = params.get("exemptions", DEFAULT_EXEMPTIONS)

    transcription = sim_data.process.transcription
    transcription_regulation = sim_data.process.transcription_regulation

    # Build index lookups
    all_tu_ids = list(transcription.rna_data["id"])  # e.g. "TU00352[c]"
    all_tf_ids = transcription_regulation.tf_ids  # e.g. "PC00010"

    # tu_id in flat files has no compartment suffix; rna_data IDs have "[c]"
    tu_id_to_idx = {rna_id[:-3]: i for i, rna_id in enumerate(all_tu_ids)}
    tf_id_to_idx = {tf_id: j for j, tf_id in enumerate(all_tf_ids)}

    basal_prob = transcription_regulation.basal_prob
    deltaI = transcription_regulation.delta_prob["deltaI"].copy()
    deltaJ = transcription_regulation.delta_prob["deltaJ"].copy()
    deltaV = transcription_regulation.delta_prob["deltaV"].copy()

    # Mask of entries to KEEP (all True initially)
    keep = np.ones(len(deltaI), dtype=bool)

    for tf_active_id, tu_ids in exemptions.items():
        if tf_active_id not in tf_id_to_idx:
            print(
                f"tf_tu_exemption: TF '{tf_active_id}' not found in tf_ids, skipping."
            )
            continue

        tf_j = tf_id_to_idx[tf_active_id]

        for tu_id in tu_ids:
            if tu_id not in tu_id_to_idx:
                print(f"tf_tu_exemption: TU '{tu_id}' not found in rna_data, skipping.")
                continue

            tu_i = tu_id_to_idx[tu_id]

            # Find the entry for this (TU, TF) pair in delta_prob
            entry_mask = (deltaI == tu_i) & (deltaJ == tf_j)
            n_entries = entry_mask.sum()

            if n_entries == 0:
                print(
                    f"tf_tu_exemption: No delta_prob entry found for "
                    f"TU '{tu_id}' x TF '{tf_active_id}', skipping."
                )
                continue

            # Absorb delta back into basal_prob to preserve net expression.
            # At runtime: effective = basal_prob[tu] + delta_prob[tu, tf] * p_bound
            # For 0CS TFs p_bound ≈ 1, so effective = basal_prob + delta_v.
            # After: effective = (basal_prob + delta_v) + 0 * p_bound = same.
            delta_v = deltaV[entry_mask].sum()
            basal_prob[tu_i] += delta_v
            basal_prob[tu_i] = max(0.0, basal_prob[tu_i])

            # Mark this entry for removal
            keep &= ~entry_mask

            print(
                f"tf_tu_exemption: Removed lexA delta ({delta_v:.4f}) from "
                f"TU '{tu_id}', absorbed into basal_prob (new value: "
                f"{basal_prob[tu_i]:.4f})."
            )

    # Rebuild delta_prob arrays with exempt entries removed
    transcription_regulation.basal_prob = basal_prob
    transcription_regulation.delta_prob = {
        "deltaI": deltaI[keep],
        "deltaJ": deltaJ[keep],
        "deltaV": deltaV[keep],
        "shape": transcription_regulation.delta_prob["shape"],
    }

    return sim_data
