from ecoli.library.schema import bulk_name_to_idx
import pandas as pd

# %%
# Loading Simdata to obtain the bulk id labels
from ecoli.library.sim_data import LoadSimData

sim_data_default = "out/plasmid_in_ecoli2/parca/kb/simData.cPickle"
sim_data = LoadSimData(sim_data_default).sim_data

bulk_molecule_ids = sim_data.internal_state.bulk_molecules.bulk_data[
    "id"
].tolist()  # Model common name with compartments
# %%


df1 = pd.DataFrame(bulk_molecule_ids, columns=["bulk_id"])
df1.to_csv("bulk_molecule_ids1.csv", index=False)
# %%
dntp_ids = bulk_name_to_idx(
    ["DATP[c]", "DCTP[c]", "DGTP[c]", "TTP[c]"], bulk_molecule_ids
)
polymerized_dntp_ids = bulk_name_to_idx(
    [
        "polymerized_DATP[c]",
        "polymerized_DCTP[c]",
        "polymerized_DGTP[c]",
        "polymerized_TTP[c]",
    ],
    bulk_molecule_ids,
)
ntp_ids = bulk_name_to_idx(["ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]"], bulk_molecule_ids)
polymerized_ntp_ids = bulk_name_to_idx(
    [
        "polymerized_ATP[c]",
        "polymerized_CTP[c]",
        "polymerized_GTP[c]",
        "polymerized_UTP[c]",
    ],
    bulk_molecule_ids,
)
# %%
RNAP_subunits_ids = bulk_name_to_idx(
    ["EG10893-MONOMER[c]", "RPOC-MONOMER[c]", "RPOB-MONOMER[c]"], bulk_molecule_ids
)

# %%
RNase_ids = bulk_name_to_idx(["EG10861[c]", "EG10860[c]"], bulk_molecule_ids)

# %%
rna_data = sim_data.process.transcription.rna_data

mrna_tu_ids = rna_data["id"][rna_data["is_mRNA"]].tolist()

cistron_data = sim_data.process.transcription.cistron_data

mrna_cistron_ids = cistron_data["id"][cistron_data["is_mRNA"]].tolist()

# %%
replisome_monomer_subunit_ids = bulk_name_to_idx(
    ["CPLX0-3621[c]", "EG10239-MONOMER[c]", "EG11500-MONOMER[c]", "EG11412-MONOMER[c]"],
    bulk_molecule_ids,
)

replisome_trimer_subunit_ids = bulk_name_to_idx(
    ["CPLX0-2361[c]", "CPLX0-3761[c]"], bulk_molecule_ids
)

# %%
full_plasmids = sim_data.internal_state.unique_molecule.unique_molecule_definitions.get(
    "full_plasmid"
)
x = full_plasmids["_entryState"].sum()
