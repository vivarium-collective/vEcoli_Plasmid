from ecoli.library.schema import bulk_name_to_idx
import pandas as pd

# %%
# Loading Simdata to obtain the bulk id labels
from ecoli.library.sim_data import LoadSimData

sim_data_default = "out/plasmid/parca/kb/simData.cPickle"
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
rna_data = sim_data.process.transcription.rna_data

mrna_tu_ids = rna_data["id"][rna_data["is_mRNA"]].tolist()

cistron_data = sim_data.process.transcription.cistron_data

mrna_cistron_ids = cistron_data["id"][cistron_data["is_mRNA"]].tolist()

mrna_df = pd.DataFrame(mrna_cistron_ids, columns=["cistron_id"])

# %%
rnaA_id = bulk_name_to_idx(["PD03831[c]", "MONOMER0-160[c]"], bulk_molecule_ids)

# %%
nucleotide_pathways = {
    "adenosine deoxyribonucleotides de novo biosynthesis II": {
        "RXN0-747",
        "ADPREDUCT-RXN",
        "DADPKIN-RXN",
        "RXN0-745",
    },
    "guanosine deoxyribonucleotides de novo biosynthesis II": {
        "RXN0-748",
        "GDPREDUCT-RXN",
        "DGDPKIN-RXN",
        "RXN0-746",
    },
    "pyrimidine deoxyribonucleotides de novo biosynthesis II": {
        "RXN0-723",
        "DCTP-DEAM-RXN",
        "RXN0-724",
        "DUTP-PYROP-RXN",
        "THYMIDYLATESYN-RXN",
        "DTMPKI-RXN",
        "DTDPKIN-RXN",
    },
    "pyrimidine deoxyribonucleotides de novo biosynthesis I": {
        "RXN-12195",
        "CDPREDUCT-RXN",
        "DCDPKIN-RXN",
        "UDPREDUCT-RXN",
        "DUDPKIN-RXN",
        "DUTP-PYROP-RXN",
        "THYMIDYLATESYN-RXN",
        "DTMPKI-RXN",
        "DTDPKIN-RXN",
    },
}

keys = list(sim_data.process.metabolism.reaction_stoich.keys())


def base_id(k):
    return k.split(" (")[0]


grouped = {
    pathway: {rxn: i for i, k in enumerate(keys) if (rxn := base_id(k)) in reactions}
    for pathway, reactions in nucleotide_pathways.items()
}
# %%
replisome_subunits = bulk_name_to_idx(
    [
        "CPLX0-3621[c]",
        "EG10239-MONOMER[c]",
        "EG11500-MONOMER[c]",
        "EG11412-MONOMER[c]",
        "CPLX0-2361[c]",
        "CPLX0-3761[c]",
    ],
    bulk_molecule_ids,
)

# %%
adp_id = bulk_name_to_idx("ADP[c]", bulk_molecule_ids)
phosphatidylglycerol = bulk_name_to_idx("CPD-8260[c]", bulk_molecule_ids)
