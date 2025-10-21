# %%
import os
from reconstruction.spreadsheets import JsonWriter

# Path to flat dir
FLAT_DIR = "/Users/rashmidissasekara/Documents/code/vEcoli/reconstruction/ecoli/flat"
out_file = os.path.join(FLAT_DIR, "plasmid_dna_sites.tsv")

# Define rows
rows = [
    {
        "id": "P-ori",
        "common_name": "ori",
        "synonyms": ["origin of replication"],  # list, JsonWriter will encode
        "type": "origin-of-replication",
        "left_end_pos": 2534,
        "right_end_pos": 3122,
        "direction": None,  # will become JSON "null"
    }
]

# Write with JsonWriter
with open(out_file, "w", encoding="utf-8") as f:
    writer = JsonWriter(f, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)

# print(f"Wrote plasmid_dna_sites.tsv at {out_file}")
