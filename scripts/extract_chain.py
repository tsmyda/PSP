#!/usr/bin/env python3
# Usage:
#   python extract_chain.py <entry_id> <in_dir> <out_dir>
# Example:
#   python scripts/extract_chain.py 1abc_1 data/structures/train/pdb data/chains/train/
#
# Uses gemmi instead of Biopython: same selection logic, but large assemblies
# (ribosomes, capsids; up to ~500 MB mmCIF) fit in memory and parse much faster.
import sys
from pathlib import Path
import gemmi

entry_id, in_dir, out_dir = sys.argv[1:4]
pdb_id, entity_id = entry_id.rsplit("_", 1)
path = f"{in_dir}/{pdb_id}.cif"

doc = gemmi.cif.read(path)
block = doc.sole_block()
# Chains are named by auth_asym_id (pdb_strand_id), not label_asym_id
scheme = block.find("_pdbx_poly_seq_scheme.", ["pdb_strand_id", "entity_id"])
chains = {row.str(0) for row in scheme if row.str(1) == entity_id}

st = gemmi.make_structure_from_block(block)
del doc, block
for model in st:
    for name in {ch.name for ch in model} - chains:
        model.remove_chain(name)

Path(out_dir).mkdir(parents=True, exist_ok=True)
st.make_mmcif_document().write_file(f"{out_dir}/{entry_id}.cif")
