#!/usr/bin/env python3
# Usage:
#   python make_fasta.py <chain.cif>... [--out-dir DIR] [--observed]
# Example:
#   uv run scripts/make_fasta.py data/chains/test/*.cif --out-dir data/fasta/test
#
# Writes one <entry_id>.fasta per input file. By default the full SEQRES sequence
# of the entity; with --observed only residues that have coordinates (first chain).
import argparse
from pathlib import Path

import gemmi


def entity_sequence(st, entity_id, observed):
    if observed:
        # extract_chain.py keeps only chains of this entity
        chain = st[0][0]
        return gemmi.one_letter_code([res.name for res in chain.get_polymer()]), chain.name
    entity = st.get_entity(entity_id)
    return gemmi.one_letter_code([gemmi.Entity.first_mon(m) for m in entity.full_sequence]), None


def main():
    parser = argparse.ArgumentParser(description="Extract FASTA from chain mmCIF files")
    parser.add_argument("cif", type=Path, nargs="+")
    parser.add_argument("--out-dir", type=Path, default=Path("data/fasta"))
    parser.add_argument("--observed", action="store_true", help="only residues with coordinates")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    for path in args.cif:
        entry_id = path.stem
        entity_id = entry_id.rsplit("_", 1)[1]
        st = gemmi.read_structure(str(path))
        st.setup_entities()
        sequence, chain = entity_sequence(st, entity_id, args.observed)
        header = f">{entry_id} entity={entity_id}" + (f" chain={chain}" if chain else "")
        (args.out_dir / f"{entry_id}.fasta").write_text(f"{header}\n{sequence.upper()}\n")
    print(f"{len(args.cif)} FASTA files -> {args.out_dir}")


if __name__ == "__main__":
    main()
