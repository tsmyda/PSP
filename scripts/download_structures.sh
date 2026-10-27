#!/usr/bin/env bash
# Usage: scripts/download_structures.sh train pdb
#        scripts/download_structures.sh eval  af

SPLIT=$1   # train / eval / test
DB=$2      # pdb / af

cd "$(dirname "$0")/.."
OUT="data/structures/${SPLIT}/${DB}"
mkdir -p "$OUT"

if [[ "$DB" == "pdb" ]]; then
    awk -F'_' '{print "https://files.rcsb.org/download/"$1".cif"}' "splits/${SPLIT}_pdb.txt" \
        | sort -u | wget -nc -i - -P "$OUT"
elif [[ "$DB" == "af" ]]; then
    python3 -c "
import re
for line in open('splits/${SPLIT}_af.txt'):
    m = re.match(r'AF_AF(.+?)(F\d+)_\d+$', line.strip())
    if m:
        print(f'https://alphafold.ebi.ac.uk/files/AF-{m.group(1)}-{m.group(2)}-model_v6.cif')
" | sort -u | wget -nc -i - -P "$OUT"
fi
