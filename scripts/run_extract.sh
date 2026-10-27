#!/usr/bin/env bash
# Extract chains for all splits into data/chains/<split>/, skipping existing files.
# Files >100 MB run with fewer workers so RAM doesn't run out.
# Usage: scripts/run_extract.sh
cd "$(dirname "$0")/.."
PY="${PYTHON:-.venv/bin/python}"
LOG=data/extract_errors.log
: > "$LOG"
for big in 0 1; do
  for s in train eval test; do
    while read -r e; do
      [ -z "$e" ] && continue
      f="data/structures/$s/pdb/${e%_*}.cif"
      size=$(stat -c %s "$f" 2>/dev/null || echo 0)
      if (( (size > 100000000) == big )); then echo "$s $e"; fi
    done < "splits/${s}_pdb.txt"
  done | xargs -P $([ $big = 1 ] && echo 2 || echo 8) -L1 sh -c \
    '[ -f data/chains/$0/$1.cif ] || '"$PY"' scripts/extract_chain.py $1 data/structures/$0/pdb data/chains/$0 2>/dev/null || echo "$0 $1" >> '"$LOG"
done
for s in train eval test; do echo "$s: $(ls data/chains/$s | wc -l)"; done
echo "errors: $(wc -l < "$LOG")"
