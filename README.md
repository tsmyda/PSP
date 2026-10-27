# PSP: Protein Structure Prediction

## Layout

```
├── scripts/
│   ├── prepare_splits.py       # cluster filtering, sampling, train/eval/test split
│   ├── download_structures.sh  # download mmCIF files from RCSB / AFDB
│   ├── extract_chain.py        # cut out the chains of one entity (gemmi)
│   └── run_extract.sh          # run extract_chain.py for all splits in parallel
├── splits/                     # protein IDs per split (tracked in git)
│   ├── {train,eval,test}_pdb.txt
│   └── excluded.txt            # IDs removed after sampling (no atomic coordinates)
└── data/                       # large files, not tracked (see .gitignore)
    ├── clusters-by-entity-50.txt
    ├── structures/<split>/pdb/ # full PDB entries (~26 GB)
    └── chains/<split>/         # one mmCIF per protein (~2.7 GB)
```

## Dataset

| Split | Clusters | Proteins |
|-------|----------|----------|
| train | 400      | 3998     |
| eval  | 50       | 499      |
| test  | 50       | 499      |

Clusters come from RCSB at 50% sequence identity. Whole clusters are assigned to one split, so no test protein has a homolog above 50% identity in the training set.

## Reproduce

Requires [uv](https://docs.astral.sh/uv/) and `wget`.

```bash
uv sync
mkdir -p data
wget -P data https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-50.txt

uv run scripts/prepare_splits.py --team pdb           # writes splits/
for s in train eval test; do scripts/download_structures.sh $s pdb; done
scripts/run_extract.sh                                # writes data/chains/
```
