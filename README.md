# PSP: Protein Structure Prediction

## Layout

```
├── scripts/
│   ├── prepare_splits.py       # cluster filtering, sampling, train/eval/test split
│   ├── download_structures.sh  # download mmCIF files from RCSB / AFDB
│   ├── extract_chain.py        # cut out the chains of one entity (gemmi)
│   ├── run_extract.sh          # run extract_chain.py for all splits in parallel
│   ├── make_fasta.py           # chain mmCIF -> FASTA (SEQRES or observed residues)
│   └── predict.py              # dummy predictor: FASTA -> CA-only PDB + PNG
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
uv sync --extra cpu     # or --extra cu128 for NVIDIA GPU; on macOS both give CPU + MPS
mkdir -p data
wget -P data https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-50.txt

uv run scripts/prepare_splits.py --team pdb           # writes splits/
for s in train eval test; do scripts/download_structures.sh $s pdb; done
scripts/run_extract.sh                                # writes data/chains/
```

## Predict

```bash
uv run scripts/make_fasta.py data/chains/test/*.cif --out-dir data/fasta/test
uv run scripts/predict.py data/fasta/test/10ZN_3.fasta --out-dir predictions/
```

`make_fasta.py` writes the full SEQRES sequence of the entity; `--observed` keeps only residues with coordinates (matches the ground-truth CA set).

`predict.py` reads one sequence from a FASTA file and writes `predictions/<id>.pdb` (one CA `ATOM` record per residue) and `predictions/<id>.png` (3D CA trace). The model is a placeholder that returns a random walk with 3.8 Å steps.

- Sequence is one-hot encoded into `(D, N)`: `D = --max-len` (default 1024), `N = 22` (20 amino acids + `X` + `PAD`). Longer sequences are truncated, shorter ones are padded and masked.
- `--device auto` uses CUDA, then MPS, then CPU. `--seed` makes the output reproducible.
