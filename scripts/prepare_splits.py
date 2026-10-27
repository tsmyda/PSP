#!/usr/bin/env python3
"""
prepare_splits.py
Usage: python3 scripts/prepare_splits.py --team pdb   (for Team 1)
       python3 scripts/prepare_splits.py --team af    (for Team 2)

Reads data/clusters-by-entity-50.txt, writes splits/{train,eval,test}_<team>.txt.
"""

import argparse
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CLUSTERS = ROOT / "data" / "clusters-by-entity-50.txt"
SPLITS_DIR = ROOT / "splits"
# IDs removed after sampling (e.g. entries without atomic coordinates)
EXCLUDED = SPLITS_DIR / "excluded.txt"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Prepare dataset splits for PSP project."
    )
    parser.add_argument(
        "--team",
        choices=["pdb", "af"],
        required=True,
        help="Choose 'pdb' for Team 1 (PDB) or 'af' for Team 2 (AFDB)",
    )
    parser.add_argument(
        "--seed", type=int, default=42, help="Random seed for reproducibility"
    )
    return parser.parse_args()


def is_team_member(protein_id, team):
    if team == "pdb":
        return not (protein_id.startswith("AF_") or protein_id.startswith("MA_"))
    elif team == "af":
        return protein_id.startswith("AF_")
    return False


def main():
    args = parse_args()
    random.seed(args.seed)

    # 1. Parse clusters and filter entries by team source
    valid_clusters = []

    with open(CLUSTERS, "r") as f:
        for line in f:
            members = line.strip().split()
            # Keep only members belonging to the designated team source
            team_members = [m for m in members if is_team_member(m, args.team)]

            # Step 3 requirement: Select clusters with at least 10 valid members
            if len(team_members) >= 10:
                valid_clusters.append(team_members)

    print(f"Total valid clusters found for team '{args.team}': {len(valid_clusters)}")

    # Step 4 requirement: Shuffle clusters
    random.shuffle(valid_clusters)

    # Step 5 & 6 requirements: Sample 10 members per cluster and assign to 4000/500/500 splits
    # Each cluster contributes 10 members, so we need:
    # Train: 400 clusters = 4,000 proteins
    # Eval:   50 clusters =   500 proteins
    # Test:   50 clusters =   500 proteins
    total_needed_clusters = 400 + 50 + 50  # 500 clusters total

    if len(valid_clusters) < total_needed_clusters:
        print(
            f"Error: Not enough clusters available! Found {len(valid_clusters)}, need {total_needed_clusters}."
        )
        sys.exit(1)

    train_ids, eval_ids, test_ids = [], [], []

    for idx, cluster in enumerate(valid_clusters[:total_needed_clusters]):
        sampled_members = random.sample(cluster, 10)

        if idx < 400:
            train_ids.extend(sampled_members)
        elif idx < 450:
            eval_ids.extend(sampled_members)
        else:
            test_ids.extend(sampled_members)

    # Drop excluded IDs after sampling, so the remaining selection is unchanged
    excluded = set()
    if EXCLUDED.exists():
        excluded = {
            line.strip()
            for line in EXCLUDED.read_text().splitlines()
            if line.strip() and not line.startswith("#")
        }

    # Step 7 requirement: Save IDs to text files
    suffix = "pdb" if args.team == "pdb" else "af"

    splits = {
        f"train_{suffix}.txt": train_ids,
        f"eval_{suffix}.txt": eval_ids,
        f"test_{suffix}.txt": test_ids,
    }

    for filename, ids in splits.items():
        ids = [pid for pid in ids if pid not in excluded]
        with open(SPLITS_DIR / filename, "w") as f:
            for pid in ids:
                f.write(f"{pid}\n")
        print(f"Saved {len(ids)} IDs to {filename}")


if __name__ == "__main__":
    main()
