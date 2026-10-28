#!/usr/bin/env python3
# Usage:
#   python predict.py <input.fasta> [--out-dir DIR] [--max-len D] [--device auto|cpu|cuda|mps] [--seed S]
# Example:
#   uv run scripts/predict.py data/example.fasta --out-dir predictions/
#
# Writes <out-dir>/<id>.pdb (CA-only ATOM records) and <out-dir>/<id>.png (3D CA trace).
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import skbio
import torch
from torch import nn

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
UNK, PAD = "X", "<PAD>"
VOCAB = list(AMINO_ACIDS) + [UNK, PAD]
TOKEN_INDEX = {tok: i for i, tok in enumerate(VOCAB)}

THREE_LETTER = {
    "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE",
    "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU",
    "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG",
    "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR",
    UNK: "UNK",
}

CA_CA_DISTANCE = 3.8


def read_single_fasta(path):
    records = list(skbio.io.read(str(path), format="fasta", lowercase=True))
    if len(records) != 1:
        raise ValueError(f"{path}: expected exactly one sequence, found {len(records)}")
    record = records[0]
    # Gaps and stop codons are not residues
    residues = "".join(c for c in str(record).upper() if c not in "-.*")
    return record.metadata.get("id") or Path(path).stem, residues


def encode(sequence, max_len):
    """One-hot (D, N) encoding and (D,) mask; truncates to D, pads the rest with PAD."""
    residues = [c if c in TOKEN_INDEX else UNK for c in sequence[:max_len]]
    indices = [TOKEN_INDEX[c] for c in residues] + [TOKEN_INDEX[PAD]] * (max_len - len(residues))
    one_hot = nn.functional.one_hot(torch.tensor(indices), num_classes=len(VOCAB)).float()
    mask = torch.zeros(max_len, dtype=torch.bool)
    mask[: len(residues)] = True
    return one_hot, mask, "".join(residues)


class DummyStructureModel(nn.Module):
    """Ignores sequence content; returns a random walk with CA-CA spacing for valid residues."""

    def forward(self, x, mask):
        steps = torch.randn(x.shape[0], 3, device=x.device)
        steps = CA_CA_DISTANCE * steps / steps.norm(dim=-1, keepdim=True)
        coords = torch.cumsum(steps, dim=0)
        coords = coords[mask]
        return coords - coords.mean(dim=0)


def write_pdb(path, residues, coords, chain_id="A"):
    lines = []
    for i, (res, (x, y, z)) in enumerate(zip(residues, coords.tolist()), start=1):
        # wwPDB ATOM record, columns 1-80
        lines.append(
            f"{'ATOM':<6}{i:>5} {' CA ':<4}{'':1}{THREE_LETTER[res]:>3} {chain_id:1}{i:>4}{'':1}   "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}          {'C':>2}{'':2}"
        )
    last = len(residues)
    lines.append(f"{'TER':<6}{last + 1:>5}      {THREE_LETTER[residues[-1]]:>3} {chain_id:1}{last:>4}")
    lines.append("END")
    Path(path).write_text("\n".join(lines) + "\n")


def plot_trace(path, coords, title):
    xyz = coords.T
    fig = plt.figure(figsize=(7, 7))
    ax = fig.add_subplot(projection="3d")
    ax.plot(*xyz, color="lightgray", linewidth=1)
    points = ax.scatter(*xyz, c=range(len(coords)), cmap="viridis", s=12)
    fig.colorbar(points, ax=ax, shrink=0.6, pad=0.12, label="residue index")
    ax.set(xlabel="X [Å]", ylabel="Y [Å]", zlabel="Z [Å]", title=title)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def pick_device(name):
    if name != "auto":
        return torch.device(name)
    if torch.cuda.is_available():
        return torch.device("cuda")
    if torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def main():
    parser = argparse.ArgumentParser(description="Dummy CA-only structure predictor")
    parser.add_argument("fasta", type=Path)
    parser.add_argument("--out-dir", type=Path, default=Path("predictions"))
    parser.add_argument("--max-len", type=int, default=1024)
    parser.add_argument("--device", default="auto", choices=["auto", "cpu", "cuda", "mps"])
    parser.add_argument("--seed", type=int, default=None)
    args = parser.parse_args()

    if args.seed is not None:
        torch.manual_seed(args.seed)
    device = pick_device(args.device)

    name, sequence = read_single_fasta(args.fasta)
    if not sequence:
        raise ValueError(f"{args.fasta}: sequence is empty")
    x, mask, residues = encode(sequence, args.max_len)

    model = DummyStructureModel().to(device).eval()
    with torch.no_grad():
        coords = model(x.to(device), mask.to(device)).cpu()
    assert coords.shape == (len(residues), 3)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    pdb_path = args.out_dir / f"{name}.pdb"
    png_path = args.out_dir / f"{name}.png"
    write_pdb(pdb_path, residues, coords)
    plot_trace(png_path, coords.numpy(), f"{name} (L={len(residues)}, dummy prediction)")

    truncated = f", truncated from {len(sequence)}" if len(sequence) > args.max_len else ""
    print(f"{name}: L={len(residues)}{truncated}, device={device}")
    print(f"  {pdb_path}\n  {png_path}")


if __name__ == "__main__":
    main()
