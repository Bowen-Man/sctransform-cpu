"""Write ordered v2 fit-cell and fit-gene titles without fitting the model."""
from __future__ import annotations

import argparse
from pathlib import Path

from scipy import sparse

from sctransform_cpu import make_selection_manifest


def _lines(path: Path) -> list[str]:
    return path.read_text(encoding="utf-8").splitlines()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("data", type=Path, help="directory containing counts.npz and title files")
    parser.add_argument("output", type=Path)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--n-cells", type=int, default=2000)
    parser.add_argument("--n-genes", type=int, default=2000)
    parser.add_argument("--min-cells", type=int, default=5)
    args = parser.parse_args()
    manifest = make_selection_manifest(
        sparse.load_npz(args.data / "counts.npz"),
        gene_names=_lines(args.data / "gene_names.txt"),
        cell_names=_lines(args.data / "cell_names.txt"),
        seed=args.seed,
        n_cells=args.n_cells,
        n_genes=args.n_genes,
        min_cells=args.min_cells,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    manifest.write(args.output)


if __name__ == "__main__":
    main()
