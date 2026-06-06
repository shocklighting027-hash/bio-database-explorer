"""Command-line interface."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

from .blast import identify_sequence
from .pipeline import explore_gene_family


def _read_fasta(path: Path) -> str:
    lines = path.read_text(encoding="utf-8").splitlines()
    return "".join(line.strip() for line in lines if not line.startswith(">"))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="bio-db")
    subparsers = parser.add_subparsers(dest="command", required=True)

    family = subparsers.add_parser("family", help="Explore a human gene family")
    family.add_argument("prefix", help="Gene-family prefix, for example TMTC")
    family.add_argument("--motif", default="DW")
    family.add_argument("--email", default=os.getenv("NCBI_EMAIL"))
    family.add_argument("--api-key", default=os.getenv("NCBI_API_KEY"))
    family.add_argument("--output", type=Path, default=Path("results/gene_family.csv"))

    blast = subparsers.add_parser("blast", help="Identify a nucleotide FASTA with NCBI BLAST")
    blast.add_argument("fasta", type=Path)
    blast.add_argument("--max-hits", type=int, default=10)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.command == "family":
        if not args.email:
            raise SystemExit("Set NCBI_EMAIL or pass --email.")
        table = explore_gene_family(
            family=args.prefix,
            email=args.email,
            api_key=args.api_key,
            motif=args.motif,
        )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(args.output, index=False)
        summary = table.drop(columns=["nucleotide_sequence", "protein_sequence"])
        print(summary.to_string(index=False))
        print(f"\nSaved {len(table)} records to {args.output}")
        return

    sequence = _read_fasta(args.fasta)
    hits = identify_sequence(sequence, max_hits=args.max_hits)
    for index, hit in enumerate(hits, start=1):
        print(
            f"{index:>2}. {hit.accession} | {hit.identity_percent:.1f}% identity | "
            f"E={hit.e_value:.2g}\n    {hit.title}"
        )


if __name__ == "__main__":
    main()
