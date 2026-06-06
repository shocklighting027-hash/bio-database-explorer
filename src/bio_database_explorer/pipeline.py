"""End-to-end gene-family exploration pipeline."""

from __future__ import annotations

import re

import pandas as pd

from .analysis import gc_content, motif_summary
from .clients import NCBIClient, UniProtClient

GENE_IN_PARENTHESES = re.compile(r"\(([A-Z][A-Z0-9-]+)\)")


def extract_gene_symbol(title: str) -> str:
    """Extract a gene symbol from an NCBI FASTA description."""
    matches = GENE_IN_PARENTHESES.findall(title.upper())
    if not matches:
        raise ValueError(f"Could not extract a gene symbol from: {title}")
    return matches[-1]


def explore_gene_family(
    family: str,
    email: str,
    api_key: str | None = None,
    motif: str = "DW",
) -> pd.DataFrame:
    """Collect MANE transcripts and reviewed proteins for a human gene family."""
    family = family.strip().upper().rstrip("*")
    if not family or not family.replace("-", "").isalnum():
        raise ValueError("Family prefix must contain letters, numbers, or hyphens")

    ncbi = NCBIClient(email=email, api_key=api_key)
    uniprot = UniProtClient()
    ids = ncbi.search_mane_transcripts(f"{family}*")

    rows: list[dict[str, object]] = []
    seen_genes: set[str] = set()
    for record_id in ids:
        nucleotide = ncbi.fetch_nucleotide(record_id)
        gene = extract_gene_symbol(nucleotide.title)
        if gene in seen_genes:
            continue
        seen_genes.add(gene)

        protein = uniprot.fetch_human_gene(gene)
        motif_result = motif_summary(protein.sequence if protein else "", motif=motif)
        rows.append(
            {
                "gene": gene,
                "ncbi_accession": nucleotide.accession,
                "nucleotide_length": len(nucleotide.sequence),
                "gc_fraction": gc_content(nucleotide.sequence),
                "uniprot_accession": protein.accession if protein else None,
                "protein_name": protein.protein_name if protein else None,
                "protein_length": protein.length if protein else None,
                "motif": motif,
                "motif_count": motif_result.count,
                "motif_positions_0_based": ",".join(map(str, motif_result.positions)),
                "motif_contexts": ";".join(motif_result.contexts),
                "nucleotide_sequence": nucleotide.sequence,
                "protein_sequence": protein.sequence if protein else None,
            }
        )

    columns = [
        "gene",
        "ncbi_accession",
        "nucleotide_length",
        "gc_fraction",
        "uniprot_accession",
        "protein_name",
        "protein_length",
        "motif",
        "motif_count",
        "motif_positions_0_based",
        "motif_contexts",
        "nucleotide_sequence",
        "protein_sequence",
    ]
    return pd.DataFrame(rows, columns=columns).sort_values("gene").reset_index(drop=True)

