import pytest

from bio_database_explorer.pipeline import extract_gene_symbol


def test_extract_gene_symbol_from_ncbi_title() -> None:
    title = (
        "NM_152588.3 Homo sapiens transmembrane O-mannosyltransferase "
        "targeting cadherins 2 (TMTC2), transcript variant 1, mRNA"
    )
    assert extract_gene_symbol(title) == "TMTC2"


def test_extract_gene_symbol_requires_parenthesized_symbol() -> None:
    with pytest.raises(ValueError, match="Could not extract"):
        extract_gene_symbol("record without a gene symbol")

