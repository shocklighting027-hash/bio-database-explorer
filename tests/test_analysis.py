import pytest

from bio_database_explorer.analysis import gc_content, motif_summary


def test_gc_content() -> None:
    assert gc_content("AAGC") == pytest.approx(0.5)


def test_gc_content_rejects_invalid_symbols() -> None:
    with pytest.raises(ValueError, match="Unsupported"):
        gc_content("ACGTX")


def test_motif_summary_finds_all_overlapping_occurrences() -> None:
    result = motif_summary("DWDW", motif="DW", flank=1)

    assert result.count == 2
    assert result.positions == (0, 2)
    assert result.contexts == ("DWD", "WDW")


def test_motif_summary_handles_missing_protein() -> None:
    assert motif_summary("").count == 0

