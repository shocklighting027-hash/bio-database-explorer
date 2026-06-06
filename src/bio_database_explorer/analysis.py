"""Pure sequence-analysis helpers."""

from dataclasses import dataclass

from Bio.SeqUtils import gc_fraction


@dataclass(frozen=True)
class MotifSummary:
    """Locations and sequence contexts for a protein motif."""

    motif: str
    count: int
    positions: tuple[int, ...]
    contexts: tuple[str, ...]


def gc_content(sequence: str) -> float:
    """Return the GC fraction in the inclusive range 0..1."""
    normalized = sequence.strip().upper()
    if not normalized:
        raise ValueError("Nucleotide sequence must not be empty")
    invalid = set(normalized) - set("ACGTUN")
    if invalid:
        raise ValueError(f"Unsupported nucleotide symbols: {sorted(invalid)}")
    return float(gc_fraction(normalized))


def motif_summary(protein: str, motif: str = "DW", flank: int = 4) -> MotifSummary:
    """Find every overlapping motif occurrence and its surrounding context."""
    if not protein:
        return MotifSummary(motif=motif, count=0, positions=(), contexts=())
    if not motif:
        raise ValueError("Motif must not be empty")
    if flank < 0:
        raise ValueError("Flank must be non-negative")

    protein = protein.upper()
    motif = motif.upper()
    positions: list[int] = []
    start = 0
    while True:
        position = protein.find(motif, start)
        if position == -1:
            break
        positions.append(position)
        start = position + 1

    contexts = tuple(
        protein[max(0, position - flank) : position + len(motif) + flank]
        for position in positions
    )
    return MotifSummary(
        motif=motif,
        count=len(positions),
        positions=tuple(positions),
        contexts=contexts,
    )

