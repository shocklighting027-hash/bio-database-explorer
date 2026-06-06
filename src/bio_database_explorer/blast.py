"""NCBI BLAST helpers kept separate because remote BLAST can take minutes."""

from dataclasses import dataclass

from Bio.Blast import NCBIWWW, NCBIXML


@dataclass(frozen=True)
class BlastHit:
    title: str
    accession: str
    length: int
    e_value: float
    score: float
    identities: int
    alignment_length: int

    @property
    def identity_percent(self) -> float:
        return 100 * self.identities / self.alignment_length


def identify_sequence(
    sequence: str,
    database: str = "nt",
    program: str = "blastn",
    max_hits: int = 10,
    e_value: float = 1e-10,
) -> list[BlastHit]:
    """Submit a sequence to remote NCBI BLAST and return compact top hits."""
    normalized = "".join(sequence.split()).upper()
    if not normalized:
        raise ValueError("Sequence must not be empty")

    handle = NCBIWWW.qblast(
        program,
        database,
        normalized,
        hitlist_size=max_hits,
        expect=e_value,
    )
    try:
        record = NCBIXML.read(handle)
    finally:
        handle.close()

    hits: list[BlastHit] = []
    for alignment in record.alignments:
        if not alignment.hsps:
            continue
        hsp = alignment.hsps[0]
        hits.append(
            BlastHit(
                title=alignment.title,
                accession=alignment.accession,
                length=alignment.length,
                e_value=hsp.expect,
                score=hsp.score,
                identities=hsp.identities,
                alignment_length=hsp.align_length,
            )
        )
    return hits

