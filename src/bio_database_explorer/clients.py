"""Small clients for the NCBI and UniProt public APIs."""

from __future__ import annotations

import io
import time
from dataclasses import dataclass

import requests
from Bio import Entrez, SeqIO
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


@dataclass(frozen=True)
class NucleotideRecord:
    accession: str
    title: str
    sequence: str


@dataclass(frozen=True)
class UniProtRecord:
    accession: str
    entry_name: str
    protein_name: str
    gene: str
    length: int
    sequence: str


def _retrying_session() -> requests.Session:
    retry = Retry(
        total=4,
        backoff_factor=0.6,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset({"GET"}),
    )
    session = requests.Session()
    session.headers["User-Agent"] = "bio-database-explorer/0.1"
    session.mount("https://", HTTPAdapter(max_retries=retry))
    return session


class NCBIClient:
    """Fetch human MANE Select transcripts through Biopython Entrez."""

    def __init__(self, email: str, api_key: str | None = None) -> None:
        if not email or "@" not in email:
            raise ValueError("A valid contact email is required by NCBI")
        Entrez.email = email
        Entrez.api_key = api_key
        Entrez.tool = "bio-database-explorer"
        self.delay = 0.12 if api_key else 0.4

    def search_mane_transcripts(
        self,
        gene_query: str,
        organism: str = "Homo sapiens",
        limit: int = 100,
    ) -> list[str]:
        term = (
            f"{gene_query}[Gene Name] AND {organism}[Organism] "
            "AND mRNA[Filter] AND MANE Select[Keyword]"
        )
        with Entrez.esearch(db="nucleotide", term=term, retmax=limit) as handle:
            result = Entrez.read(handle)
        return list(result.get("IdList", []))

    def fetch_nucleotide(self, record_id: str) -> NucleotideRecord:
        time.sleep(self.delay)
        with Entrez.efetch(
            db="nucleotide",
            id=record_id,
            rettype="fasta",
            retmode="text",
        ) as handle:
            record = SeqIO.read(io.StringIO(handle.read()), "fasta")
        return NucleotideRecord(
            accession=record.id,
            title=record.description,
            sequence=str(record.seq),
        )


class UniProtClient:
    """Fetch reviewed human protein records from the UniProt REST API."""

    endpoint = "https://rest.uniprot.org/uniprotkb/search"

    def __init__(self, timeout: float = 20.0) -> None:
        self.timeout = timeout
        self.session = _retrying_session()

    def fetch_human_gene(self, gene: str) -> UniProtRecord | None:
        params = {
            "query": f"(gene_exact:{gene}) AND (organism_id:9606) AND (reviewed:true)",
            "fields": "accession,id,protein_name,gene_primary,length,sequence",
            "format": "tsv",
            "size": 1,
        }
        response = self.session.get(self.endpoint, params=params, timeout=self.timeout)
        response.raise_for_status()
        lines = response.text.strip().splitlines()
        if len(lines) < 2:
            return None

        values = lines[1].split("\t")
        if len(values) != 6:
            raise ValueError("UniProt returned an unexpected table format")
        accession, entry_name, protein_name, primary_gene, length, sequence = values
        return UniProtRecord(
            accession=accession,
            entry_name=entry_name,
            protein_name=protein_name,
            gene=primary_gene,
            length=int(length),
            sequence=sequence,
        )

