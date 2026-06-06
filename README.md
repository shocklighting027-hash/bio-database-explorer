# Bio Database Explorer

A small, reproducible bioinformatics project that identifies an unknown
nucleotide sequence and explores the corresponding human gene family using
public biological databases.

The original notebook started as a course assignment. This version turns it
into a reusable Python package with a command-line interface, tests, API error
handling, and a concise demonstration notebook.

## Research question

1. Which gene is the unknown nucleotide sequence associated with?
2. Which human genes belong to the same family?
3. What are their MANE Select transcript and reviewed protein sequences?
4. What is the GC content of each transcript?
5. Where does the `DW` amino-acid motif occur in each protein?

The original BLAST run identified the query as **TMTC2**. The family workflow
then retrieved human **TMTC1-TMTC4** records from NCBI and UniProt.

## Data sources

- [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/) for sequence identification
- [NCBI Nucleotide](https://www.ncbi.nlm.nih.gov/nuccore/) for MANE Select
  transcripts
- [UniProt](https://www.uniprot.org/) for reviewed human protein records

Database contents change over time, so accessions and returned records may
differ between runs.

## Project structure

```text
bio-database-explorer/
|-- data/
|   `-- mysterious_sequence.fasta
|-- notebooks/
|   `-- 01_gene_family_exploration.ipynb
|-- results/
|-- src/bio_database_explorer/
|   |-- analysis.py
|   |-- blast.py
|   |-- clients.py
|   |-- cli.py
|   `-- pipeline.py
|-- tests/
|-- .env.example
`-- pyproject.toml
```

## Quick start

Python 3.11 or newer is required.

```bash
python -m venv .venv
```

Activate the environment on Windows:

```powershell
.\.venv\Scripts\Activate.ps1
```

Install the project:

```bash
python -m pip install -e ".[notebook,dev]"
```

Set the email required by the
[NCBI Entrez usage policy](https://www.ncbi.nlm.nih.gov/books/NBK25497/):

```powershell
$env:NCBI_EMAIL="your.email@example.com"
```

Explore the TMTC family and save the result:

```bash
bio-db family TMTC --output results/tmtc_family.csv
```

Run remote BLAST for the included query:

```bash
bio-db blast data/mysterious_sequence.fasta --max-hits 5
```

Remote BLAST can take several minutes. Avoid repeatedly submitting the same
sequence during development.

## Notebook

Open `notebooks/01_gene_family_exploration.ipynb` after installing the
`notebook` dependencies. Network-dependent cells are clearly marked and are
not executed automatically.

## Output

The family pipeline creates one row per gene with:

- NCBI and UniProt accessions
- nucleotide and protein sequence lengths
- GC fraction
- all motif positions, reported as zero-based offsets
- local sequence context around each motif
- full nucleotide and protein sequences

Generated result files are ignored by Git because they can be recreated from
the source databases.

## Quality improvements over the original notebook

- reusable modules instead of duplicated notebook functions
- exact, reviewed, human-only UniProt queries
- request timeouts, retries, and NCBI rate limiting
- no personal email committed to source code
- all motif occurrences are counted
- pure analysis functions covered by tests
- automated linting and tests with GitHub Actions

## Tests

```bash
ruff check .
pytest
```

Tests do not call external APIs. This keeps the test suite fast and prevents
CI from putting unnecessary load on NCBI or UniProt.

## Limitations

- Gene symbols are extracted from the conventional parenthesized symbol in
  NCBI FASTA descriptions.
- The family search currently targets human MANE Select mRNA records.
- API responses are live and therefore not guaranteed to remain identical.
- BLAST results should be interpreted biologically rather than treated as a
  definitive annotation on their own.

## License

MIT

