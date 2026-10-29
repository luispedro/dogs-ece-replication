# ECE replication catalogue

Evidence-based catalogue of replication systems in the 185 circular
extrachromosomal elements (ECEs) of the Shanghai dog gut long-read catalogue
([Zenodo 16356977](https://zenodo.org/records/16356977)). See `PLAN.md` for
the design and `REPORT.md` for results.

## Running

Everything runs inside the Pixi environment (never the system Python):

```bash
pixi install
pixi run fetch-refs                 # download inputs + references, pin SHA-256 in data/manifest.tsv
pixi run test                       # unit tests (gene boundary handling, origin stats, evidence rules)
pixi run jug execute jugfile.py &   # start one or more workers
pixi run jug execute jugfile.py &
pixi run jug status jugfile.py
```

`fetch-refs` is idempotent: once a checksum is recorded in `data/manifest.tsv`
it is treated as pinned and a changed download is an error. The first Jug task
re-verifies all files against the manifest.

Jug caches results by task arguments, not code. Every task receives
`analysis_version` from `config/analysis.yaml` (and the checksums of the files
it reads, including `config/replication_rules.yaml`), so bump the version
after changing code.

## Layout

| Path | Contents |
|---|---|
| `config/analysis.yaml` | paths, tool parameters, seeds, analysis version |
| `config/replication_rules.yaml` | family → role/mechanism mapping, specificities, ambiguities, threshold sets, combination rules |
| `config/reference_proteins.tsv` | curated characterised replication/maintenance/transfer proteins (DIAMOND reference set) |
| `config/validation_elements.tsv` | characterised replicons used to check the rules (incl. difficult negatives) |
| `data/manifest.tsv` | URLs, versions, published MD5s and pinned SHA-256s of every input/reference |
| `jugfile.py` | task graph |
| `src/ece_replication/` | `genes` (circular-aware Pyrodigal), `search` (HMMER/DIAMOND/MMseqs2), `annotate`, `origins`, `evidence`, `cluster`, `compare`, `plots`, `analysis` |
| `results/` | output tables (see `REPORT.md`) |
| `figures/` | Matplotlib figures (PNG + PDF) |

## References used

* Pfam 38.2 (`Pfam-A.hmm.gz`, gathering thresholds from `Pfam-A.hmm.dat.gz`)
* PHROGs v4 — HMMER3 profiles rebuilt with `hmmbuild` from the published MSAs
  (the PHROG server only distributes HH-suite/MMseqs2 profiles; its TLS
  certificate is invalid, so integrity rests on the pinned SHA-256)
* 60 characterised proteins from UniProt (release recorded in the manifest)
* 23 characterised replicons from NCBI nuccore (versioned accessions)
