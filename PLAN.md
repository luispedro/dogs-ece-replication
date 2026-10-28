# ECE replication analysis: operational plan

Operationalise the analysis as an **evidence-based catalogue of replication systems**, followed by comparisons between ECE types and discovery of recurring replication modules. The initial analysis can use the assembled ECE sequences alone; raw-read analysis should be a targeted follow-up rather than a prerequisite.

## 1. Start with the compact ECE dataset

From the [Zenodo record](https://zenodo.org/records/16356977), the essential inputs are:

- `SHD1_EC.fna.gz`: the 185 non-redundant circular elements.
- `SHD1_EC_props.tsv.gz`: identifiers, source samples/contigs, categories, and notes.

The metadata distinguishes **58 plasmids, 30 viruses, and 97 uncategorized elements**. Six uncategorized elements have markers for both plasmids and viruses; retain these as a separate subgroup in comparisons, without automatically calling them phage–plasmids.

First, validate identifiers and compute length, GC content, ambiguous-base counts, and sequence checksums. Preserve the original categories throughout.

**Important caveat:** the [original methods](https://www.biorxiv.org/content/10.1101/2025.09.17.676595v1.full-text) assigned plasmid labels using geNomad hallmark genes. Consequently, enrichment of replication markers among “plasmids” may partly reflect how that category was defined, rather than an independent biological finding.

## 2. Define what we can infer

Separate three levels of inference:

| Level | Question | Example output |
|---|---|---|
| Protein/module | What replication-associated machinery is encoded? | Initiator, helicase, primase, polymerase |
| Replicon family | Does it resemble a characterized replication system? | A particular Rep family or RNA-based origin |
| Mechanism | What replication process does the evidence support? | Rolling-circle, theta, strand-displacement, unresolved |

For viruses, allow multiple candidate mechanisms: some switch replication modes during their life cycle. “Viral replication module detected” is useful even when the precise mechanism cannot be assigned.

Every assignment should include its supporting genes/regions, reference matches, and an explicit confidence tier. **No recognizable Rep protein does not mean no replication**: host-dependent and RNA-initiated systems exist.

## 3. Build the evidence pipeline

### A. Predict genes, accounting for circular boundaries

Use **Pyrodigal/Prodigal** in metagenomic mode, with explicit recovery of genes spanning the arbitrary FASTA boundary—for example, by calling an alternative rotation and reconciling coordinates and duplicates.

Output proteins, coding sequences, and coordinates. Circular-boundary handling is particularly important here because missing one initiator could change the interpretation of an entire element.

### B. Annotate replication and neighbouring functions

Use two complementary searches:

- **HMMER against pinned protein-family profiles**, including relevant Pfam and viral profiles.
- **DIAMOND against a curated reference set** of experimentally characterized replication proteins and systems.

Search for:

- Replication initiators and rolling-circle nickases.
- Helicases, primases, polymerases, and replication-associated accessory proteins.
- Partitioning and maintenance proteins.
- Relaxases and conjugation machinery.
- Viral structural and packaging proteins.

Keep replication, transfer, and maintenance annotations distinct. In particular, **a relaxase hit alone is not evidence of vegetative rolling-circle replication**.

The main scientific investment is a small, auditable mapping from reference families to mechanisms, recording literature support, thresholds, coverage requirements, and known ambiguities.

### C. Investigate candidate origins

For elements with candidate replication genes, inspect surrounding noncoding regions for:

- Direct repeats/iterons.
- Inverted repeats and potential stem–loops.
- AT-rich regions.
- Characterized origin sequences or RNA-based replication regions.

Treat these as supporting evidence, not stand-alone mechanism classifiers. Repeats and AT-rich regions are common; assess their enrichment against composition-matched background sequences.

### D. Combine evidence transparently

Start with a rule-based interpretation rather than a trained classifier:

- **Strong support:** a well-characterized replication-system match with compatible domain architecture and genomic context.
- **Tentative:** family-level or incomplete evidence.
- **Conflicting:** incompatible or multiple candidate systems.
- **Unresolved:** insufficient evidence.

Validate the rules on characterized reference elements, including difficult negatives such as transfer-associated nickases, and check sensitivity to annotation thresholds.

## 4. Compare types and discover patterns

Produce an **ECE × replication-feature matrix**, then ask:

1. Which systems occur in each original ECE category?
2. Does mechanism or module composition vary with element size?
3. Do the dual-marker elements combine viral and plasmid-like modules?
4. Do uncategorized elements share replication modules with labelled elements?
5. Are there recurring, currently uncharacterized gene neighbourhoods?

Use **MMseqs2** to cluster candidate replication proteins and compare their neighbouring genes. This provides a discovery route beyond known family annotations.

For comparisons, report counts, effect sizes, and uncertainty; use exact/permutation tests where counts are small and correct for multiple testing. Check whether results persist after grouping related ECEs, since non-redundancy does not imply evolutionary independence.

Suggested Matplotlib figures:

- Mechanism/support proportions by ECE category, including unresolved cases.
- Clustered replication-feature heatmap.
- Element length versus candidate mechanism.
- Gene-neighbourhood diagrams for representative modules.
- Replication-protein similarity network, coloured by ECE category.

## 5. Make it reproducible with Pixi and Jug

A minimal proposed structure:

```text
pixi.toml
pixi.lock
config/
    analysis.yaml
    replication_rules.yaml
data/
    manifest.tsv
jugfile.py
src/ece_replication/
tests/
results/
figures/
```

Pixi would manage Python, Jug, Matplotlib, pandas, SciPy, Pyrodigal, HMMER, DIAMOND, and MMseqs2. Reference databases need their own version/checksum manifest: the software lockfile does not pin them. All Python and analysis commands should run through the Pixi environment, not system Python.

The Jug task graph would be:

```text
download + verify inputs/references
                  ↓
        validate catalogue
                  ↓
     predict genes per element
                  ↓
      profile/reference searches
                  ↓
    origin/context evidence per element
                  ↓
         evidence assignments
                  ↓
 clustering → comparisons → plots
```

Use per-element tasks where sensible and batch searches where tool overhead dominates. Task arguments should include input checksums, configuration, database versions, and an explicit analysis version: Jug does not automatically invalidate cached results when implementation code changes.

## Recommended first milestone

**Get all 185 elements through gene prediction, targeted annotation, and an auditable evidence table before adding elaborate origin prediction or downloading raw reads.**

That first pass should deliver:

- One row per ECE with candidate systems, confidence, and unresolved/conflicting flags.
- A gene-level evidence table.
- Three initial figures: category proportions, feature heatmap, and representative neighbourhoods.
- A shortlist of unresolved or unusual elements for deeper investigation.

The initial conclusion would be **“these sequences support these replication-system hypotheses,”** not direct proof of replication dynamics. Targeted read validation or experiments can then address the most interesting hypotheses.
