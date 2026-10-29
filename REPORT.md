# ECE replication: first-milestone report

**Scope.** This report covers the first milestone in `PLAN.md`. All 185 circular ECEs went through gene prediction (with recovery of genes that cross the sequence start/end), targeted annotation, origin screening and rule-based evidence assignment. Raw reads were not used. The conclusions are of the form *"these sequences support these replication-system hypotheses"*. They do not show how the elements actually replicate.

Everything is reproducible with `pixi run fetch-refs && pixi run jug execute jugfile.py` (see `README.md`). Every number below comes from a file in `results/`.

---

## 1. Inputs and catalogue validation

* Inputs were `SHD1_EC.fna.gz` and `SHD1_EC_props.tsv.gz` from Zenodo 16356977. Both match the published MD5s, and the SHA-256s are pinned in `data/manifest.tsv`.
* There are 185 elements: 58 plasmids, 30 viruses, 91 uncategorized, and 6 uncategorized elements with both plasmid and virus markers. The original category is kept everywhere. The dual-marker six form their own subgroup (`uncategorized_dual_marker`). The FASTA and metadata identifiers match one-to-one, the FASTA headers agree with the metadata, and there are no ambiguous bases and no duplicate sequences.
* Lengths range from 1.3 to 425 kb (median 19.6 kb). GC ranges from 24% to 66% (median 37%). Median length is 30 kb for plasmids, 37 kb for viruses, 10.5 kb for uncategorized elements and 58 kb for the dual-marker subgroup (`results/catalogue.tsv`).
* The `_Size_` tag in the working header matches the sequence length for only 35/185 elements (median absolute difference 10 bp). It probably records the length before polishing or dereplication. It is recorded in the catalogue but not used.

**Caveat carried through all comparisons.** In the original study the "plasmid" label was assigned from geNomad plasmid hallmark genes, and those include relaxases and replication proteins. Plasmid-versus-other enrichments of exactly those features (most obviously relaxases: 84% of plasmids against 26% of all others) are therefore partly built into the labels and are not independent findings.

## 2. Pipeline summary

| Step | Implementation | Output |
|---|---|---|
| Gene prediction | Pyrodigal 3.7.1 in metagenomic mode on the sequence as given plus three rotations (0.5, 0.25, 0.75). Complete genes from a rotated pass are added when they cross the original start/end or replace a gene that was partial in the unrotated pass. Duplicates are resolved by stop codon. | 8,141 genes. 212 genes were partial in the unrotated pass. 132 were rescued from rotations, and 113 of those cross the start/end (`gene_prediction_stats.tsv`). |
| Profile search | HMMER 3.4 `hmmsearch` against Pfam 38.2 (gathering thresholds) and PHROGs v4. The PHROG HMMER profiles were rebuilt with `hmmbuild` from the published MSAs. A fixed `-Z` keeps E-values stable. | `work/hmm/*.domtbl` |
| Reference search | DIAMOND 2.2.8 (`--ultra-sensitive`) against 60 experimentally characterised proteins from UniProt 2026_03 (`config/reference_proteins.tsv`). These cover initiators, replisome proteins, partition proteins, relaxases as difficult negatives, and viral packaging proteins. | `gene_features.tsv` |
| Family → mechanism mapping | `config/replication_rules.yaml`. 166 Pfam families and 45 PHROG annotations map to a role (replication / maintenance / transfer / viral / integration), the mechanisms each can support, a specificity (high/medium/low), a literature note and known ambiguities. One DNA motif is included: the ColE1 RNAI/RNAII origin region. | |
| Origin screening | Up to 1 kb of intergenic DNA on each side of each candidate initiator, plus two distant control genes per element. Three statistics are computed: k-mer direct repeats (iteron-like), the longest inverted-repeat stem, and the maximum AT fraction. Each is tested against 200 first-order Markov backgrounds, with BH correction. | `origins.tsv` |
| Evidence combination | `evidence.py`, a rule-based combination into tiers: strong / tentative / conflicting / unresolved. | `element_evidence.tsv` (one row per ECE), `initiator_calls.tsv`, `genes.tsv` (gene-level evidence) |
| Discovery | MMseqs2 clustering of all proteins (30% identity, 50% coverage), initiator similarity network, initiator neighbourhoods, recurrent unannotated gene pairs, and unannotated clusters enriched in unresolved ECEs. | `results/discovery/` |
| Statistics | Permutation chi-square tests (20,000 permutations), Fisher tests with Haldane-corrected odds ratios, Wilson confidence intervals and BH q-values. Everything is repeated with one representative per group of related ECEs (protein-family Jaccard ≥ 0.3, giving 158 groups). | `results/comparisons/` |

**Tier rules (default thresholds):**

* **Strong:** a complete, high-specificity initiator. It needs corroboration from a characterised reference with the same mechanism (DIAMOND, ≥25% identity and ≥50% coverage) or from significant direct repeats next to it. The genome context must be compatible: an RCR initiator on an element above 25 kb is only accepted if a viral module is present.
* **Tentative:** family-level, incomplete or uncorroborated evidence, or a replication module whose mechanism cannot be assigned (`unassigned`).
* **Conflicting:** initiators with no mechanism in common.
* **Unresolved:** no replication-associated gene at all.

The rules also handle three special cases:

* A relaxase never counts as replication evidence.
* Rep_trans and Rol_Rep_N hits are ignored when the element carries a conjugation apparatus, because these families include the MobT relaxases of Tn916-like ICEs.
* Multi-mode viral initiators (lambda O-like: theta, then rolling circle) are not treated as conflicts.

## 3. Rule validation on characterised replicons

The same pipeline was run on 23 characterised replicons (`results/validation/`). Reference proteins from the element under test were excluded, so for example pT181 cannot be matched to its own RepC.

| Expected | Elements | Result (default thresholds) |
|---|---|---|
| Rolling circle (RCR) plasmids | pT181, pC194, pMV158, pE194, pUB110 | 5/5 RCR (3 strong) |
| RCR viruses | phiX174, M13, PCV2, P2 | 4/4 RCR (PCV2 strong) |
| Theta plasmids | F, RK2, R6K, R1, pSC101, pIP501, pAMβ1, pSM19035 | 8/8 theta (6 strong) |
| RNA-primed theta | ColE1 | theta, tentative (ColE1 origin motif + Rop; no Rep protein, as expected) |
| Strand displacement | RSF1010 | strand displacement, tentative (RepC + RepB′ primase). Its MobA relaxase is correctly not used. |
| Multi-mode virus | lambda | rolling circle\|theta |
| Difficult negative | Tn916 (excised ICE) | not called RCR. The MobT relaxase is flagged `rep_trans_in_transfer_context`. |
| Cointegrate | pSK41 | flagged **conflicting**: a RepA_N theta initiator plus Rep_1 RCR genes from integrated small plasmids |
| Exploratory | pBI143 (Bacteroides) | RepL-family initiator, mechanism unassigned |

Results: all 21 elements with a defined expectation were classified correctly (this includes lambda's dual mode and Tn916's expected "no vegetative replicon"). The pSK41 cointegrate was flagged as conflicting, and pBI143 is exploratory. The outcome is identical under the strict and relaxed thresholds. The relaxase-carrying difficult negatives (F TraI, RK2 TraI, pMV158 MobM, pE194 Pre, RSF1010 MobA, Tn916 MobT) never produced a replication call.

**This was a development set, not a held-out test.** Two rules were changed after seeing it:

* Rol_Rep_N was downgraded and flagged as a transfer look-alike because it matched the Tn916 MobT protein.
* RepB-RCR_reg (PF10723) was demoted because, despite its name, it hits the IncFII CopB proteins of the theta plasmids R1 and F.

A third change came from the ECE data rather than the validation set. Origin features were tested next to candidate initiators and next to control genes. Only direct repeats were enriched near initiators: 26% of initiator flanks against 7.5% of control flanks at P<0.05 (OR 4.3, Fisher P = 4×10⁻⁵). Inverted repeats (39% vs 41%) and AT-rich windows (40% vs 38%) were not enriched. Only direct repeats are therefore allowed to upgrade a call; hairpins and AT-rich windows are reported but not used.

Moving from the strict to the relaxed threshold set changes the number of strong ECE calls only from 34 to 42 (§4.2).

## 4. Results

### 4.1 Overall catalogue (`fig1_category_proportions.png`)

| Subgroup | n | Theta | Rolling circle | RC\|theta | Unassigned module | None detected | Strong | Tentative | Unresolved |
|---|---|---|---|---|---|---|---|---|---|
| Plasmid | 58 | 26 | 3 | 1 | 15 | 13 | 19 | 26 | 13 |
| Virus | 30 | 5 | 8 | 0 | 13 | 4 | 4 | 22 | 4 |
| Uncategorized | 91 | 18 | 2 | 0 | 29 | 42 | 14 | 35 | 42 |
| Dual marker | 6 | 4 | 1 | 1 | 0 | 0 | 2 | 4 | 0 |
| **All** | 185 | 53 | 14 | 2 | 57 | 59 | 39 | 87 | 59 |

No element was conflicting under the default rules. Only 68/185 ECEs (37%) received a specific mechanism hypothesis; another 57 carry replication genes whose mechanism cannot be assigned.

### Q1. Which systems occur in each category?

Mechanism class depends on category: permutation χ² = 62.8, P = 4.5×10⁻⁴, Cramér's V = 0.34. The result holds with one ECE per related group (n = 158, P = 2×10⁻⁴, V = 0.39).

* **Plasmids are dominated by theta-type initiators, overwhelmingly iteron-binding Rep_3 (Rep3_N/Rep3_C).** Theta calls cover 26/58 plasmids (45%, 95% CI 33–58%) against 21% of other ECEs (OR 3.0, q = 0.007). Rep_3 is by far the most common initiator family in the whole catalogue: its largest similarity component contains 50 proteins from all four subgroups (`fig5_replication_network.png`). RCR plasmids are rare: 3/58. These are small elements with Rep_1 or HUH (REP_ORF2-G2P) initiators, e.g. SHD1_EC.005, a 4.4 kb plasmid with a pC194-family Rep at 48% identity to pC194 Rep.
* **Viruses are enriched for RCR initiators:** 8/30 (27%, CI 14–44%) against 4% of other ECEs, OR 8.7, q = 0.002. With one ECE per related group the OR is 9.5 and q = 0.002. The viral RCR set includes:
  * a CRESS-DNA-like *Viral_Rep* element (SHD1_EC.087, 9.2 kb; similar to PCV1 Rep);
  * small Rep_1/HUH elements (SHD1_EC.069, SHD1_EC.085);
  * a recurrent cluster of five 35–46 kb phages (PC01786; SHD1_EC.063/065/067/068/079). Each encodes an HUH-endonuclease Rep (REP_ORF2-G2P domain, the domain shared by phiX174 and P2 gene A) next to a DNA polymerase and SSB (`fig4_neighbourhoods.png`, top row). Their genome size and replication module are consistent with P2-like RCR replication of a dsDNA phage. All five are tentative.
* The largest single group of viruses (13/30) carries replisome components (helicase, primase, polymerase, SSB) but no recognisable initiator. That is a "viral replication module detected, mechanism unassigned" result, as the plan anticipated.
* **Five "viruses" carry Rep_3 theta initiators** (SHD1_EC.060, .074, .076, .077, .086). Three of them (SHD1_EC.060, .074, .076) also carry ParA-like partition genes. These are phage-plasmid candidates.
* No RSF1010-type strand-displacement replicon and no protein-primed system was found in the ECEs.

### Q2. Does mechanism or module composition vary with size? (`fig3_length_vs_mechanism.png`)

Yes, strongly (Kruskal–Wallis P = 5×10⁻⁹):

| Class | Median length |
|---|---|
| No replication gene detected | 9.3 kb |
| Rolling circle | 11.6 kb |
| Theta | 38.0 kb |
| Unassigned module | 34.6 kb |

* Theta elements are larger than RCR elements (rank-biserial 0.40, q = 0.036). RCR calls fall into two groups: small plasmids and ssDNA-virus-like elements under about 15 kb, and the 35–46 kb phage cluster.
* Elements with no replication gene detected are mostly small (75% are under 13.5 kb). For these, the absence of a recognised initiator is more likely to reflect an unknown small Rep family than missing replication (see Q5).
* Module richness grows with size. The number of distinct module types correlates with log length at Spearman ρ = 0.71; replication roles alone at ρ = 0.65; maintenance roles at ρ = 0.69.

### Q3. Do the dual-marker elements combine viral and plasmid-like modules? (`comparisons/dual_marker_elements.tsv`)

**Mostly not, at the level of replication.**

* All six have a mechanism-assigned replicon: four theta (Rep_3), one RCR (Rep_1; SHD1_EC.091, 5.9 kb), and one lambda O-like multi-mode initiator (SHD1_EC.093, 94 kb).
* All six carry plasmid-like transfer or maintenance genes: relaxase, T4CP/T4SS or ParA.
* Only one, SHD1_EC.089 (67 kb), also carries detectable viral structural genes (tail).
* SHD1_EC.093 is the most phage-plasmid-like replicon: a lambda O-type phage initiator together with ParA, T4SS and relaxase genes, but no detected viral structural module.

So most dual-marker elements look like plasmids carrying isolated virus-associated genes, rather than full phage–plasmid hybrids. SHD1_EC.089 and SHD1_EC.093 are the best phage–plasmid candidates, together with the Rep_3-carrying "viruses" listed under Q1.

### Q4. Do uncategorized elements share replication modules with labelled elements?

**Partly.**

* 28/91 uncategorized ECEs (excluding dual-marker) carry an initiator gene. Nine of them share an initiator protein cluster (≥30% identity) with a plasmid or virus. In addition, eight multi-element initiator clusters mix uncategorized or dual-marker elements with labelled ones (`discovery/initiator_clusters.tsv`). Most of these shared clusters are Rep_3; the PC01786 phage HUH-Rep cluster also contains one uncategorized element (SHD1_EC.148).
* At the family level, the uncategorized elements use the same systems as the labelled ones: Rep_3 theta (18) and a few RCR (2). The difference is that they are much more often **unresolved**: 42/91 (46%) against 19% of other ECEs, OR 3.8, q = 5×10⁻⁴; with one ECE per group OR 5.0, q = 2×10⁻⁴.
* They are also small: median 10.5 kb.

### Q5. Recurrent uncharacterised neighbourhoods and candidate novel initiators

* **PC03265, a candidate unrecognised initiator family.** This 440–600 aa protein cluster has no Pfam or PHROG annotation above threshold. It occurs in 17 ECEs: 15 unresolved, 2 resolved (Fisher P = 7×10⁻⁷, q = 6×10⁻⁶; `discovery/unresolved_enriched_clusters.tsv`).
  * Its hosts are small (5–13 kb) mobilisable plasmids and uncategorized elements. Most of them carry a MobA/VirD2-type relaxase plus a NikA-like protein (`fig4_neighbourhoods.png`, row 7) and nothing else recognisable.
  * Four members have sub-threshold hits (E ≈ 10⁻⁴, HMM positions 34–90) to Pfam **Replicase** (PF03090), the iteron-binding RepA family related to Rep_3.
  * The flanking DNA has iteron-like direct repeats in 7/28 regions (q < 0.1; 6 of the 7 upstream of the gene), against 2/25 control regions in the same elements. That difference is suggestive but not significant (Fisher P = 0.15).
  * Together, this makes PC03265 the strongest candidate for a divergent Rep_3/Replicase-like theta initiator. It alone would resolve roughly a quarter of the unresolved ECEs.
* **Three small unannotated clusters (PC01009, PC03796, PC03802; 100–120 aa)** occur together in unresolved plasmids and uncategorized elements (e.g. SHD1_EC.039, .127, .130). They are weaker candidates: q = 0.06, with no origin signal.
* **27 recurrent modules of adjacent unannotated genes** (each pair present in at least 3 ECEs) (`discovery/recurrent_uncharacterised_modules.tsv`), are mostly specific to the Bacteroidales-like Rep_3 plasmid families, e.g. UM002/UM003 in 6 plasmids. UM013 links three unresolved plasmids that also carry PC03265.

**Direct DNA similarity to characterised replicons** (`discovery/nucleotide_similarity_to_characterised_replicons.tsv`):

* SHD1_EC.055 (97.5% identity over R1's RepA region) and SHD1_EC.026 (F/R1-like regions) are Enterobacterial IncF-type plasmids.
* SHD1_EC.000 matches ColE1 over 2.5 kb (85% identity) and carries the ColE1 RNAI/RNAII origin motif. It is the only ECE assigned RNA-primed, ColE1-type theta replication.
* SHD1_EC.131 and SHD1_EC.175 share Tn916-like regions. Consistently, their Rep_trans proteins sit in a conjugative context and are flagged as probable MobT relaxases.

### 4.2 Robustness

* **Thresholds** (`sensitivity_counts.tsv`, `figS1_threshold_sensitivity.png`): 161/185 elements keep the same mechanism and tier under the strict, default and relaxed sets. Changes are mostly between strong and tentative, plus a few unassigned calls that drop to unresolved under strict thresholds.
* **Relatedness:** every category effect reported above has the same direction and similar magnitude with one ECE per related group (`comparisons/*_grouped.tsv`).

## 5. Shortlist for deeper investigation (`results/shortlist.tsv`)

38 elements are marked high priority; 88 more are routine, mostly "replisome components without initiator".

1. **PC03265 carriers** (17 ECEs, e.g. SHD1_EC.009, .011, .012, .028, .042, .052, .096, .103, .104, .181). Next steps: test the candidate initiator with structure prediction and remote homology (HHpred/Foldseek against Rep_3/RepA structures), and map iterons upstream of it.
2. **Phage-plasmid candidates:**
   * viruses with partition genes and/or Rep_3: SHD1_EC.060, .065, .067, .068, .074, .076, .080, .081, .082;
   * dual-marker elements: SHD1_EC.089 and SHD1_EC.093 (lambda O-like initiator, no viral structural genes).
3. **The 35–46 kb HUH-Rep phage cluster (PC01786).** These are candidates for RCR replication of a dsDNA phage. Read-level coverage or strand-asymmetry analysis could test this.
4. **MobT versus RCR:** SHD1_EC.131 and SHD1_EC.175 (Tn916-like).
5. **Viruses with structural modules but no assigned replication system,** e.g. SHD1_EC.058, .059, .066, .071 to .075, .083. These are probably large phages that use unrecognised or host-dependent initiation.

## 6. Limitations

* **Category labels are circular** for the plasmid-hallmark features (see §1).
* **No initiator detected is not evidence of no replication.** RNA-primed systems are only detectable through the one ColE1 motif included, and host-dependent or novel initiators (e.g. PC03265) are invisible to family searches.
* **The DIAMOND reference set is small and biased** towards Firmicutes and Enterobacteria plasmids. Only 12 of the 39 strong calls rest on a characterised reference; the rest rest on iteron-like repeats next to a complete Rep_3 initiator.
* **The validation set guided two rule changes** (§3) and is small. A held-out set, for example from PLSDB with literature-curated replication types, would be the next check.
* **Origin statistics are crude:** exact k-mer repeats, perfect stems, and first-order Markov backgrounds. Only direct repeats behaved as informative in this dataset.
* **Mechanism remains a hypothesis.** Testing it needs read-level data (e.g. coverage skew around candidate origins, ssDNA intermediates) or experiments.

## 7. Outputs

| File | Content |
|---|---|
| `results/element_evidence.tsv` | One row per ECE: mechanism, tier, supporting genes, initiator families, reference and origin support, module roles, flags, related-group id |
| `results/genes.tsv` | Gene-level table: coordinates (including `wraps_origin`), best Pfam/PHROG, mapped roles, evidence text, protein cluster |
| `results/gene_features.tsv` | Every rule-mapped hit for each threshold set (source, family, score, coverage, identity) |
| `results/initiator_calls.tsv` | Gene-level initiator summaries |
| `results/origins.tsv`, `results/origin_motifs.tsv` | Origin statistics for candidate and control flanks; ColE1 motif hits |
| `results/feature_matrix.tsv` | ECE × replication-feature binary matrix |
| `results/comparisons/` | Category tests (all ECEs and one per related group), length tests, dual-marker table |
| `results/discovery/` | Initiator clusters, neighbourhoods, conserved neighbours, recurrent unannotated modules, unresolved-enriched clusters and their origin statistics, nucleotide similarity to characterised replicons |
| `results/validation/` | Validation evidence and evaluation |
| `results/sensitivity_*.tsv` | Threshold sensitivity |
| `results/shortlist.tsv` | Prioritised elements for follow-up |
| `results/provenance/` | Tool versions and input checksums |
| `figures/fig1`–`fig5`, `figS1` | Category proportions, clustered feature heatmap, length vs. mechanism, representative neighbourhoods, initiator similarity network, threshold sensitivity (PNG + PDF) |
