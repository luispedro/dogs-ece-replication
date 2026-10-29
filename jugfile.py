"""Jug task graph for the ECE replication catalogue.

Run (from the project root, always through pixi):

    pixi run fetch-refs                       # once: download + pin references
    pixi run jug execute jugfile.py &         # start as many workers as wanted
    pixi run jug status jugfile.py
"""
from jug import TaskGenerator, bvalue

from ece_replication import pipeline as P
from ece_replication.util import load_config, read_fasta

cfg = load_config()
V = cfg["analysis_version"]
paths, refs, inp = cfg["paths"], cfg["references"], cfg["inputs"]
S = cfg["search"]

verify_inputs = TaskGenerator(P.verify_inputs)
validate_catalogue = TaskGenerator(P.validate_catalogue)
validation_catalogue = TaskGenerator(P.validation_catalogue)
predict_genes = TaskGenerator(P.predict_genes)
collect_genes = TaskGenerator(P.collect_genes)
combine_proteins = TaskGenerator(P.combine_proteins)
hmmsearch_all = TaskGenerator(P.hmmsearch_all)
subset_hits = TaskGenerator(P.subset_hits)
diamond_refs = TaskGenerator(P.diamond_refs)
tool_versions = TaskGenerator(P.tool_versions)

checksums = verify_inputs(V)
versions = tool_versions(V)

# ---- catalogues ---------------------------------------------------------
ece_catalogue = validate_catalogue(inp["fasta"], inp["props"], inp["expected_counts"],
                                   inp["dual_marker_note"], checksums, V)
val_catalogue = validation_catalogue(refs["validation_fasta"], refs["validation_table"],
                                     checksums["validation_elements.fna"], V)

# Sequences are read at import time only to enumerate per-element tasks; each
# task also receives the verified checksum so it cannot run on unverified input.
DATASETS = {
    "ece": (inp["fasta"], "SHD1_EC.fna.gz"),
    "validation": (refs["validation_fasta"], "validation_elements.fna"),
}

gene_sets, annotations = {}, {}
for dataset, (fasta, key) in DATASETS.items():
    seqs = {name: seq.upper() for name, _, seq in read_fasta(fasta)}
    per_element = [predict_genes(name, seq, cfg["gene_prediction"], (V, checksums[key]))
                   for name, seq in seqs.items()]
    gene_sets[dataset] = collect_genes(dataset, per_element, paths["work_dir"], V)

# Profile searches are batched: reading the large HMM databases dominates the
# cost, so all proteins (ECE + validation) are searched in one pass per
# database. A fixed -Z keeps E-values independent of the batch composition.
combined = combine_proteins([gene_sets[d] for d in DATASETS], paths["work_dir"], V)
hmm_hits = {db_name: hmmsearch_all(db_name, refs[f"{db_name}_hmm"], checksums[db_key], combined, S["threads"],
                                   S["hmmsearch_evalue"], paths["work_dir"], V)
            for db_name, db_key in [("pfam", "Pfam-A.hmm.gz"), ("phrog", "phrogs_v4.hmm.gz")]}
for dataset in DATASETS:
    gs = gene_sets[dataset]
    annotations[dataset] = {
        "pfam": subset_hits(hmm_hits["pfam"], gs, V),
        "phrog": subset_hits(hmm_hits["phrog"], gs, V),
        "diamond": diamond_refs(refs["reference_proteins_fasta"], checksums["reference_proteins.faa"], gs,
                                S["threads"], S["diamond_mode"], S["diamond_evalue"], paths["work_dir"], dataset, V),
    }

# ---- rules and annotation ----------------------------------------------
from ece_replication import analysis as A  # noqa: E402
from ece_replication import plots as PL  # noqa: E402
from ece_replication.util import file_digest  # noqa: E402

load_rules = TaskGenerator(P.load_rules)
load_reference_table = TaskGenerator(P.load_reference_table)
load_pfam_meta = TaskGenerator(P.load_pfam_meta)
load_phrog_annot = TaskGenerator(P.load_phrog_annot)
build_features = TaskGenerator(P.build_features)
general_annotation = TaskGenerator(P.general_annotation)
origin_candidates = TaskGenerator(P.origin_candidates)
element_origins = TaskGenerator(P.element_origins)
combine_origins = TaskGenerator(P.combine_origins)
assign_evidence = TaskGenerator(P.assign_evidence)
gene_evidence_table = TaskGenerator(P.gene_evidence_table)
scan_origin_motifs = TaskGenerator(P.scan_origin_motifs)
enriched_cluster_origins = TaskGenerator(P.enriched_cluster_origins)

RULES_PATH = "config/replication_rules.yaml"
rules = load_rules(RULES_PATH, file_digest(RULES_PATH), V)
ref_table = load_reference_table(refs["reference_proteins_table"], file_digest(refs["reference_proteins_table"]), V)
pfam_meta = load_pfam_meta(refs["pfam_dat"], checksums["Pfam-A.hmm.dat.gz"], V)
phrog_annot = load_phrog_annot(refs["phrog_annot"], checksums["phrog_annot_v4.tsv"], V)

catalogues = {"ece": ece_catalogue, "validation": val_catalogue}
features, evidence_sets, gene_tables, origins_tables, motif_tables = {}, {}, {}, {}, {}
for dataset, (fasta, key) in DATASETS.items():
    hits = annotations[dataset]
    gs = gene_sets[dataset]
    feats = {th: build_features(gs, hits["pfam"], hits["phrog"], hits["diamond"], rules, ref_table, pfam_meta,
                                phrog_annot, th, dataset == "validation", V)
             for th in cfg["threshold_sets"]}
    features[dataset] = feats
    general = general_annotation(gs, hits["pfam"], hits["phrog"], rules, pfam_meta, phrog_annot, V)
    gene_tables[dataset] = gene_evidence_table(gs, feats["default"], general, rules, V)
    cands = origin_candidates(gs, feats["default"], V)
    seqs = {name: seq.upper() for name, _, seq in read_fasta(fasta)}
    per_elem = [element_origins(name, seq, gs, cands, cfg["origins"], V) for name, seq in seqs.items()]
    origins_tables[dataset] = combine_origins(per_elem, V)
    motif_tables[dataset] = scan_origin_motifs(fasta, checksums[key], rules, V)
    evidence_sets[dataset] = {th: assign_evidence(catalogues[dataset], feats[th], origins_tables[dataset],
                                                  motif_tables[dataset], rules, th, V)
                              for th in cfg["threshold_sets"]}

# ---- discovery, comparisons, validation, outputs ------------------------
protein_clusters = TaskGenerator(A.protein_clusters)
replication_network = TaskGenerator(A.replication_network)


@TaskGenerator
def element_groups(gene_set, fam, min_jaccard, version):
    from ece_replication import cluster
    return cluster.element_groups(gene_set["genes"], fam, min_jaccard)


@TaskGenerator
def discovery(gene_tab, feats, fam, evid, ncfg, version):
    from ece_replication import cluster
    focal = sorted(set(feats[feats.role == "initiator"].gene_id))
    nb = cluster.neighbourhoods(gene_tab, fam, focal, ncfg["flank_genes"])
    return {"neighbourhoods": nb,
            "unresolved_enriched": cluster.unresolved_enriched_clusters(gene_tab, fam, evid, ncfg["min_elements_for_recurrent"]),
            "conserved_neighbours": cluster.conserved_neighbours(nb, 2),
            "unknown_modules": cluster.recurrent_unknown_modules(gene_tab, fam, ncfg["min_elements_for_recurrent"])}


run_comparisons = TaskGenerator(A.run_comparisons)
shared_modules = TaskGenerator(A.shared_modules)
evaluate_validation = TaskGenerator(A.evaluate_validation)

fam = protein_clusters(gene_sets["ece"], S, paths["work_dir"], V)
groups = element_groups(gene_sets["ece"], fam, cfg["comparisons"]["group_jaccard"], V)
ece_ev = evidence_sets["ece"]["default"]
network = replication_network(gene_sets["ece"], features["ece"]["default"], S, paths["work_dir"], V)
disc = discovery(gene_tables["ece"], features["ece"]["default"], fam, evidence_sets["ece"]["default"]["elements"],
                 cfg["neighbourhood"], V)
comparison = run_comparisons(ece_catalogue, features["ece"]["default"], ece_ev["elements"], groups[0],
                             cfg["comparisons"], V)
novel_origins = enriched_cluster_origins(disc["unresolved_enriched"], gene_sets["ece"], fam, inp["fasta"],
                                         cfg["origins"], 0.1, V)
shared = shared_modules(ece_ev["elements"], features["ece"]["default"], fam, V)


@TaskGenerator
def replicon_similarity(query_fna, target_fna, checksum_q, checksum_t, work_dir, threads, version):
    """Nucleotide similarity of ECEs to characterised replicons (e.g. known origin regions)."""
    import gzip
    import shutil
    from ece_replication import search
    from ece_replication.util import project_path
    q = project_path(work_dir) / "ece" / "ece_sequences.fna"
    q.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(project_path(query_fna), "rb") as src, open(q, "wb") as dst:
        shutil.copyfileobj(src, dst)
    return search.mmseqs_nucleotide_search(q, target_fna, project_path(work_dir) / "ece" / "replicon_similarity.tsv",
                                           threads=threads)


replicon_hits = replicon_similarity(inp["fasta"], refs["validation_fasta"], checksums["SHD1_EC.fna.gz"],
                                    checksums["validation_elements.fna"], paths["work_dir"], S["threads"], V)
validation_eval = evaluate_validation([evidence_sets["validation"][th]["elements"] for th in cfg["threshold_sets"]],
                                      V)


@TaskGenerator
def write_results(ece_catalogue, gene_set, gene_tab, feats_by_th, ev_by_th, origins_df, fam, groups, comparison,
                  shared, disc, validation_eval, val_ev_default, val_feats, replicon_hits, motifs, novel_origins,
                  versions, checksums, out_dir, version):
    import pandas as pd
    from ece_replication import analysis as A
    ev = ev_by_th["default"]["elements"]
    gene_tab = gene_tab.assign(protein_cluster=gene_tab.gene_id.map(fam))
    sens = pd.concat([ev_by_th[t]["elements"][["element", "subgroup", "mechanism", "tier"]].assign(threshold_set=t)
                      for t in ev_by_th])
    sens_counts = sens.groupby(["threshold_set", "subgroup", "tier"]).size().rename("n").reset_index()
    wide = sens.pivot(index="element", columns="threshold_set", values="mechanism")
    tw = sens.pivot(index="element", columns="threshold_set", values="tier")
    stable = (wide.nunique(axis=1) == 1) & (tw.nunique(axis=1) == 1)
    sens_elem = wide.add_prefix("mechanism_").join(tw.add_prefix("tier_")).assign(stable=stable).reset_index()
    matrix = A.feature_matrix(ece_catalogue, feats_by_th["default"], ev)
    unknown = disc["unknown_modules"]
    sl = A.shortlist(ev, comparison, unknown, fam, gene_tab, disc["unresolved_enriched"])
    ev_out = ev.assign(related_group=ev.element.map(groups[0]))
    tables = {
        "catalogue": ece_catalogue,
        "gene_prediction_stats": gene_set["stats"],
        "genes": gene_tab,
        "gene_features": pd.concat([feats_by_th[t] for t in feats_by_th]),
        "element_evidence": ev_out,
        "initiator_calls": ev_by_th["default"]["initiators"],
        "origins": origins_df,
        "feature_matrix": matrix.reset_index().rename(columns={"index": "element"}),
        "sensitivity_counts": sens_counts,
        "sensitivity_by_element": sens_elem,
        "element_relatedness": groups[1],
        "comparisons/feature_omnibus": comparison["feature_omnibus"],
        "comparisons/feature_by_category": comparison["feature_by_category"],
        "comparisons/feature_omnibus_grouped": comparison["feature_omnibus_grouped"],
        "comparisons/feature_by_category_grouped": comparison["feature_by_category_grouped"],
        "comparisons/omnibus": comparison["omnibus"],
        "comparisons/length_summary": comparison["length_summary"],
        "comparisons/length_pairs": comparison["length_pairs"],
        "comparisons/length_spearman": comparison["spearman"],
        "comparisons/dual_marker_elements": A.dual_marker_table(comparison, feats_by_th["default"]),
        "discovery/initiator_clusters": shared,
        "discovery/initiator_neighbourhoods": disc["neighbourhoods"],
        "discovery/conserved_neighbours": disc["conserved_neighbours"],
        "discovery/recurrent_uncharacterised_modules": unknown,
        "discovery/unresolved_enriched_clusters": disc["unresolved_enriched"],
        "discovery/unresolved_enriched_cluster_origins": novel_origins,
        "origin_motifs": motifs,
        "shortlist": sl,
        "validation/evaluation": validation_eval,
        "discovery/nucleotide_similarity_to_characterised_replicons": replicon_hits,
        "validation/element_evidence": val_ev_default,
        "validation/gene_features": val_feats,
        "provenance/tool_versions": pd.DataFrame(sorted(versions.items()), columns=["tool", "version"]),
        "provenance/input_checksums": pd.DataFrame(sorted(checksums.items()), columns=["file", "sha256"]),
    }
    return A.write_tables(out_dir, tables)


@TaskGenerator
def make_figures(ece_catalogue, feats, ev_by_th, gene_tab, fam, shared, disc, network, rules, out_dir, version):
    import pandas as pd
    from ece_replication import analysis as A, plots as PL
    ev = ev_by_th["default"]["elements"]
    out = []
    out.append(PL.category_proportions(ev, f"{out_dir}/fig1_category_proportions.png"))
    m = A.feature_matrix(ece_catalogue, feats, ev)
    m = m[[c for c in m.columns if c.startswith(("role:", "pfam:"))]]
    from ece_replication.annotate import module_of
    colmod = {c: ("replication" if c.startswith("pfam:") else module_of(c.split(":", 1)[1], rules)) for c in m.columns}
    out.append(PL.feature_heatmap(m, ev, f"{out_dir}/fig2_feature_heatmap.png", colmod))
    out.append(PL.length_vs_mechanism(ev, f"{out_dir}/fig3_length_vs_mechanism.png"))
    ex = A.neighbourhood_examples(ev, shared, fam, gene_tab, disc["unknown_modules"], disc["unresolved_enriched"])
    out.append(PL.neighbourhoods(ex, gene_tab, ev, f"{out_dir}/fig4_neighbourhoods.png"))
    sub = ev.set_index("element").subgroup
    gt = gene_tab.set_index("gene_id")
    nodes = pd.DataFrame({"gene_id": network["ids"]})
    nodes["subgroup"] = nodes.gene_id.map(gt.element).map(sub)
    nodes["label"] = nodes.gene_id.map(gt.primary_family).fillna("")
    out.append(PL.replication_network(nodes, network["edges"], f"{out_dir}/fig5_replication_network.png"))
    sens = pd.concat([ev_by_th[t]["elements"].groupby("tier").size().rename("n").reset_index().assign(threshold_set=t)
                      for t in ev_by_th])
    out.append(PL.sensitivity(sens, f"{out_dir}/figS1_threshold_sensitivity.png"))
    return out


results_written = write_results(ece_catalogue, gene_sets["ece"], gene_tables["ece"], features["ece"],
                                evidence_sets["ece"], origins_tables["ece"], fam, groups, comparison, shared, disc,
                                validation_eval, evidence_sets["validation"]["default"]["elements"],
                                features["validation"]["default"], replicon_hits, motif_tables["ece"], novel_origins,
                                versions, checksums,
                                paths["results_dir"], V)
figures_written = make_figures(ece_catalogue, features["ece"]["default"], evidence_sets["ece"], gene_tables["ece"],
                               fam, shared, disc, network, rules, paths["figures_dir"], V)
