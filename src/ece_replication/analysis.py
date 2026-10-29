"""ECE-level analyses: matrices, comparisons, discovery, validation, outputs."""
import numpy as np
import pandas as pd

from . import cluster, compare, plots, search
from .util import project_path, write_fasta

MAIN_ROLES = ["initiator", "helicase", "primase", "primase_helicase", "polymerase", "helicase_loader",
              "replisome_organiser", "ssb", "copy_control", "replication_accessory",
              "partition_atpase", "partition_cbp", "partition_other", "toxin", "antitoxin", "multimer_resolution",
              "relaxase", "mob_accessory", "coupling_protein", "t4ss",
              "terminase", "portal", "capsid", "tail", "prohead_protease", "head_tail", "viral_structural", "lysis",
              "integrase", "transposase"]


def protein_clusters(gene_set, scfg, work_dir, version):
    cl = search.mmseqs_cluster(gene_set["faa"], project_path(work_dir) / "ece" / "mmseqs_all",
                               min_seq_id=scfg["mmseqs_cluster_min_seq_id"],
                               coverage=scfg["mmseqs_cluster_coverage"], threads=scfg["threads"])
    return cluster.family_map(cl)


def replication_network(gene_set, feats, scfg, work_dir, version):
    ids = sorted(set(feats[feats.role == "initiator"].gene_id))
    g = gene_set["genes"].set_index("gene_id").loc[ids]
    faa = write_fasta(zip(g.index, g.aa_seq), project_path(work_dir) / "ece" / "initiators.faa")
    hits = search.mmseqs_allvsall(faa, project_path(work_dir) / "ece" / "initiators_allvsall.tsv",
                                  evalue=scfg["mmseqs_network_evalue"], threads=scfg["threads"])
    return {"ids": ids, "edges": hits}


def feature_matrix(catalogue, feats, evid):
    """Binary ECE x feature matrix: roles, mechanism-informative families, mechanism and tier indicators."""
    m = pd.DataFrame(index=catalogue.element)
    for role in MAIN_ROLES:
        es = set(feats[feats.role == role].element)
        m[f"role:{role}"] = m.index.isin(es)
    fam = feats[(feats.module == "replication") & (feats.source == "pfam")]
    for f, g in fam.groupby("family_name"):
        m[f"pfam:{f}"] = m.index.isin(set(g.element))
    ev = evid.set_index("element")
    mclass = pd.Series(plots.mechanism_class(evid), index=evid.element)
    for k in sorted(set(mclass)):
        m[f"mechanism:{k}"] = (mclass == k).reindex(m.index).values
    for t in ["strong", "tentative", "conflicting", "unresolved"]:
        m[f"tier:{t}"] = (ev.tier == t).reindex(m.index).values
    return m.astype(bool)


def module_presence(feats, catalogue):
    out = pd.DataFrame(index=catalogue.element)
    defs = {
        "has_initiator": feats.role == "initiator",
        "has_replisome_component": feats.role.isin(["helicase", "primase", "primase_helicase", "polymerase",
                                                    "helicase_loader", "replisome_organiser", "ssb"]),
        "has_partition": feats.role.isin(["partition_atpase", "partition_cbp", "partition_other"]),
        "has_toxin_antitoxin": feats.role.isin(["toxin", "antitoxin"]),
        "has_relaxase": feats.role == "relaxase",
        "has_conjugation_apparatus": feats.role.isin(["coupling_protein", "t4ss"]),
        "has_viral_structural": feats.module == "viral",
        "has_integrase": feats.role == "integrase",
        "has_transposase": feats.role == "transposase",
    }
    for k, mask in defs.items():
        out[k] = out.index.isin(set(feats[mask].element))
    return out


def run_comparisons(catalogue, feats, evid, groups, ccfg, version):
    ev = evid.copy()
    ev["mclass"] = plots.mechanism_class(ev)
    mods = module_presence(feats, catalogue)
    df = ev.merge(mods, left_on="element", right_index=True)
    df["group"] = df.element.map(groups)
    for k in ["rolling_circle", "theta", "strand_displacement", "rolling_circle|theta", "unassigned",
              "none_detected", "conflicting"]:
        df[f"mech_{k}"] = df.mclass == k
    df["mech_any_assigned"] = df.mclass.isin(["rolling_circle", "theta", "strand_displacement",
                                              "rolling_circle|theta", "protein_primed"])
    df["tier_strong"] = df.tier == "strong"
    df["tier_unresolved"] = df.tier == "unresolved"
    mech_feats = [c for c in df.columns if c.startswith("mech_") and df[c].any()]
    feats_bin = mech_feats + ["tier_strong", "tier_unresolved"] + list(mods.columns)
    omni, per_cat = compare.binary_feature_tests(df, feats_bin, "subgroup", ccfg["n_permutations"], ccfg["seed"])
    omni_mech = compare.permutation_chi2(df.subgroup.values, df.mclass.values, ccfg["n_permutations"], ccfg["seed"])
    omni_tier = compare.permutation_chi2(df.subgroup.values, df.tier.values, ccfg["n_permutations"], ccfg["seed"] + 1)
    # one representative (lowest identifier) per group of related ECEs
    rep = df.sort_values("element").drop_duplicates("group")
    omni_g, per_cat_g = compare.binary_feature_tests(rep, feats_bin, "subgroup", ccfg["n_permutations"],
                                                     ccfg["seed"])
    omni_mech_g = compare.permutation_chi2(rep.subgroup.values, rep.mclass.values, ccfg["n_permutations"],
                                           ccfg["seed"])
    len_summary, len_pairs, kw = compare.length_tests(df, "mclass")
    df["log10_length"] = np.log10(df.length)
    df["n_replication_roles"] = df.replication_roles.map(lambda s: len([x for x in s.split(";") if x]))
    df["n_maintenance_roles"] = df.maintenance_roles.map(lambda s: len([x for x in s.split(";") if x]))
    df["n_transfer_roles"] = df.transfer_roles.map(lambda s: len([x for x in s.split(";") if x]))
    df["n_modules"] = df[["has_initiator", "has_partition", "has_toxin_antitoxin", "has_relaxase",
                          "has_conjugation_apparatus", "has_viral_structural", "has_integrase"]].sum(1)
    spear = compare.spearman_table(df, "log10_length", ["n_replication_roles", "n_maintenance_roles",
                                                        "n_transfer_roles", "n_modules"])
    omnibus = pd.DataFrame([
        {"test": "subgroup x mechanism class (all ECEs)", "chi2": omni_mech[0], "p_perm": omni_mech[1],
         "cramers_v": omni_mech[2], "n": len(df)},
        {"test": "subgroup x support tier (all ECEs)", "chi2": omni_tier[0], "p_perm": omni_tier[1],
         "cramers_v": omni_tier[2], "n": len(df)},
        {"test": "subgroup x mechanism class (one per related group)", "chi2": omni_mech_g[0],
         "p_perm": omni_mech_g[1], "cramers_v": omni_mech_g[2], "n": len(rep)},
        {"test": "log10 length ~ mechanism class (Kruskal-Wallis)", "chi2": np.nan, "p_perm": kw,
         "cramers_v": np.nan, "n": len(df)},
    ])
    return {"table": df, "feature_omnibus": omni, "feature_by_category": per_cat,
            "feature_omnibus_grouped": omni_g, "feature_by_category_grouped": per_cat_g,
            "omnibus": omnibus, "length_summary": len_summary, "length_pairs": len_pairs,
            "spearman": spear, "n_groups": int(rep.shape[0])}


def dual_marker_table(comp, feats):
    df = comp["table"]
    d = df[df.subgroup == "uncategorized_dual_marker"].copy()
    d["plasmid_like_modules"] = d[["has_partition", "has_relaxase", "has_conjugation_apparatus"]].any(axis=1)
    d["viral_modules"] = d.has_viral_structural
    cols = ["element", "length", "mechanism", "tier", "initiator_families", "replication_roles",
            "maintenance_roles", "transfer_roles", "viral_roles", "integration_roles", "plasmid_like_modules",
            "viral_modules", "flags"]
    return d[cols]


def shared_modules(evid, feats, fam, version):
    """Initiator protein clusters and the ECE categories that carry them."""
    ini = feats[feats.role == "initiator"][["gene_id", "element", "family_name", "source"]].copy()
    ini["cluster"] = ini.gene_id.map(fam)
    sub = evid.set_index("element").subgroup
    ini["subgroup"] = ini.element.map(sub)
    rows = []
    for c, g in ini.groupby("cluster"):
        es = g.drop_duplicates("element")
        cnt = es.subgroup.value_counts()
        rows.append({"cluster": c, "n_elements": len(es), "elements": ";".join(sorted(es.element)),
                     **{f"n_{k}": int(cnt.get(k, 0)) for k in plots.CAT_COLORS},
                     "families": ";".join(sorted(set(g[g.source == "pfam"].family_name))) or
                     ";".join(sorted(set(g.family_name)))[:120]})
    out = pd.DataFrame(rows).sort_values("n_elements", ascending=False)
    labelled = out.n_plasmid + out.n_virus
    unc = out.n_uncategorized + out.n_uncategorized_dual_marker
    out["shared_uncat_labelled"] = (labelled > 0) & (unc > 0)
    return out


def evaluate_validation(evs, version):
    rows = []
    for ev in evs:
        for r in ev.itertuples():
            exp = set(r.expected_mechanism.split(";"))
            got = set(r.mechanism.split("|")) if r.tier != "conflicting" else {"conflicting"}
            if exp == {"none"}:
                outcome = "correct" if r.mechanism in ("none_detected", "unassigned") else "false_positive"
            elif exp == {"unknown"}:
                outcome = "exploratory"
            elif got == exp:
                outcome = "correct"
            elif r.tier == "conflicting" and exp <= set(r.mechanism.split("|")):
                outcome = "flagged_conflict"
            elif r.mechanism in ("none_detected", "unassigned"):
                outcome = "missed"
            elif got & exp:
                outcome = "partial"
            else:
                outcome = "wrong"
            rows.append({"element": r.element, "threshold_set": r.threshold_set, "expected": r.expected_mechanism,
                         "assigned": r.mechanism, "tier": r.tier, "outcome": outcome,
                         "initiator_families": r.initiator_families, "reference_support": r.reference_support,
                         "flags": r.flags, "purpose": r.purpose})
    return pd.DataFrame(rows)


HIGH_PRIORITY = ("incompatible initiators", "dual plasmid/virus markers", "phage-plasmid candidate",
                 "candidate novel initiator", "recurrent uncharacterised module", "MobT vs RCR",
                 "RCR-type initiator on a large element")


def shortlist(evid, comp, unknown_modules, fam, gene_tab, enriched=None, q_max=0.1):
    """Unresolved / unusual elements for deeper investigation, with a priority tier."""
    df = comp["table"]
    rows = []
    novel = {}
    if enriched is not None and len(enriched):
        for r in enriched[enriched.q_bh <= q_max].itertuples():
            for e in r.elements.split(";"):
                novel.setdefault(e, []).append(r.cluster)
    um_elems = {}
    for r in unknown_modules.itertuples():
        for e in r.elements.split(";"):
            um_elems.setdefault(e, []).append(r.module_id)
    for r in df.itertuples():
        reasons = []
        if r.tier == "unresolved":
            reasons.append("no replication-associated gene detected")
        if r.tier == "conflicting":
            reasons.append("incompatible initiators")
        if "rcr_initiator_on_large_element" in r.flags:
            reasons.append("RCR-type initiator on a large element")
        if "viral_module_without_assigned_replication" in r.flags:
            reasons.append("viral structural module without assigned replication mechanism")
        if r.subgroup == "uncategorized_dual_marker":
            reasons.append("dual plasmid/virus markers")
        if r.subgroup == "virus" and r.has_partition:
            reasons.append("virus carrying partition genes (phage-plasmid candidate)")
        if "rep_trans_in_transfer_context" in r.flags:
            reasons.append("Rep_trans-like protein in conjugative context (MobT vs RCR)")
        if "replisome_components_only" in r.flags:
            reasons.append("replisome components without an initiator")
        if r.tier == "unresolved" and r.element in um_elems:
            reasons.append("carries recurrent uncharacterised module " + ",".join(um_elems[r.element]))
        if r.element in novel:
            reasons.append("carries unannotated cluster enriched in unresolved ECEs (candidate novel initiator) " +
                           ",".join(novel[r.element]))
        if reasons:
            rows.append({"element": r.element, "subgroup": r.subgroup, "length": r.length,
                         "mechanism": r.mechanism, "tier": r.tier, "reasons": "; ".join(reasons),
                         "initiator_families": r.initiator_families, "flags": r.flags})
    out = pd.DataFrame(rows)
    out["n_reasons"] = out.reasons.str.count("; ") + 1
    out["priority"] = np.where(out.reasons.map(lambda s: any(k in s for k in HIGH_PRIORITY)), "high", "routine")
    return out.sort_values(["priority", "n_reasons", "subgroup", "element"], ascending=[True, False, True, True])


def neighbourhood_examples(evid, shared, fam, gene_tab, unknown_modules, enriched=None, n_clusters=5, n_unknown=1,
                           n_enriched=2):
    tier_rank = {"strong": 0, "tentative": 1, "conflicting": 2, "unresolved": 3}
    ev = evid.set_index("element")
    gt = gene_tab.assign(cluster=gene_tab.gene_id.map(fam))
    ex = []
    for r in shared.head(n_clusters).itertuples():
        members = gt[(gt.cluster == r.cluster) & gt.element.isin(r.elements.split(";"))].copy()
        members["rank"] = members.element.map(lambda e: tier_rank[ev.loc[e, "tier"]])
        best = members.sort_values(["rank", "element"]).iloc[0]
        ex.append((f"Initiator cluster {r.cluster} ({r.families or 'unclassified'}; n={r.n_elements})",
                   best.gene_id))
    for r in unknown_modules.head(n_unknown).itertuples():
        fams = set(r.families.split(";"))
        e = sorted(r.elements.split(";"))[0]
        g = gt[(gt.element == e) & gt.cluster.isin(fams)].sort_values("gene_index").iloc[0]
        ex.append((f"Recurrent uncharacterised module {r.module_id} (n={r.n_elements} ECEs)", g.gene_id))
    if enriched is not None:
        for r in enriched.head(n_enriched).itertuples():
            ex.append((f"Unannotated cluster {r.cluster} enriched in unresolved ECEs "
                       f"({r.n_unresolved_elements} unresolved, {r.n_other_elements} other)", r.example_gene))
    return ex


def write_tables(out_dir, tables):
    out = project_path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    written = []
    for name, df in tables.items():
        p = out / f"{name}.tsv"
        p.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(p, sep="\t", index=False, float_format="%.4g")
        written.append(str(p))
    return written
