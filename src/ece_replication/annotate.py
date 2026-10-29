"""Apply the family -> role/mechanism rules to search hits, producing gene-level evidence."""
import gzip

import numpy as np
import pandas as pd

from .util import project_path

SPEC_RANK = {"low": 1, "medium": 2, "high": 3}


def load_pfam_meta(dat_path):
    meta, cur = {}, {}
    with gzip.open(project_path(dat_path), "rt") as f:
        for line in f:
            if line.startswith("#=GF"):
                cur[line[5:7]] = line[10:].strip()
            elif line.startswith("//"):
                ga = cur.get("GA", "0; 0;").replace(";", " ").split()
                meta[cur["ID"]] = {"acc": cur["AC"].split(".")[0], "desc": cur.get("DE", ""),
                                   "ga_seq": float(ga[0]), "ga_dom": float(ga[1]), "type": cur.get("TP", "")}
                cur = {}
    return meta


def load_phrog_annot(path):
    df = pd.read_csv(project_path(path), sep="\t", dtype=str, keep_default_na=False)
    return {f"phrog_{r.phrog}": (r.annot, r.category) for r in df.itertuples()}


def module_of(role, rules):
    for mod, roles in rules["modules"].items():
        if role in roles:
            return mod
    return "other"


def pfam_pass(df, meta, th):
    """Boolean mask of (protein, model) rows passing the Pfam threshold set."""
    ga_seq = df.model.map(lambda m: meta.get(m, {}).get("ga_seq", np.inf))
    ga_dom = df.model.map(lambda m: meta.get(m, {}).get("ga_dom", np.inf))
    mode = th["pfam_score"]
    if mode == "gathering":
        ok = (df.full_score >= ga_seq) & (df.best_dom_score >= ga_dom)
    elif mode.startswith("gathering_plus_"):
        extra = float(mode.rsplit("_", 1)[1])
        ok = (df.full_score >= ga_seq + extra) & (df.best_dom_score >= ga_dom + extra)
    elif mode.startswith("evalue_"):
        ok = df.best_dom_ievalue <= float(mode.split("_", 1)[1])
    else:
        raise ValueError(mode)
    return ok & (df.hmm_coverage >= th["pfam_min_hmm_coverage"])


def build_features(genes, pfam, phrog, diamond, rules, refs, pfam_meta, phrog_annot, th, exclude_self=False):
    """Long table: one row per (gene, evidence source, family) that maps to a rule."""
    gene_elem = dict(zip(genes.gene_id, genes.element))
    feats = []
    pf_rules = rules["pfam"]
    if len(pfam):
        p = pfam[pfam_pass(pfam, pfam_meta, th)].copy()
        p["acc"] = p.model.map(lambda m: pfam_meta.get(m, {}).get("acc", ""))
        for r in p.itertuples():
            rule = pf_rules.get(r.acc)
            if rule is None:
                continue
            feats.append({"gene_id": r.protein, "source": "pfam", "family_id": r.acc, "family_name": r.model,
                          "role": rule["role"], "mechanisms": ";".join(rule.get("mechanisms", [])) or "-",
                          "specificity": rule["specificity"], "score": r.full_score,
                          "evalue": r.best_dom_ievalue, "coverage": r.hmm_coverage,
                          "identity": np.nan, "transfer_lookalike": bool(rule.get("transfer_lookalike", False)),
                          "detail": rule.get("family", "")})
    if len(phrog):
        q = phrog[(phrog.full_evalue <= th["phrog_evalue"]) & (phrog.hmm_coverage >= th["phrog_min_hmm_coverage"])]
        q = q.sort_values("full_score", ascending=False).drop_duplicates("protein")
        for r in q.itertuples():
            annot, category = phrog_annot.get(r.model, ("", ""))
            rule = rules["phrog"]["annot"].get(annot) or rules["phrog"]["category"].get(category)
            if rule is None:
                continue
            feats.append({"gene_id": r.protein, "source": "phrog", "family_id": r.model, "family_name": annot,
                          "role": rule["role"], "mechanisms": ";".join(rule.get("mechanisms", [])) or "-",
                          "specificity": rule["specificity"], "score": r.full_score, "evalue": r.full_evalue,
                          "coverage": r.hmm_coverage, "identity": np.nan, "transfer_lookalike": False,
                          "detail": f"{annot} [{category}]"})
    if len(diamond):
        d = diamond[(diamond.evalue <= th["diamond_evalue"]) & (diamond.pident >= th["diamond_min_identity"]) &
                    (diamond.qcov >= th["diamond_min_coverage"]) & (diamond.scov >= th["diamond_min_coverage"])]
        d = d.merge(refs, left_on="reference", right_on="accession")
        if exclude_self:
            d = d[d.source_element != d.protein.map(gene_elem)]
        d = d.sort_values("bitscore", ascending=False).drop_duplicates("protein")
        for r in d.itertuples():
            system = (r.pident >= th["diamond_system_identity"]) and (min(r.qcov, r.scov) >= th["diamond_system_coverage"])
            role = r.role
            if role == "relaxase_primase":
                role = "relaxase"
            feats.append({"gene_id": r.protein, "source": "diamond", "family_id": r.accession, "family_name": r.label,
                          "role": role, "mechanisms": r.mechanisms if r.mechanisms != "-" else "-",
                          "specificity": "high" if system else "medium", "score": r.bitscore, "evalue": r.evalue,
                          "coverage": min(r.qcov, r.scov), "identity": r.pident, "transfer_lookalike": False,
                          "detail": f"{r.label} ({'system' if system else 'family'}-level; {r.pident:.0f}% id)"})
    cols = ["gene_id", "source", "family_id", "family_name", "role", "mechanisms", "specificity", "score",
            "evalue", "coverage", "identity", "transfer_lookalike", "detail"]
    f = pd.DataFrame(feats, columns=cols)
    f["element"] = f.gene_id.map(gene_elem)
    f["module"] = f.role.map(lambda r: module_of(r, rules))
    return f


def best_general_annotation(genes, pfam, phrog, pfam_meta, phrog_annot, th):
    """Best Pfam (any family) and PHROG description per gene, for neighbourhood/novelty analysis."""
    out = pd.DataFrame({"gene_id": genes.gene_id})
    if len(pfam):
        p = pfam[pfam_pass(pfam, pfam_meta, th)].sort_values("full_score", ascending=False)
        names = p.groupby("protein").model.apply(lambda s: ";".join(s.head(4)))
        out["pfam_hits"] = out.gene_id.map(names).fillna("")
        best = p.drop_duplicates("protein").set_index("protein").model
        out["best_pfam"] = out.gene_id.map(best).fillna("")
    else:
        out["pfam_hits"] = out["best_pfam"] = ""
    if len(phrog):
        q = phrog[(phrog.full_evalue <= th["phrog_evalue"]) & (phrog.hmm_coverage >= th["phrog_min_hmm_coverage"])]
        q = q.sort_values("full_score", ascending=False).drop_duplicates("protein").set_index("protein")
        out["phrog"] = out.gene_id.map(q.model).fillna("")
        out["phrog_annot"] = out.phrog.map(lambda m: phrog_annot.get(m, ("", ""))[0] if m else "")
        out["phrog_category"] = out.phrog.map(lambda m: phrog_annot.get(m, ("", ""))[1] if m else "")
    else:
        out["phrog"] = out["phrog_annot"] = out["phrog_category"] = ""

    def characterised(r):
        pf = [x for x in r.pfam_hits.split(";") if x and not x.startswith(("DUF", "UPF"))]
        ph = r.phrog_annot and r.phrog_annot != "unknown function"
        return bool(pf) or bool(ph)

    out["has_function"] = [characterised(r) for r in out.itertuples()]
    return out


def gene_table(genes, features, general, rules):
    """One row per gene with its mapped evidence summarised."""
    g = genes.drop(columns=["nt_seq", "aa_seq"], errors="ignore").merge(general, on="gene_id", how="left")
    if len(features):
        f = features.copy()
        f["rank"] = f.specificity.map(SPEC_RANK) * 10 + f.source.map({"diamond": 3, "pfam": 2, "phrog": 1})
        f = f.sort_values("rank", ascending=False)
        best = f.drop_duplicates("gene_id").set_index("gene_id")
        agg = f.groupby("gene_id").agg(
            roles=("role", lambda s: ";".join(dict.fromkeys(s))),
            evidence=("detail", lambda s: " | ".join(dict.fromkeys(s))),
            sources=("source", lambda s: ";".join(sorted(set(s)))),
        )
        g = g.merge(agg, left_on="gene_id", right_index=True, how="left")
        g["primary_role"] = g.gene_id.map(best.role)
        g["primary_module"] = g.gene_id.map(best.module)
        g["primary_specificity"] = g.gene_id.map(best.specificity)
        g["primary_mechanisms"] = g.gene_id.map(best.mechanisms)
        g["primary_family"] = g.gene_id.map(best.family_name)
    for c in ["roles", "evidence", "sources", "primary_role", "primary_module", "primary_specificity",
              "primary_mechanisms", "primary_family"]:
        if c not in g:
            g[c] = ""
        g[c] = g[c].fillna("")
    return g
