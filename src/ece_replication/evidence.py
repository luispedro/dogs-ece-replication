"""Rule-based combination of gene-level features into element-level replication hypotheses.

Tiers
-----
strong      complete, high-specificity initiator for a single compatible
            mechanism set, corroborated by a characterised reference match or
            by origin-like DNA next to it, and with compatible context.
tentative   family-level / incomplete / uncorroborated evidence, or a
            replication module whose mechanism cannot be assigned.
conflicting initiators supporting incompatible mechanisms.
unresolved  no replication-associated gene detected. (This is NOT evidence
            that the element does not replicate: RNA-primed and host-dependent
            systems leave no recognisable initiator.)
"""
import pandas as pd

from .annotate import SPEC_RANK

INFORMATIVE = {"rolling_circle", "theta", "strand_displacement", "protein_primed"}


def _mechs(s):
    return {m for m in str(s).split(";") if m and m != "-"}


def gene_initiator_calls(feats, th):
    """Summarise initiator evidence per gene."""
    ini = feats[feats.role == "initiator"]
    calls = []
    for gid, g in ini.groupby("gene_id"):
        top = g.specificity.map(SPEC_RANK).max()
        topf = g[g.specificity.map(SPEC_RANK) == top]
        mech = set().union(*[_mechs(m) for m in topf.mechanisms]) & (INFORMATIVE | {"host_like", "unknown"})
        informative = mech & INFORMATIVE
        if informative:
            mech = informative
        # a mechanism-informative lower-tier feature can refine an 'unknown' top-tier call
        if not informative:
            lower = set().union(*[_mechs(m) for m in g.mechanisms]) & INFORMATIVE
            if lower:
                mech = lower
                top = min(top, SPEC_RANK["medium"])
        complete = bool(((g.source == "pfam") & (g.coverage >= th["initiator_full_coverage"])).any() or
                        ((g.source == "diamond") & (g.coverage >= th["diamond_system_coverage"])).any() or
                        ((g.source == "phrog") & (g.coverage >= th["initiator_full_coverage"])).any())
        dmd = g[g.source == "diamond"]
        ref_mech = set().union(*[_mechs(m) for m in dmd.mechanisms]) if len(dmd) else set()
        calls.append({
            "gene_id": gid, "element": g.element.iloc[0],
            "specificity": {v: k for k, v in SPEC_RANK.items()}[top],
            "mechanisms": ";".join(sorted(mech)) or "unknown",
            "complete": complete,
            "reference_support": bool(ref_mech & mech & INFORMATIVE),
            "reference_hits": " | ".join(dmd.detail),
            "transfer_lookalike": bool(g.transfer_lookalike.any()) and not bool(ref_mech & INFORMATIVE),
            "families": ";".join(dict.fromkeys(g.family_name)),
        })
    return pd.DataFrame(calls, columns=["gene_id", "element", "specificity", "mechanisms", "complete",
                                        "reference_support", "reference_hits", "transfer_lookalike", "families"])


def _origin_support(origins, gene_id, mechs, q, allowed=("dr", "ir", "at")):
    if origins is None or not len(origins):
        return False, ""
    o = origins[(origins.gene_id == gene_id) & (origins.control == False)]  # noqa: E712
    if not len(o):
        return False, ""
    notes, ok = [], False
    for r in o.itertuples():
        hits = []
        if r.dr_q <= q:
            hits.append(f"DR x{r.dr_copies}")
        if r.ir_q <= q:
            hits.append(f"IR stem {r.ir_stem}")
        if r.at_q <= q:
            hits.append(f"AT {r.at_max:.2f}")
        if hits:
            notes.append(f"{r.side}: " + ", ".join(hits))
        # iterons support theta; hairpins (dso/nick-site) support RCR; any feature supports others
        dr = "dr" in allowed and r.dr_q <= q
        ir = "ir" in allowed and r.ir_q <= q
        if "theta" in mechs and dr:
            ok = True
        if "rolling_circle" in mechs and ir:
            ok = True
        if mechs & {"strand_displacement", "protein_primed"} and (dr or ir):
            ok = True
    return ok, "; ".join(notes)


def assign_element(element, length, feats, calls, origins, comb, th_name, motifs=None):
    roles = set(feats.role)
    modules = {m: sorted(set(feats[feats.module == m].role)) for m in
               ["replication", "maintenance", "transfer", "viral", "integration"]}
    n_viral = feats[feats.module == "viral"].gene_id.nunique()
    viral_context = n_viral >= comb["viral_module_min_genes"]
    transfer_context = bool(roles & set(comb["transfer_context_roles"]))
    flags = []

    c = calls[calls.element == element]
    if len(c):
        mask = (c.mechanisms.map(lambda m: bool(_mechs(m) & INFORMATIVE)).astype(bool) &
                (c.specificity.map(SPEC_RANK) >= SPEC_RANK["medium"]))
        informative = c[mask.values]
    else:
        informative = c
    # ICE/MobT look-alikes are not counted as initiators when a transfer apparatus is present
    if transfer_context and informative.transfer_lookalike.any():
        flags.append("rep_trans_in_transfer_context")
        informative = informative[~informative.transfer_lookalike]

    motif_rows = motifs[motifs.element == element] if motifs is not None and len(motifs) else None
    motif_note = ""
    if motif_rows is not None and len(motif_rows):
        motif_note = "; ".join(f"{r.motif} at {r.position}{r.strand} ({r.mismatches} mm)" for r in motif_rows.itertuples())
    mech_sets = [_mechs(m) & INFORMATIVE for m in informative.mechanisms]
    support_genes = list(informative.gene_id)
    origin_notes = []
    if mech_sets:
        inter = set.intersection(*mech_sets)
        union = set.union(*mech_sets)
        if not inter:
            tier, mechanism = "conflicting", "|".join(sorted(union))
            flags.append("incompatible_initiators")
        else:
            mechanism = "|".join(sorted(inter))
            strong_genes = []
            for r in informative.itertuples():
                if r.specificity != "high" or not r.complete:
                    continue
                o_ok, o_note = _origin_support(origins, r.gene_id, inter, comb["origin_q_threshold"],
                                               tuple(comb.get("origin_support_features", ["dr", "ir", "at"])))
                if o_note:
                    origin_notes.append(f"{r.gene_id}: {o_note}")
                if r.reference_support or o_ok:
                    strong_genes.append(r.gene_id)
            context_ok = True
            # large dsDNA phages such as P2 replicate by RCR, so size only counts against
            # an RCR call when there is no viral module
            if inter == {"rolling_circle"} and length > comb["rcr_context_max_length"] and not viral_context:
                flags.append("rcr_initiator_on_large_element")
                context_ok = False
            if inter <= {"theta", "rolling_circle"} and len(inter) == 2 and not viral_context:
                flags.append("multimode_viral_initiator_without_viral_module")
            tier = "strong" if strong_genes and context_ok else "tentative"
            if len(informative) > 1:
                flags.append(f"{len(informative)}_compatible_initiators")
    else:
        # mechanism-informative accessories (e.g. IncQ primase, Rop, CopG-like regulators)
        if len(feats):
            acc_mask = ((feats.module == "replication") & (feats.role != "initiator") &
                        feats.mechanisms.map(lambda m: bool(_mechs(m) & INFORMATIVE)).astype(bool) &
                        (feats.specificity.map(SPEC_RANK) >= SPEC_RANK["medium"]))
            acc = feats[acc_mask.values]
        else:
            acc = feats
        rep = feats[feats.module == "replication"]
        if motif_rows is not None and len(motif_rows):
            acc = pd.concat([acc, pd.DataFrame({"gene_id": [f"origin:{r.motif}@{r.position}" for r in
                                                           motif_rows.itertuples()],
                                                "mechanisms": list(motif_rows.mechanisms)})], ignore_index=True)
            flags.append("characterised_origin_sequence")
        if len(acc):
            sets = [_mechs(m) & INFORMATIVE for m in acc.mechanisms]
            inter = set.intersection(*sets)
            mechanism = "|".join(sorted(inter)) if inter else "|".join(sorted(set.union(*sets)))
            tier = "tentative" if inter else "conflicting"
            support_genes = sorted(set(acc.gene_id))
            flags.append("no_initiator_mechanism_from_accessory")
        elif len(rep):
            tier, mechanism = "tentative", "unassigned"
            support_genes = sorted(set(rep.gene_id))
            if len(c):
                flags.append("initiator_of_unknown_mechanism")
            else:
                flags.append("replisome_components_only")
        else:
            tier, mechanism = "unresolved", "none_detected"
            if modules["transfer"]:
                flags.append("transfer_genes_without_replication_genes")
    if viral_context and mechanism in ("none_detected", "unassigned"):
        flags.append("viral_module_without_assigned_replication")
    return {
        "element": element, "threshold_set": th_name,
        "mechanism": mechanism, "tier": tier,
        "n_initiator_genes": int(len(c)), "n_informative_initiators": int(len(informative)),
        "support_genes": ";".join(support_genes),
        "initiator_families": " | ".join(dict.fromkeys(informative.families)) if len(informative) else
        " | ".join(dict.fromkeys(c.families)),
        "reference_support": " | ".join(x for x in informative.reference_hits if x) if len(informative) else "",
        "origin_support": " ; ".join(origin_notes),
        "origin_motifs": motif_note,
        **{f"{m}_roles": ";".join(v) for m, v in modules.items()},
        "n_viral_genes": int(n_viral),
        "viral_context": viral_context, "transfer_context": transfer_context,
        "flags": ";".join(flags),
    }


def assign_all(catalogue, feats, origins, comb, th, th_name, motifs=None):
    calls = gene_initiator_calls(feats, th)
    rows = []
    for r in catalogue.itertuples():
        f = feats[feats.element == r.element]
        rows.append(assign_element(r.element, r.length, f, calls, origins, comb, th_name, motifs))
    ev = pd.DataFrame(rows)
    keep = [c for c in ["element", "category", "subgroup", "note", "length", "gc", "expected_mechanism", "purpose"]
            if c in catalogue.columns]
    return catalogue[keep].merge(ev, on="element"), calls
