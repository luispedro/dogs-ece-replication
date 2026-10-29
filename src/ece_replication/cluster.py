"""Protein-family clustering, element relatedness and gene-neighbourhood discovery."""
from collections import Counter, defaultdict
from itertools import combinations

import networkx as nx
import pandas as pd


def family_map(clusters):
    """gene_id -> family id (representative gene id, prefixed)."""
    reps = {r: f"PC{i:05d}" for i, r in enumerate(sorted(clusters.representative.unique()))}
    return dict(zip(clusters.protein, clusters.representative.map(reps)))


def element_groups(genes, fam, min_jaccard):
    """Connected components of elements linked by protein-family Jaccard >= min_jaccard."""
    sets = genes.assign(fam=genes.gene_id.map(fam)).groupby("element").fam.apply(set).to_dict()
    G = nx.Graph()
    G.add_nodes_from(sets)
    inv = defaultdict(set)
    for e, s in sets.items():
        for f in s:
            inv[f].add(e)
    pairs = set()
    for es in inv.values():
        if len(es) > 1:
            pairs.update(combinations(sorted(es), 2))
    rows = []
    for a, b in pairs:
        j = len(sets[a] & sets[b]) / len(sets[a] | sets[b])
        rows.append((a, b, j))
        if j >= min_jaccard:
            G.add_edge(a, b, weight=j)
    groups = {}
    for i, comp in enumerate(sorted(nx.connected_components(G), key=lambda c: (-len(c), sorted(c)[0]))):
        for e in comp:
            groups[e] = f"G{i:03d}"
    return groups, pd.DataFrame(rows, columns=["element_a", "element_b", "jaccard"])


def neighbourhoods(gene_tab, fam, focal_ids, flank):
    """Rows: focal gene, neighbour offset (-flank..flank), neighbour gene/family/role (circular)."""
    rows = []
    by_elem = {e: g.sort_values("gene_index").reset_index(drop=True) for e, g in gene_tab.groupby("element")}
    pos = {gid: (e, i) for e, g in by_elem.items() for i, gid in enumerate(g.gene_id)}
    for gid in focal_ids:
        e, i = pos[gid]
        g = by_elem[e]
        n = len(g)
        focal_strand = g.loc[i, "strand"]
        seen = set()
        for off in range(-flank, flank + 1):
            j = (i + off) % n
            if j in seen and off != 0:
                continue
            seen.add(j)
            r = g.loc[j]
            rows.append({"focal_gene": gid, "element": e, "offset": off if focal_strand == "+" else -off,
                         "gene_id": r.gene_id, "family": fam.get(r.gene_id, ""), "strand": r.strand,
                         "same_strand": r.strand == focal_strand, "role": r.primary_role,
                         "module": r.primary_module, "best_pfam": r.best_pfam, "phrog_annot": r.phrog_annot,
                         "has_function": bool(r.has_function)})
    return pd.DataFrame(rows)


def recurrent_unknown_modules(gene_tab, fam, min_elements):
    """Adjacent pairs of protein families lacking functional annotation that recur in >= min_elements elements.

    Pairs are merged into modules (connected components of the recurrent-pair graph).
    """
    unknown_fam = gene_tab.assign(fam=gene_tab.gene_id.map(fam)).groupby("fam").has_function.any()
    unknown_fam = set(unknown_fam[~unknown_fam].index)
    pair_elems = defaultdict(set)
    for e, g in gene_tab.groupby("element"):
        g = g.sort_values("gene_index")
        fams = [fam[x] for x in g.gene_id]
        n = len(fams)
        if n < 2:
            continue
        for i in range(n if n > 2 else 1):
            a, b = fams[i], fams[(i + 1) % n]
            if a != b and a in unknown_fam and b in unknown_fam:
                pair_elems[tuple(sorted((a, b)))].add(e)
    G = nx.Graph()
    for (a, b), es in pair_elems.items():
        if len(es) >= min_elements:
            G.add_edge(a, b, elements=es)
    mods = []
    for i, comp in enumerate(sorted(nx.connected_components(G), key=lambda c: -len(c))):
        es = set()
        for a, b, d in G.subgraph(comp).edges(data=True):
            es |= d["elements"]
        mods.append({"module_id": f"UM{i:03d}", "families": ";".join(sorted(comp)), "n_families": len(comp),
                     "elements": ";".join(sorted(es)), "n_elements": len(es)})
    return pd.DataFrame(mods, columns=["module_id", "families", "n_families", "elements", "n_elements"])


def conserved_neighbours(nb, min_elements):
    """For each focal family, neighbour families present around it in >= min_elements distinct elements."""
    rows = []
    nb = nb.copy()
    focal_fam = nb[nb.offset == 0].set_index("focal_gene").family
    nb["focal_family"] = nb.focal_gene.map(focal_fam)
    for ff, g in nb[nb.offset != 0].groupby("focal_family"):
        n_focal_elems = g.element.nunique()
        cnt = g.groupby("family").element.nunique()
        for f, c in cnt.items():
            if c >= min_elements:
                sub = g[g.family == f]
                rows.append({"focal_family": ff, "neighbour_family": f, "n_elements": int(c),
                             "n_focal_elements": int(n_focal_elems),
                             "median_offset": float(sub.offset.median()),
                             "neighbour_role": Counter(sub.role).most_common(1)[0][0],
                             "neighbour_pfam": Counter(sub.best_pfam).most_common(1)[0][0],
                             "neighbour_has_function": bool(sub.has_function.any())})
    return pd.DataFrame(rows)


def unresolved_enriched_clusters(gene_tab, fam, evid, min_elements):
    """Protein clusters without functional annotation that are over-represented in ECEs lacking any
    recognisable replication gene: candidates for uncharacterised replication initiators."""
    from scipy.stats import fisher_exact
    unresolved = set(evid[evid.tier == "unresolved"].element)
    n_u, n_r = len(unresolved), evid.element.nunique() - len(unresolved)
    gt = gene_tab.assign(fam=gene_tab.gene_id.map(fam))
    func = gt.groupby("fam").has_function.any()
    rows = []
    for f, g in gt.groupby("fam"):
        if func[f]:
            continue
        es = set(g.element)
        k_u = len(es & unresolved)
        if k_u < min_elements:
            continue
        k_r = len(es - unresolved)
        p = fisher_exact([[k_u, n_u - k_u], [k_r, n_r - k_r]], alternative="greater")[1]
        rows.append({"cluster": f, "n_unresolved_elements": k_u, "n_other_elements": k_r,
                     "median_length_aa": float(g.length_aa.median()), "fisher_p_greater": p,
                     "elements": ";".join(sorted(es)), "example_gene": sorted(g.gene_id)[0]})
    out = pd.DataFrame(rows, columns=["cluster", "n_unresolved_elements", "n_other_elements", "median_length_aa",
                                      "fisher_p_greater", "elements", "example_gene"])
    if len(out):
        from scipy.stats import false_discovery_control
        out["q_bh"] = false_discovery_control(out.fisher_p_greater.values, method="bh")
    return out.sort_values("fisher_p_greater")
