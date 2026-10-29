"""Matplotlib figures."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import networkx as nx  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import FancyArrow, Patch  # noqa: E402
from scipy.cluster.hierarchy import leaves_list, linkage  # noqa: E402
from scipy.spatial.distance import pdist  # noqa: E402

from .util import project_path  # noqa: E402

CAT_COLORS = {"plasmid": "#1f77b4", "virus": "#d62728", "uncategorized": "#7f7f7f",
              "uncategorized_dual_marker": "#9467bd"}
CAT_LABEL = {"plasmid": "Plasmid", "virus": "Virus", "uncategorized": "Uncategorized",
             "uncategorized_dual_marker": "Uncat., dual marker"}
MECH_COLORS = {"rolling_circle": "#e6550d", "theta": "#3182bd", "strand_displacement": "#31a354",
               "rolling_circle|theta": "#9e9ac8", "protein_primed": "#fdd0a2", "conflicting": "#000000",
               "unassigned": "#bdbdbd", "none_detected": "#f0f0f0"}
TIER_COLORS = {"strong": "#08519c", "tentative": "#6baed6", "conflicting": "#252525", "unresolved": "#f0f0f0"}
MODULE_COLORS = {"replication": "#e6550d", "maintenance": "#31a354", "transfer": "#756bb1",
                 "viral": "#d62728", "integration": "#8c6d31", "other": "#bdbdbd", "": "#ffffff"}


def mechanism_class(ev):
    return np.where(ev.tier == "conflicting", "conflicting", ev.mechanism)


def _save(fig, path):
    path = project_path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    return str(path)


def category_proportions(ev, path):
    ev = ev.assign(mclass=mechanism_class(ev))
    cats = [c for c in CAT_COLORS if c in set(ev.subgroup)]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    for ax, col, colors, title in [(axes[0], "mclass", MECH_COLORS, "Candidate mechanism"),
                                   (axes[1], "tier", TIER_COLORS, "Support tier")]:
        tab = pd.crosstab(ev.subgroup, ev[col]).reindex(cats).fillna(0)
        order = [k for k in colors if k in tab.columns] + [k for k in tab.columns if k not in colors]
        frac = tab[order].div(tab.sum(1), axis=0)
        bottom = np.zeros(len(cats))
        for k in order:
            ax.bar(range(len(cats)), frac[k].values, bottom=bottom, color=colors.get(k, "#cccccc"),
                   edgecolor="black", linewidth=0.4, label=k.replace("_", " "))
            bottom += frac[k].values
        ax.set_xticks(range(len(cats)))
        ax.set_xticklabels([f"{CAT_LABEL[c].replace(', ', chr(10))}\n(n={int(tab.loc[c].sum())})" for c in cats],
                           fontsize=8)
        ax.set_title(title)
        ax.legend(fontsize=7, frameon=False, bbox_to_anchor=(1.0, 1.0), loc="upper left")
    axes[0].set_ylabel("Fraction of elements")
    fig.tight_layout()
    return _save(fig, path)


def feature_heatmap(matrix, ev, path, column_modules=None):
    """matrix: elements x binary features (index = element); rows clustered (Jaccard, average linkage),
    columns grouped by module (``column_modules``: column -> module) and clustered within module."""
    m = matrix.loc[:, matrix.sum(0) > 0]
    X = m.values.astype(bool)
    order = leaves_list(linkage(np.nan_to_num(pdist(X, "jaccard"), nan=0.0), "average")) if len(m) > 2 \
        else np.arange(len(m))
    column_modules = column_modules or {c: "" for c in m.columns}
    cols, bounds = [], []
    for mod in ["replication", "maintenance", "transfer", "viral", "integration", "other", ""]:
        cs = [c for c in m.columns if column_modules.get(c, "other") == mod]
        if not cs:
            continue
        if len(cs) > 2:
            sub = m[cs].values.T.astype(bool)
            cs = [cs[i] for i in leaves_list(linkage(np.nan_to_num(pdist(sub, "jaccard"), nan=0.0), "average"))]
        cols.extend(cs)
        bounds.append((len(cols), mod))
    m = m.iloc[order][cols]
    sub = ev.set_index("element").loc[m.index]
    fig = plt.figure(figsize=(max(8, 0.28 * m.shape[1] + 3), max(8, 0.06 * len(m) + 2)))
    gs = fig.add_gridspec(1, 3, width_ratios=[0.6, 0.6, max(8, 0.28 * m.shape[1])], wspace=0.02)
    ax_c = fig.add_subplot(gs[0])
    ax_m = fig.add_subplot(gs[1])
    ax = fig.add_subplot(gs[2])
    ax.imshow(m.values, aspect="auto", cmap="Greys", interpolation="nearest", vmin=0, vmax=1)
    ax.set_xticks(range(m.shape[1]))
    ax.set_xticklabels(m.columns, rotation=90, fontsize=7)
    ax.set_yticks([])
    prev = 0
    for b, mod in bounds:
        ax.axvline(b - 0.5, color="#d62728", linewidth=0.8)
        ax.text((prev + b) / 2 - 0.5, -1.5, mod, ha="center", va="bottom", fontsize=8, color="#d62728")
        prev = b
    for a, vals, cmap in [(ax_c, sub.subgroup, CAT_COLORS), (ax_m, mechanism_class(sub.reset_index()), MECH_COLORS)]:
        rgb = np.array([[matplotlib.colors.to_rgb(cmap.get(v, "#cccccc"))] for v in vals])
        a.imshow(rgb, aspect="auto", interpolation="nearest")
        a.set_xticks([])
        a.set_yticks([])
    ax_c.set_title("cat.", fontsize=7)
    ax_m.set_title("mech.", fontsize=7)
    handles = [Patch(color=v, label=CAT_LABEL[k]) for k, v in CAT_COLORS.items()]
    handles += [Patch(color=v, label=k.replace("_", " ")) for k, v in MECH_COLORS.items()]
    ax.legend(handles=handles, fontsize=7, frameon=False, bbox_to_anchor=(1.01, 1), loc="upper left")
    return _save(fig, path)


def length_vs_mechanism(ev, path, seed=0):
    ev = ev.assign(mclass=mechanism_class(ev))
    order = [k for k in MECH_COLORS if k in set(ev.mclass)]
    rng = np.random.default_rng(seed)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    data = [np.log10(ev[ev.mclass == k].length) for k in order]
    ax.boxplot(data, positions=range(len(order)), widths=0.6, showfliers=False,
               medianprops=dict(color="black"), boxprops=dict(color="#555555"))
    for i, k in enumerate(order):
        sub = ev[ev.mclass == k]
        ax.scatter(i + rng.uniform(-0.22, 0.22, len(sub)), np.log10(sub.length), s=14,
                   c=[CAT_COLORS[c] for c in sub.subgroup], edgecolor="black", linewidth=0.3, zorder=3)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([f"{k.replace('_', ' ')}\n(n={(ev.mclass == k).sum()})" for k in order], fontsize=8)
    ticks = [1e3, 3e3, 1e4, 3e4, 1e5, 3e5]
    ax.set_yticks(np.log10(ticks))
    ax.set_yticklabels([f"{t / 1000:g} kb" for t in ticks])
    ax.set_ylabel("Element length")
    ax.legend(handles=[Patch(color=v, label=CAT_LABEL[k]) for k, v in CAT_COLORS.items()],
              fontsize=7, frameon=False, loc="upper left", bbox_to_anchor=(1, 1))
    fig.tight_layout()
    return _save(fig, path)


def _short(r):
    if r.get("primary_family"):
        return str(r["primary_family"])[:18]
    if r.get("best_pfam"):
        return str(r["best_pfam"])[:18]
    return ""


def neighbourhoods(examples, gene_tab, ev, path, flank_nt=6000):
    """examples: list of (title, focal_gene_id)."""
    gt = gene_tab.set_index("gene_id")
    cat = ev.set_index("element")
    fig, axes = plt.subplots(len(examples), 1, figsize=(12, 1.25 * len(examples) + 0.6), squeeze=False)
    for ax, (title, gid) in zip(axes[:, 0], examples):
        f = gt.loc[gid]
        e = f.element
        L = int(cat.loc[e, "length"])
        g = gene_tab[gene_tab.element == e].copy()
        centre = (f.start + (f.length_nt // 2)) % L
        # coordinates relative to focal gene centre on the circle, oriented so focal gene points right
        sign = 1 if f.strand == "+" else -1

        def rel(p):
            x = (p - centre + L / 2) % L - L / 2
            return x * sign

        for r in g.itertuples():
            s, en = rel(r.start), rel(r.start + r.length_nt - 1)
            lo, hi = min(s, en), max(s, en)
            if hi < -flank_nt or lo > flank_nt or hi - lo > L / 2:
                continue
            fwd = (r.strand == "+") == (sign == 1)
            x0, dx = (lo, hi - lo) if fwd else (hi, -(hi - lo))
            color = MODULE_COLORS.get(r.primary_module or "", "#ffffff")
            if not r.primary_module and not r.has_function:
                color = "#ffffff"
            elif not r.primary_module:
                color = "#d9d9d9"
            ax.add_patch(FancyArrow(x0, 0, dx, 0, width=0.5, head_width=0.5,
                                    head_length=min(250, abs(dx) * 0.4), length_includes_head=True,
                                    facecolor=color, edgecolor="black", linewidth=2 if r.gene_id == gid else 0.8))
            lab = _short(r._asdict())
            if lab:
                ax.text((lo + hi) / 2, 0.45, lab, ha="center", va="bottom", fontsize=6, rotation=20)
        ax.set_xlim(-flank_nt, flank_nt)
        ax.set_ylim(-0.6, 1.4)
        ax.axis("off")
        ax.set_title(f"{title}: {e} ({CAT_LABEL[cat.loc[e, 'subgroup']]}, {L / 1000:.1f} kb; "
                     f"{cat.loc[e, 'mechanism']} / {cat.loc[e, 'tier']})", fontsize=8, loc="left")
    handles = [Patch(facecolor=v, edgecolor="black", label=k) for k, v in MODULE_COLORS.items() if k]
    handles += [Patch(facecolor="#d9d9d9", edgecolor="black", label="other annotated"),
                Patch(facecolor="#ffffff", edgecolor="black", label="no annotation")]
    fig.legend(handles=handles, ncol=8, fontsize=7, frameon=False, loc="lower center")
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    return _save(fig, path)


def replication_network(nodes, edges, path, seed=7):
    """nodes: DataFrame gene_id, subgroup, label; edges: query, target, bits.

    Connected components are laid out separately and tiled (largest first);
    singletons are summarised in the title rather than drawn.
    """
    G = nx.Graph()
    for r in nodes.itertuples():
        G.add_node(r.gene_id, subgroup=r.subgroup, label=r.label)
    for r in edges.itertuples():
        if r.query != r.target and r.query in G and r.target in G:
            G.add_edge(r.query, r.target)
    comps = sorted((c for c in nx.connected_components(G) if len(c) > 1), key=lambda c: (-len(c), sorted(c)[0]))
    n_single = sum(1 for c in nx.connected_components(G) if len(c) == 1)
    ncol = 6
    nrow = max(1, int(np.ceil(len(comps) / ncol)))
    fig, axes = plt.subplots(nrow, ncol, figsize=(2.3 * ncol, 2.3 * nrow + 0.8), squeeze=False)
    for ax in axes.ravel():
        ax.axis("off")
    for ax, comp in zip(axes.ravel(), comps):
        H = G.subgraph(comp)
        pos = nx.spring_layout(H, seed=seed, iterations=200)
        nx.draw_networkx_edges(H, pos, ax=ax, alpha=0.3, width=0.6)
        nx.draw_networkx_nodes(H, pos, ax=ax, node_size=45,
                               node_color=[CAT_COLORS.get(H.nodes[n]["subgroup"], "#cccccc") for n in H],
                               edgecolors="black", linewidths=0.4)
        labs = pd.Series([H.nodes[n]["label"] for n in comp])
        labs = labs[labs != ""]
        title = labs.value_counts().index[0] if len(labs) else "unlabelled"
        ax.set_title(f"{title[:22]} (n={len(comp)})", fontsize=7)
        ax.margins(0.15)
    fig.legend(handles=[Patch(color=v, label=CAT_LABEL[k]) for k, v in CAT_COLORS.items()],
               fontsize=8, frameon=False, loc="lower center", ncol=4)
    fig.suptitle(f"Candidate initiator proteins: MMseqs2 similarity components (E<=1e-3); "
                 f"{len(comps)} components, {n_single} singletons not shown", fontsize=9)
    fig.tight_layout(rect=(0, 0.04, 1, 0.97))
    return _save(fig, path)


def sensitivity(tab, path):
    fig, ax = plt.subplots(figsize=(7, 3.5))
    piv = tab.pivot_table(index="threshold_set", columns="tier", values="n", fill_value=0)
    piv = piv.reindex([t for t in ["strict", "default", "relaxed"] if t in piv.index])
    left = np.zeros(len(piv))
    for t in [k for k in TIER_COLORS if k in piv.columns]:
        ax.barh(piv.index, piv[t], left=left, color=TIER_COLORS[t], edgecolor="black", label=t)
        left += piv[t].values
    ax.legend(fontsize=7, frameon=False, bbox_to_anchor=(1, 1), loc="upper left")
    ax.set_xlabel("Number of ECEs")
    fig.tight_layout()
    return _save(fig, path)
