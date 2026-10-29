"""Statistical comparisons between ECE categories (small-count friendly)."""
import numpy as np
import pandas as pd
from scipy import stats

CATEGORY_ORDER = ["plasmid", "virus", "uncategorized", "uncategorized_dual_marker"]


def wilson(k, n, z=1.959964):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def odds_ratio_ci(a, b, c, d, z=1.959964):
    """Haldane-Anscombe corrected OR and Wald CI for [[a, b], [c, d]]."""
    a, b, c, d = (x + 0.5 for x in (a, b, c, d))
    lor = np.log(a * d / (b * c))
    se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    return float(np.exp(lor)), float(np.exp(lor - z * se)), float(np.exp(lor + z * se))


def _chi2_stat(labels, values, ncat, nval):
    t = np.zeros((ncat, nval))
    np.add.at(t, (labels, values), 1)
    exp = t.sum(1, keepdims=True) * t.sum(0, keepdims=True) / t.sum()
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.nansum(np.where(exp > 0, (t - exp) ** 2 / exp, 0.0))


def permutation_chi2(labels, values, n_perm, seed):
    """Omnibus association between two categorical vectors; P from label permutations."""
    lab_codes, lab_u = pd.factorize(pd.Series(labels))
    val_codes, val_u = pd.factorize(pd.Series(values))
    if len(val_u) < 2 or len(lab_u) < 2:
        return np.nan, 1.0, 0.0
    obs = _chi2_stat(lab_codes, val_codes, len(lab_u), len(val_u))
    rng = np.random.default_rng(seed)
    ge = sum(_chi2_stat(rng.permutation(lab_codes), val_codes, len(lab_u), len(val_u)) >= obs - 1e-12
             for _ in range(n_perm))
    n = len(labels)
    v = np.sqrt(obs / (n * (min(len(lab_u), len(val_u)) - 1)))
    return float(obs), float((ge + 1) / (n_perm + 1)), float(v)


def bh(p):
    p = np.asarray(p, dtype=float)
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    if ok.any():
        out[ok] = stats.false_discovery_control(p[ok], method="bh")
    return out


def binary_feature_tests(df, features, group_col, n_perm, seed, categories=CATEGORY_ORDER):
    """Per binary feature: per-category proportion (Wilson CI), category-vs-rest OR + Fisher P,
    omnibus permutation chi-square across categories; BH q-values over all tests."""
    cats = [c for c in categories if c in set(df[group_col])]
    per_cat, omni = [], []
    for i, f in enumerate(features):
        x = df[f].astype(bool).values
        chi, p, v = permutation_chi2(df[group_col].values, x, n_perm, seed + i)
        omni.append({"feature": f, "chi2": chi, "p_perm": p, "cramers_v": v, "n": len(x),
                     "n_with_feature": int(x.sum())})
        for c in cats:
            m = (df[group_col] == c).values
            a, b = int((x & m).sum()), int((~x & m).sum())
            cc, d = int((x & ~m).sum()), int((~x & ~m).sum())
            lo, hi = wilson(a, a + b)
            orr, olo, ohi = odds_ratio_ci(a, b, cc, d)
            per_cat.append({"feature": f, "category": c, "n_category": a + b, "k": a,
                            "proportion": a / (a + b) if a + b else np.nan, "ci_low": lo, "ci_high": hi,
                            "rest_proportion": cc / (cc + d) if cc + d else np.nan,
                            "odds_ratio": orr, "or_ci_low": olo, "or_ci_high": ohi,
                            "fisher_p": stats.fisher_exact([[a, b], [cc, d]])[1]})
    omni = pd.DataFrame(omni)
    per_cat = pd.DataFrame(per_cat)
    if len(omni):
        omni["q_bh"] = bh(omni.p_perm)
    if len(per_cat):
        per_cat["q_bh"] = bh(per_cat.fisher_p)
    return omni, per_cat


def length_tests(df, class_col, min_n=3):
    """Kruskal-Wallis of log10 length across classes with >= min_n members, plus pairwise MWU (BH)."""
    d = df.copy()
    d["log10_length"] = np.log10(d.length)
    keep = d[class_col].value_counts()
    keep = keep[keep >= min_n].index
    d = d[d[class_col].isin(keep)]
    groups = {k: g.log10_length.values for k, g in d.groupby(class_col)}
    summary = d.groupby(class_col).length.describe()[["count", "min", "25%", "50%", "75%", "max"]].reset_index()
    if len(groups) < 2:
        return summary, pd.DataFrame(), np.nan
    kw = stats.kruskal(*groups.values()).pvalue
    pairs = []
    ks = sorted(groups)
    for i in range(len(ks)):
        for j in range(i + 1, len(ks)):
            a, b = groups[ks[i]], groups[ks[j]]
            u = stats.mannwhitneyu(a, b, alternative="two-sided")
            pairs.append({"class_a": ks[i], "class_b": ks[j], "n_a": len(a), "n_b": len(b),
                          "median_a": float(10 ** np.median(a)), "median_b": float(10 ** np.median(b)),
                          "rank_biserial": float(2 * u.statistic / (len(a) * len(b)) - 1),  # >0: a larger
                          "p": float(u.pvalue)})
    pairs = pd.DataFrame(pairs)
    pairs["q_bh"] = bh(pairs.p)
    return summary, pairs, float(kw)


def spearman_table(df, x, ys):
    rows = []
    for y in ys:
        r = stats.spearmanr(df[x], df[y])
        rows.append({"x": x, "y": y, "rho": float(r.statistic), "p": float(r.pvalue), "n": len(df)})
    out = pd.DataFrame(rows)
    out["q_bh"] = bh(out.p)
    return out
