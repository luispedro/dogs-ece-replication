"""Pipeline step functions (pure functions of their arguments; wrapped as Jug tasks in jugfile.py).

Every step takes ``version`` (the analysis version string) so that bumping
the version in config/analysis.yaml invalidates cached results.  Steps that
consume a file also take the file's checksum for the same reason.
"""
import hashlib
from pathlib import Path

import pandas as pd

from . import catalogue, genes, search
from .fetch_refs import verify
from .util import project_path, read_fasta, write_fasta


def verify_inputs(version):
    return verify()


def load_sequences(fasta_path, checksum, version):
    return {name: seq.upper() for name, _, seq in read_fasta(fasta_path)}


def validate_catalogue(fasta_path, props_path, expected, dual_note, checksums, version):
    cat, _ = catalogue.build_catalogue(fasta_path, props_path, expected, dual_note)
    return cat


def validation_catalogue(fasta_path, table_path, checksum, version):
    meta = pd.read_csv(project_path(table_path), sep="\t", dtype=str)
    rows = []
    for name, desc, seq in read_fasta(fasta_path):
        st = catalogue.sequence_stats(seq)
        m = meta[meta.element == name].iloc[0].to_dict()
        rows.append({"element": name, "accession": desc, "category": m["expected_class"],
                     "subgroup": m["expected_class"], "expected_mechanism": m["expected_mechanism"],
                     "purpose": m["purpose"], **st})
    return pd.DataFrame(rows)


def predict_genes(element, seq, gp_cfg, version):
    return genes.predict_circular(element, seq, rotations=tuple(gp_cfg["rotations"]),
                                  min_gene=gp_cfg["min_gene_length_nt"])


def collect_genes(dataset, results, work_dir, version):
    """Merge per-element predictions; write proteins/CDS FASTA. Returns (genes df, stats df, faa path, sha)."""
    gl, sl = [], []
    for g, s in results:
        gl.extend(g)
        sl.append(s)
    gdf = pd.DataFrame(gl)
    sdf = pd.DataFrame(sl)
    out = project_path(work_dir) / dataset
    out.mkdir(parents=True, exist_ok=True)
    faa = write_fasta(zip(gdf.gene_id, gdf.aa_seq), out / "proteins.faa")
    write_fasta(zip(gdf.gene_id, gdf.nt_seq), out / "cds.fna")
    sha = hashlib.sha256(Path(faa).read_bytes()).hexdigest()
    return {"genes": gdf, "stats": sdf, "faa": faa, "faa_sha256": sha}


def combine_proteins(gene_sets, work_dir, version):
    """Single FASTA with the proteins of all datasets (identifiers are unique across datasets)."""
    gdf = pd.concat([gs["genes"] for gs in gene_sets], ignore_index=True)
    if gdf.gene_id.duplicated().any():
        raise ValueError("duplicate gene identifiers across datasets")
    faa = write_fasta(zip(gdf.gene_id, gdf.aa_seq), project_path(work_dir) / "all_proteins.faa")
    return {"faa": faa, "faa_sha256": hashlib.sha256(Path(faa).read_bytes()).hexdigest()}


def hmmsearch_all(db_name, hmm_path, db_sha, combined, cpu, evalue, work_dir, version):
    out = project_path(work_dir) / "hmm" / f"{db_name}.domtbl"
    search.hmmsearch(hmm_path, combined["faa"], out, cpu=cpu, evalue=evalue)
    return search.summarise_domtbl(search.parse_domtbl(out))


def subset_hits(hits, gene_set, version):
    return hits[hits.protein.isin(set(gene_set["genes"].gene_id))].reset_index(drop=True)


def diamond_refs(ref_faa, ref_sha, gene_set, threads, mode, evalue, work_dir, dataset, version):
    out = project_path(work_dir) / dataset / "diamond_refs.tsv"
    search.diamond_blastp(ref_faa, gene_set["faa"], out, threads=threads, mode=mode, evalue=evalue)
    return search.parse_diamond(out)


def tool_versions(version):
    return {
        "hmmsearch": search.run(["hmmsearch", "-h"]).stdout.splitlines()[1].strip("# "),
        "diamond": search.tool_version(["diamond", "version"]),
        "mmseqs": search.tool_version(["mmseqs", "version"]),
        "pyrodigal": __import__("pyrodigal").__version__,
        "jug": __import__("jug").__version__,
        "pandas": pd.__version__,
    }


# ---------------------------------------------------------------------------
# annotation and evidence
# ---------------------------------------------------------------------------
from . import annotate, evidence, origins as orig  # noqa: E402


def load_rules(path, checksum, version):
    from .util import load_yaml
    return load_yaml(path)


def load_reference_table(path, checksum, version):
    return pd.read_csv(project_path(path), sep="\t", dtype=str)


def load_pfam_meta(path, checksum, version):
    return annotate.load_pfam_meta(path)


def load_phrog_annot(path, checksum, version):
    return annotate.load_phrog_annot(path)


def build_features(gene_set, pfam, phrog, diamond, rules, refs, pfam_meta, phrog_annot, th_name, exclude_self,
                   version):
    th = rules["thresholds"][th_name]
    f = annotate.build_features(gene_set["genes"], pfam, phrog, diamond, rules, refs, pfam_meta, phrog_annot, th,
                                exclude_self=exclude_self)
    f["threshold_set"] = th_name
    return f


def general_annotation(gene_set, pfam, phrog, rules, pfam_meta, phrog_annot, version):
    return annotate.best_general_annotation(gene_set["genes"], pfam, phrog, pfam_meta, phrog_annot,
                                            rules["thresholds"]["default"])


def origin_candidates(gene_set, feats, version):
    """Per element: replication-module genes to test (initiators if any, else other replication genes)."""
    rep = feats[feats.module == "replication"]
    out = {e: [] for e in gene_set["genes"].element.unique()}
    for e, g in rep.groupby("element"):
        ini = sorted(set(g[g.role == "initiator"].gene_id))
        out[e] = ini if ini else sorted(set(g.gene_id))
    return out


def element_origins(element, seq, gene_set, candidates_by_element, ocfg, version):
    import zlib
    import numpy as np
    candidates = candidates_by_element.get(element, [])
    gdf = gene_set["genes"]
    gl = gdf[gdf.element == element].to_dict("records")
    L = len(seq)
    if not candidates or not gl:
        return pd.DataFrame()
    inter = genes.intergenic_regions(gl, L)
    seed = ocfg["seed"] + zlib.crc32(element.encode())
    rng = np.random.default_rng(seed)
    # controls: genes at least two positions away from every candidate, so that
    # their flanking intergenic regions are not shared with a candidate
    cand_idx = [g["gene_index"] for g in gl if g["gene_id"] in set(candidates)]
    n = len(gl)

    def far(i):
        return all(min(abs(i - j), n - abs(i - j)) >= 2 for j in cand_idx)

    others = [g["gene_id"] for g in gl if g["gene_id"] not in set(candidates) and far(g["gene_index"])]
    controls = list(rng.choice(others, size=min(2, len(others)), replace=False)) if others else []
    rows = []
    for gid, is_control in [(c, False) for c in candidates] + [(c, True) for c in controls]:
        g = next(x for x in gl if x["gene_id"] == gid)
        for side, s, e in orig.flanking_regions(g, inter, L, ocfg["flank_nt"]):
            region = orig.circular_slice(seq, s, e)
            if len(region) < ocfg["min_region_nt"]:
                continue
            st = orig.test_region(region, seq, ocfg, seed=seed + s)
            rows.append({"element": element, "gene_id": gid, "control": is_control, "side": side,
                         "start": s, "end": e, "region_length": len(region), **st})
    return pd.DataFrame(rows)


def combine_origins(frames, version):
    from scipy.stats import false_discovery_control
    frames = [f for f in frames if len(f)]
    if not frames:
        return pd.DataFrame(columns=["element", "gene_id", "control", "side", "dr_q", "ir_q", "at_q"])
    df = pd.concat(frames, ignore_index=True)
    for m in ["dr", "ir", "at"]:
        df[f"{m}_q"] = 1.0
        for ctl in (False, True):
            mask = df.control == ctl
            if mask.any():
                df.loc[mask, f"{m}_q"] = false_discovery_control(df.loc[mask, f"{m}_p"].values, method="bh")
    return df


def assign_evidence(catalogue, feats, origins_df, motifs, rules, th_name, version):
    ev, calls = evidence.assign_all(catalogue, feats, origins_df, rules["combination"],
                                    rules["thresholds"][th_name], th_name, motifs=motifs)
    return {"elements": ev, "initiators": calls}


def gene_evidence_table(gene_set, feats, general, rules, version):
    return annotate.gene_table(gene_set["genes"], feats, general, rules)


def scan_origin_motifs(fasta_path, checksum, rules, version):
    rows = []
    for name, _, seq in read_fasta(fasta_path):
        for mid, m in rules.get("origin_motifs", {}).items():
            for pos, strand, mm in orig.motif_hits(seq.upper(), m["sequence"], m["max_mismatches"]):
                rows.append({"element": name, "motif": mid, "position": pos, "strand": strand, "mismatches": mm,
                             "mechanisms": ";".join(m["mechanisms"]), "specificity": m["specificity"],
                             "family": m["family"]})
    return pd.DataFrame(rows, columns=["element", "motif", "position", "strand", "mismatches", "mechanisms",
                                       "specificity", "family"])


def enriched_cluster_origins(enriched, gene_set, fam, fasta_path, ocfg, q_max, version):
    """Origin statistics around members of unannotated clusters enriched in unresolved ECEs."""
    from scipy.stats import false_discovery_control
    top = enriched[enriched.q_bh <= q_max] if len(enriched) else enriched
    members = {}
    for r in top.itertuples():
        for gid, f in fam.items():
            if f == r.cluster:
                members.setdefault(gid.rsplit("_g", 1)[0], []).append(gid)
    seqs = {n: s.upper() for n, _, s in read_fasta(fasta_path)}
    frames = [element_origins(e, seqs[e], gene_set, {e: g}, ocfg, version) for e, g in sorted(members.items())]
    frames = [f for f in frames if len(f)]
    if not frames:
        return pd.DataFrame()
    df = pd.concat(frames, ignore_index=True)
    df["cluster"] = df.gene_id.map(fam)
    for m in ["dr", "ir", "at"]:
        df[f"{m}_q"] = 1.0
        for ctl in (False, True):
            mask = df.control == ctl
            if mask.any():
                df.loc[mask, f"{m}_q"] = false_discovery_control(df.loc[mask, f"{m}_p"].values, method="bh")
    return df
