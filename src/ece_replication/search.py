"""Thin wrappers around HMMER, DIAMOND and MMseqs2 plus output parsers.

All tools are expected on PATH from the Pixi environment.
"""
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

from .util import project_path

# Fixed database size for E-value calculation, so values are identical
# regardless of how the query set is chunked (and comparable between the
# ECE and validation runs).
HMM_Z = 10000

DOMTBL_COLS = ["target", "target_acc", "tlen", "query", "query_acc", "qlen",
               "full_evalue", "full_score", "full_bias", "dom_n", "dom_of",
               "c_evalue", "i_evalue", "dom_score", "dom_bias",
               "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to", "acc"]


def run(cmd, **kw):
    res = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if res.returncode != 0:
        raise RuntimeError(f"command failed ({res.returncode}): {' '.join(map(str, cmd))}\n{res.stderr[-4000:]}")
    return res


def tool_version(cmd):
    try:
        res = subprocess.run(cmd, capture_output=True, text=True)
        return (res.stdout + res.stderr).strip().splitlines()[0]
    except Exception as e:  # pragma: no cover
        return f"unavailable: {e}"


def hmmsearch(hmm_path, faa_path, out_path, cpu=2, evalue=1e-3):
    out_path = project_path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = str(out_path) + ".part"
    run(["hmmsearch", "--noali", "-o", os.devnull, "--domtblout", tmp,
         "-E", str(evalue), "--domE", str(evalue), "-Z", str(HMM_Z), "--domZ", str(HMM_Z),
         "--cpu", str(cpu), str(project_path(hmm_path)), str(project_path(faa_path))])
    os.replace(tmp, out_path)
    return str(out_path)


def parse_domtbl(path):
    rows = []
    with open(project_path(path)) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split(None, 22)
            if len(parts) < 22:
                continue
            rows.append(parts[:22])
    df = pd.DataFrame(rows, columns=DOMTBL_COLS)
    for c in ["tlen", "qlen", "dom_n", "dom_of", "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to"]:
        df[c] = df[c].astype(int)
    for c in ["full_evalue", "full_score", "full_bias", "c_evalue", "i_evalue", "dom_score", "dom_bias", "acc"]:
        df[c] = df[c].astype(float)
    df = df.rename(columns={"target": "protein", "query": "model", "query_acc": "model_acc",
                            "qlen": "model_len", "tlen": "protein_len"})
    return df


def _union_len(intervals):
    total, cur_s, cur_e = 0, None, None
    for s, e in sorted(intervals):
        if cur_e is None or s > cur_e + 1:
            if cur_e is not None:
                total += cur_e - cur_s + 1
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    if cur_e is not None:
        total += cur_e - cur_s + 1
    return total


def summarise_domtbl(df):
    """Collapse domain rows to one row per (protein, model) with coverage summaries."""
    if df.empty:
        return pd.DataFrame(columns=["protein", "model", "model_acc", "model_len", "protein_len", "full_evalue",
                                     "full_score", "best_dom_score", "best_dom_ievalue", "hmm_coverage",
                                     "protein_coverage", "n_domains", "ali_from", "ali_to"])
    out = []
    for (prot, model), g in df.groupby(["protein", "model"], sort=False):
        sig = g[g.i_evalue <= 1e-2]
        use = sig if len(sig) else g
        best = g.loc[g.dom_score.idxmax()]
        out.append({
            "protein": prot, "model": model, "model_acc": best.model_acc,
            "model_len": int(best.model_len), "protein_len": int(best.protein_len),
            "full_evalue": float(best.full_evalue), "full_score": float(best.full_score),
            "best_dom_score": float(best.dom_score), "best_dom_ievalue": float(best.i_evalue),
            "hmm_coverage": _union_len(list(zip(use.hmm_from, use.hmm_to))) / best.model_len,
            "protein_coverage": _union_len(list(zip(use.env_from, use.env_to))) / best.protein_len,
            "n_domains": int(len(use)),
            "ali_from": int(use.env_from.min()), "ali_to": int(use.env_to.max()),
        })
    return pd.DataFrame(out)


DIAMOND_FIELDS = ["qseqid", "sseqid", "pident", "length", "qlen", "slen", "qstart", "qend",
                  "sstart", "send", "evalue", "bitscore"]


def diamond_blastp(ref_faa, query_faa, out_path, threads=4, mode="--ultra-sensitive", evalue=1e-3):
    out_path = project_path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir="/tmp") as tmp:
        db = Path(tmp) / "ref"
        run(["diamond", "makedb", "--in", str(project_path(ref_faa)), "-d", str(db), "--quiet"])
        cmd = ["diamond", "blastp", "-d", str(db), "-q", str(project_path(query_faa)), "-o", str(out_path) + ".part",
               "--evalue", str(evalue), "--max-target-seqs", "25", "--threads", str(threads), "--quiet",
               "--outfmt", "6", *DIAMOND_FIELDS]
        if mode:
            cmd.insert(2, mode)
        run(cmd)
    os.replace(str(out_path) + ".part", out_path)
    return str(out_path)


def parse_diamond(path):
    df = pd.read_csv(project_path(path), sep="\t", names=DIAMOND_FIELDS)
    df["qcov"] = (df.qend - df.qstart + 1) / df.qlen
    df["scov"] = (df.send - df.sstart + 1) / df.slen
    return df.rename(columns={"qseqid": "protein", "sseqid": "reference"})


def mmseqs_cluster(faa_path, out_prefix, min_seq_id=0.3, coverage=0.5, threads=4):
    out_prefix = project_path(out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir="/tmp") as tmp:
        run(["mmseqs", "easy-cluster", str(project_path(faa_path)), str(out_prefix), tmp,
             "--min-seq-id", str(min_seq_id), "-c", str(coverage), "--cov-mode", "0",
             "--cluster-mode", "0", "-s", "7.5", "--threads", str(threads), "-v", "1"])
    tsv = str(out_prefix) + "_cluster.tsv"
    df = pd.read_csv(tsv, sep="\t", names=["representative", "protein"])
    for suffix in ["_all_seqs.fasta", "_rep_seq.fasta"]:
        p = Path(str(out_prefix) + suffix)
        if p.exists():
            p.unlink()
    return df


def mmseqs_allvsall(faa_path, out_path, evalue=1e-3, threads=4):
    out_path = project_path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir="/tmp") as tmp:
        run(["mmseqs", "easy-search", str(project_path(faa_path)), str(project_path(faa_path)),
             str(out_path) + ".part", tmp, "-e", str(evalue), "-s", "7.5", "--threads", str(threads), "-v", "1",
             "--format-output", "query,target,pident,alnlen,qcov,tcov,evalue,bits"])
    os.replace(str(out_path) + ".part", out_path)
    return pd.read_csv(out_path, sep="\t", names=["query", "target", "pident", "alnlen", "qcov", "tcov",
                                                  "evalue", "bits"])


def mmseqs_nucleotide_search(query_fna, target_fna, out_path, evalue=1e-5, threads=4):
    out_path = project_path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir="/tmp") as tmp:
        run(["mmseqs", "easy-search", str(project_path(query_fna)), str(project_path(target_fna)),
             str(out_path) + ".part", tmp, "--search-type", "3", "-e", str(evalue), "--threads", str(threads),
             "-v", "1", "--format-output", "query,target,pident,alnlen,qstart,qend,tstart,tend,evalue,bits"])
    os.replace(str(out_path) + ".part", out_path)
    return pd.read_csv(out_path, sep="\t", names=["query", "target", "pident", "alnlen", "qstart", "qend",
                                                  "tstart", "tend", "evalue", "bits"])


def which_all(tools=("hmmsearch", "diamond", "mmseqs")):
    return {t: shutil.which(t) for t in tools}
