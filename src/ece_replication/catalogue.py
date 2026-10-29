"""Validate the ECE catalogue and compute per-element sequence statistics."""
import hashlib

import pandas as pd

from .util import open_maybe_gz, project_path, read_fasta

IUPAC_ACGT = set("ACGT")


def sequence_stats(seq):
    s = seq.upper()
    n = len(s)
    gc = s.count("G") + s.count("C")
    acgt = sum(s.count(b) for b in "ACGT")
    return {
        "length": n,
        "gc": gc / acgt if acgt else float("nan"),
        "n_ambiguous": n - acgt,
        "seq_md5": hashlib.md5(s.encode()).hexdigest(),
        "seq_sha256": hashlib.sha256(s.encode()).hexdigest(),
    }


def subgroup(category, note, dual_note):
    if category == "uncategorized" and isinstance(note, str) and note.strip() == dual_note:
        return "uncategorized_dual_marker"
    return category


def build_catalogue(fasta_path, props_path, expected, dual_note):
    """Return (catalogue DataFrame, {id: sequence}). Raises if identifiers/counts are inconsistent."""
    seqs, headers = {}, {}
    for name, desc, seq in read_fasta(fasta_path):
        if name in seqs:
            raise ValueError(f"duplicate FASTA identifier {name}")
        seqs[name] = seq.upper()
        headers[name] = desc
    with open_maybe_gz(project_path(props_path)) as f:
        props = pd.read_csv(f, sep="\t", dtype=str, keep_default_na=False)
    if props["Element"].duplicated().any():
        raise ValueError("duplicate identifiers in metadata")
    missing_seq = set(props["Element"]) - set(seqs)
    missing_meta = set(seqs) - set(props["Element"])
    if missing_seq or missing_meta:
        raise ValueError(f"identifier mismatch: no sequence for {sorted(missing_seq)}, "
                         f"no metadata for {sorted(missing_meta)}")
    rows = []
    for r in props.itertuples(index=False):
        st = sequence_stats(seqs[r.Element])
        # FASTA description repeats metadata: "<working_header> <sample> <contig> <category> [(<note>)]"
        expected_hdr = " ".join([r.Working_header, r.Sample, r.Contig, r.Category] +
                                ([f"({r.Note})"] if r.Note else []))
        header_consistent = headers[r.Element] == expected_hdr
        # The working header carries the size of the pre-dereplication cluster; it is recorded
        # (and compared) for transparency, not used as a validity criterion.
        size_tag = r.Working_header.rsplit("_Size_", 1)
        size_tag = int(size_tag[1]) if len(size_tag) == 2 else None
        rows.append({
            "element": r.Element,
            "working_header": r.Working_header,
            "sample": r.Sample,
            "contig": r.Contig,
            "category": r.Category,
            "note": r.Note,
            "subgroup": subgroup(r.Category, r.Note, dual_note),
            **st,
            "header_consistent": header_consistent,
            "working_header_size": size_tag,
            "length_minus_header_size": (st["length"] - size_tag) if size_tag is not None else None,
        })
    cat = pd.DataFrame(rows)
    counts = {
        "total": len(cat),
        "plasmid": int((cat.category == "plasmid").sum()),
        "virus": int((cat.category == "virus").sum()),
        "uncategorized": int((cat.category == "uncategorized").sum()),
        "dual_marker": int((cat.subgroup == "uncategorized_dual_marker").sum()),
    }
    for k, v in expected.items():
        if counts[k] != v:
            raise ValueError(f"expected {v} {k} elements, found {counts[k]}")
    if not cat.header_consistent.all():
        bad = cat[~cat.header_consistent].element.tolist()
        raise ValueError(f"FASTA headers disagree with metadata for {bad}")
    if cat.seq_sha256.duplicated().any():
        raise ValueError("identical sequences in a supposedly non-redundant catalogue")
    return cat, seqs
