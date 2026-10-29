"""Gene prediction on circular elements with explicit recovery of origin-spanning genes."""
import pyrodigal

MAX_OVERLAP = 60  # Prodigal's own maximum allowed overlap between genes
ROTATIONS = (0.5, 0.25, 0.75)


def _call(seq, min_gene):
    finder = pyrodigal.GeneFinder(meta=True, closed=False, min_gene=min_gene)
    return finder.find_genes(seq.encode())


def _gene_record(g, start, end, wraps, table, source):
    aa = g.translate(translation_table=table)
    return {
        "start": start,
        "end": end,
        "strand": "+" if g.strand == 1 else "-",
        "wraps_origin": wraps,
        "start_type": g.start_type,
        "rbs_motif": g.rbs_motif or "",
        "score": round(float(g.score), 2),
        "translation_table": table,
        "gc_cont": round(float(g.gc_cont), 4),
        "prediction_pass": source,
        "nt_seq": g.sequence(),
        "aa_seq": aa.rstrip("*"),
    }


def _span_nt(start, end, L):
    return end - start + 1 if end >= start else (L - start + 1) + end


def _positions(start, end, L):
    if end >= start:
        return set(range(start, end + 1))
    return set(range(start, L + 1)) | set(range(1, end + 1))


def _stop_key(r):
    return (r["strand"], r["end"] if r["strand"] == "+" else r["start"])


def predict_circular(element, seq, rotations=ROTATIONS, min_gene=90):
    """Predict genes on a circular sequence.

    Prodigal treats the FASTA ends as hard edges, so genes crossing the
    arbitrary start/end are reported as partial (or missed).  The sequence is
    therefore also predicted in rotated form; complete genes from rotated
    passes are added when they span the original boundary or replace a
    partial gene from the unrotated pass.  Duplicates are resolved by stop
    codon position and conflicting overlaps (> 60 nt) are rejected.

    Returns (genes, stats). Coordinates are 1-based inclusive on the original
    sequence; for genes with ``wraps_origin`` the end coordinate is smaller than
    the start coordinate.
    """
    L = len(seq)
    genes_a = _call(seq, min_gene)
    table = genes_a[0].translation_table if len(genes_a) else 11
    chosen, partial_regions = [], []
    for g in genes_a:
        if g.partial_begin or g.partial_end:
            partial_regions.append(("+" if g.strand == 1 else "-", _positions(g.begin, g.end, L)))
            continue
        chosen.append(_gene_record(g, g.begin, g.end, False, g.translation_table, "unrotated"))
    keys = {_stop_key(r) for r in chosen}
    chosen_pos = [_positions(r["start"], r["end"], L) for r in chosen]

    n_rescued = n_conflict = 0
    for frac in rotations:
        h = int(L * frac)
        rot = seq[h:] + seq[:h]
        boundary = L - h  # 0-based index in `rot` of original position 0
        for g in _call(rot, min_gene):
            if g.partial_begin or g.partial_end:
                continue
            b0, e0 = g.begin - 1, g.end - 1
            start, end = (b0 + h) % L + 1, (e0 + h) % L + 1
            spans = b0 < boundary <= e0
            rec = _gene_record(g, start, end, spans, g.translation_table, f"rotation_{frac}")
            if _stop_key(rec) in keys:
                continue
            pos = _positions(start, end, L)
            replaces_partial = any(s == rec["strand"] and len(pos & p) > MAX_OVERLAP for s, p in partial_regions)
            if not (spans or replaces_partial):
                continue
            clash = [i for i, cp in enumerate(chosen_pos) if len(pos & cp) > MAX_OVERLAP]
            if any(chosen[i]["prediction_pass"] != "unrotated" and chosen[i]["strand"] == rec["strand"]
                   for i in clash):
                continue  # already recovered from an earlier rotation
            if clash and spans:
                # A gene crossing the origin was called on unbroken sequence; it wins over
                # lower-scoring same-strand calls from the unrotated pass (typically alternative
                # starts/frames forced by the edge). Overlaps on the other strand or with
                # higher-scoring genes are kept: heavily overlapping genes are real in small phages.
                drop = [i for i in clash if chosen[i]["prediction_pass"] == "unrotated"
                        and chosen[i]["strand"] == rec["strand"] and chosen[i]["score"] < rec["score"]]
                for i in sorted(drop, reverse=True):
                    keys.discard(_stop_key(chosen[i]))
                    del chosen[i]
                    del chosen_pos[i]
                n_conflict += len(drop)
            elif clash:
                n_conflict += 1
                continue
            chosen.append(rec)
            chosen_pos.append(pos)
            keys.add(_stop_key(rec))
            n_rescued += 1

    genes = sorted(chosen, key=lambda r: (r["start"], r["end"]))
    for i, r in enumerate(genes, 1):
        r["element"] = element
        r["gene_id"] = f"{element}_g{i:04d}"
        r["gene_index"] = i
        r["length_nt"] = _span_nt(r["start"], r["end"], L)
        r["length_aa"] = len(r["aa_seq"])
    covered = set().union(*chosen_pos) if chosen_pos else set()
    stats = {
        "element": element,
        "n_genes": len(genes),
        "n_origin_spanning": sum(r["wraps_origin"] for r in genes),
        "n_rescued_from_rotation": n_rescued,
        "n_partial_in_unrotated": len(partial_regions),
        "n_conflicts_resolved": n_conflict,
        "coding_density": len(covered) / L if L else 0.0,
        "translation_table": table,
    }
    return genes, stats


def intergenic_regions(genes, L):
    """Return list of (start, end) gaps on a circular sequence (1-based; end < start wraps)."""
    if not genes:
        return [(1, L)]
    covered = bytearray(L + 1)
    for g in genes:
        for p in _positions(g["start"], g["end"], L):
            covered[p] = 1
    gaps = []
    p = 1
    while p <= L:
        if not covered[p]:
            q = p
            while q + 1 <= L and not covered[q + 1]:
                q += 1
            gaps.append([p, q])
            p = q + 1
        else:
            p += 1
    if len(gaps) >= 2 and gaps[0][0] == 1 and gaps[-1][1] == L:
        last = gaps.pop()
        gaps[0][0] = last[0]  # wraps
    return [tuple(g) for g in gaps]
