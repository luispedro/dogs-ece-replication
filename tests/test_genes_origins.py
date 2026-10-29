import numpy as np
import pytest

from ece_replication import genes, origins
from ece_replication.util import read_fasta, project_path

VAL = project_path("data/refs/validation_elements.fna")


def _val(name):
    if not VAL.exists():
        pytest.skip("validation genomes not fetched")
    for n, _, s in read_fasta(VAL):
        if n == name:
            return s.upper()
    raise KeyError(name)


def test_origin_spanning_gene_recovered():
    seq = _val("pT181")
    ref, _ = genes.predict_circular("x", seq)
    longest = max(ref, key=lambda g: g["length_aa"])
    assert not longest["wraps_origin"]
    # rotate so the cut falls in the middle of the longest gene
    mid = (longest["start"] + longest["end"]) // 2
    rot = seq[mid:] + seq[:mid]
    pred, stats = genes.predict_circular("x", rot)
    assert stats["n_origin_spanning"] >= 1
    proteins = {g["aa_seq"] for g in pred}
    assert longest["aa_seq"] in proteins
    wrapped = [g for g in pred if g["aa_seq"] == longest["aa_seq"]][0]
    assert wrapped["wraps_origin"] and wrapped["end"] < wrapped["start"]
    assert wrapped["length_nt"] == longest["length_nt"]


def test_rotation_invariance_of_gene_set():
    seq = _val("pC194")
    a, _ = genes.predict_circular("x", seq)
    b, _ = genes.predict_circular("x", seq[1000:] + seq[:1000])
    pa = {g["aa_seq"] for g in a}
    pb = {g["aa_seq"] for g in b}
    assert len(pa & pb) >= 0.8 * max(len(pa), len(pb))


def test_intergenic_wraps():
    g = [{"start": 10, "end": 50, "strand": "+"}, {"start": 60, "end": 90, "strand": "-"}]
    gaps = genes.intergenic_regions(g, 100)
    assert (91, 9) in gaps and (51, 59) in gaps


def test_flanking_regions_strand():
    gaps = [(91, 9), (51, 59)]
    plus = {"start": 10, "end": 50, "strand": "+"}
    sides = {s: (a, b) for s, a, b in origins.flanking_regions(plus, gaps, 100, 1000)}
    assert sides["upstream"] == (91, 9) and sides["downstream"] == (51, 59)
    minus = {"start": 60, "end": 90, "strand": "-"}
    sides = {s: (a, b) for s, a, b in origins.flanking_regions(minus, gaps, 100, 1000)}
    assert sides["upstream"] == (91, 9) and sides["downstream"] == (51, 59)


def test_repeat_statistics_detect_planted_signals():
    rng = np.random.default_rng(0)
    bg = "".join(rng.choice(list("ACGT"), 400))
    iteron = "TGTGACAAATTGCCCTCAG"
    seq = bg[:100] + iteron * 4 + bg[100:]
    assert origins.dr_copies(origins.encode(seq), 10) >= 4
    stem = "GCGCATTAGGCC"
    rc = stem[::-1].translate(str.maketrans("ACGT", "TGCA"))
    hp = bg[:200] + stem + "TTTT" + rc + bg[200:]
    assert origins.ir_stem(origins.encode(hp), 3, 20, 30) >= len(stem)
    cfg = dict(dr_kmer=10, ir_loop_min=3, ir_loop_max=20, ir_max_stem=30, at_window=60, n_background=50)
    res = origins.test_region(seq, bg * 5, cfg, seed=1)
    assert res["dr_p"] < 0.05


def test_circular_slice():
    assert origins.circular_slice("ABCDEFGHIJ", 9, 2) == "IJAB"
    assert origins.circular_slice("ABCDEFGHIJ", 2, 4) == "BCD"


def test_phix174_gene_a_spans_origin():
    # phiX174 gene A (513 aa) crosses the GenBank origin and overlaps other genes
    g, _ = genes.predict_circular("phiX174", _val("phiX174"))
    a = [x for x in g if x["length_aa"] == 513]
    assert len(a) == 1 and a[0]["wraps_origin"]


def test_cole1_origin_motif():
    from ece_replication.util import load_yaml
    m = load_yaml("config/replication_rules.yaml")["origin_motifs"]["ColE1_RNAI_RNAII"]
    hits = origins.motif_hits(_val("ColE1"), m["sequence"], m["max_mismatches"])
    assert any(p == 753 and s == "-" for p, s, _ in hits)
    assert not origins.motif_hits(_val("pT181"), m["sequence"], m["max_mismatches"])
