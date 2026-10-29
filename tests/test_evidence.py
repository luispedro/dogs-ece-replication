import pandas as pd

from ece_replication import evidence
from ece_replication.util import load_yaml

RULES = load_yaml("config/replication_rules.yaml")
TH = RULES["thresholds"]["default"]
COMB = RULES["combination"]


def feat(gene, element, role, mechs, spec, source="pfam", cov=0.9, lookalike=False, name="X", module=None):
    mod = module or next((m for m, roles in RULES["modules"].items() if role in roles), "other")
    return {"gene_id": gene, "element": element, "source": source, "family_id": name, "family_name": name,
            "role": role, "mechanisms": mechs, "specificity": spec, "score": 100, "evalue": 1e-30,
            "coverage": cov, "identity": 40.0 if source == "diamond" else float("nan"),
            "transfer_lookalike": lookalike, "detail": name, "module": mod}


def run(rows, length=5000):
    f = pd.DataFrame(rows)
    cat = pd.DataFrame({"element": ["E"], "length": [length], "category": ["plasmid"], "subgroup": ["plasmid"]})
    ev, _ = evidence.assign_all(cat, f, None, COMB, TH, "default")
    return ev.iloc[0]


def test_strong_rcr_with_reference():
    r = run([feat("g1", "E", "initiator", "rolling_circle", "high", name="Rep_1"),
             feat("g1", "E", "initiator", "rolling_circle", "medium", source="diamond", name="pT181 RepC")])
    assert r.mechanism == "rolling_circle" and r.tier == "strong"


def test_uncorroborated_is_tentative():
    r = run([feat("g1", "E", "initiator", "theta", "high", name="Rep3_N")])
    assert r.mechanism == "theta" and r.tier == "tentative"


def test_incompatible_initiators_conflict():
    r = run([feat("g1", "E", "initiator", "rolling_circle", "high", name="Rep_1"),
             feat("g2", "E", "initiator", "theta", "high", name="Rep3_N")])
    assert r.tier == "conflicting"


def test_relaxase_alone_is_not_replication():
    r = run([feat("g1", "E", "relaxase", "-", "high", name="Mob_Pre")])
    assert r.tier == "unresolved" and r.mechanism == "none_detected"
    assert "transfer_genes_without_replication_genes" in r["flags"]


def test_rep_trans_ignored_in_conjugative_context():
    r = run([feat("g1", "E", "initiator", "rolling_circle", "medium", name="Rep_trans", lookalike=True),
             feat("g2", "E", "coupling_protein", "-", "high", name="T4CP"),
             feat("g3", "E", "t4ss", "-", "high", name="VirB4")])
    assert r.mechanism != "rolling_circle"
    assert "rep_trans_in_transfer_context" in r["flags"]


def test_multimode_viral_intersection():
    r = run([feat("g1", "E", "initiator", "theta;rolling_circle", "high", name="Phage_rep_O"),
             feat("g2", "E", "initiator", "rolling_circle", "high", name="Rep_1")])
    assert r.mechanism == "rolling_circle"


def test_large_rcr_flagged_not_strong():
    r = run([feat("g1", "E", "initiator", "rolling_circle", "high", name="Rep_1"),
             feat("g1", "E", "initiator", "rolling_circle", "medium", source="diamond", name="ref")], length=80000)
    assert r.tier == "tentative" and "rcr_initiator_on_large_element" in r["flags"]


def test_replisome_only_unassigned():
    r = run([feat("g1", "E", "primase", "-", "low", name="PriCT_1")])
    assert r.mechanism == "unassigned" and r.tier == "tentative"


def test_large_rcr_phage_allowed():
    r = run([feat("g1", "E", "initiator", "rolling_circle", "high", name="Phage_GPA"),
             feat("g1", "E", "initiator", "rolling_circle", "medium", source="diamond", name="P2 A"),
             feat("g2", "E", "portal", "-", "high", name="Phage_portal"),
             feat("g3", "E", "terminase", "-", "high", name="TerL")], length=33000)
    assert r.tier == "strong" and "rcr_initiator_on_large_element" not in r["flags"]
