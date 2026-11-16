import cobra
import pytest

from build_draft_model import build_draft_model, remap_gpr, remap_gpr_detailed

MAP = {"1": "A", "2": "B", "3": "C", "4": "D"}


def _balanced(s: str) -> bool:
    depth = 0
    for ch in s:
        depth += ch == "("
        depth -= ch == ")"
        if depth < 0:
            return False
    return depth == 0


@pytest.mark.parametrize("rule, expected", [
    ("1", "A"),
    ("1 or 2", "A or B"),
    ("1 or X", "A"),
    ("1 and 2", "A and B"),
    ("(1 and 2) or 3", "(A and B) or C"),
])
def test_fully_or_partly_mapped(rule, expected):
    assert remap_gpr(rule, MAP) == expected


@pytest.mark.parametrize("rule", ["X", "X or Y", "1 and X"])
def test_strict_policy_rejects(rule):
    assert remap_gpr(rule, MAP, complex_policy="strict") is None


def test_mixed_rule_with_unmapped_complex_member_is_balanced():
    # DNApolym-like: an unmapped complex member next to unmapped isozymes
    rule = "(1 and 2 and X) or Y or 3 or 4"
    strict = remap_gpr(rule, MAP, complex_policy="strict")
    assert strict == "C or D"
    raven = remap_gpr(rule, MAP, complex_policy="raven")
    assert _balanced(raven)
    assert raven == "(A and B and OLD_X) or C or D"


def test_raven_policy_keeps_partial_complex_with_old_placeholder():
    new, info = remap_gpr_detailed("1 and 2 and X", MAP, complex_policy="raven")
    assert new == "A and B and OLD_X"
    assert info["complex_partial"] is True


def test_raven_policy_still_needs_one_mapped_gene():
    assert remap_gpr("X and Y", MAP, complex_policy="raven") is None
    assert remap_gpr("X", MAP, complex_policy="raven") is None


def test_orthogroup_fallback_used_only_when_no_one_to_one():
    fallback = {"X": ["F1", "F2"], "1": ["SHOULD_NOT_BE_USED"]}
    new, info = remap_gpr_detailed("1 and X", MAP, fallback=fallback)
    assert new == "A and (F1 or F2)"
    assert info["tiers"] == {"tier1", "tier2"}


def test_empty_rule():
    assert remap_gpr("", MAP) == ""


def _toy_template() -> cobra.Model:
    m = cobra.Model("toy")
    for rid, rule in [("R_or", "1 or X"), ("R_and", "1 and X"),
                      ("R_none", "X"), ("R_nogpr", "")]:
        r = cobra.Reaction(rid)
        m.add_reactions([r])
        r.gene_reaction_rule = rule
    return m


def test_build_draft_strict_raven_and_evidence():
    evidence = []
    draft = build_draft_model(_toy_template(), MAP, complex_policy="raven",
                              evidence_out=evidence)
    assert {r.id for r in draft.reactions} == {"R_or", "R_and", "R_nogpr"}
    status = {e["reaction_id"]: e["status"] for e in evidence}
    assert status == {"R_or": "tier1", "R_and": "complex_partial",
                      "R_none": "removed", "R_nogpr": "no_gpr_in_template"}


def test_build_draft_lenient_unchanged_network():
    draft = build_draft_model(_toy_template(), MAP, keep_unmapped=True)
    assert len(draft.reactions) == 4
    assert draft.reactions.R_none.gene_reaction_rule == ""
