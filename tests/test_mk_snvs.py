from strkit_rust_ext import mk_snvs_dict, mk_snvs_hash


REF_SEQ = "ACACACATGGCCATAC"
Q_SEQ = "ATCCCAAA"

#         A  A       T  T       C  G       C  C       C  C       A  A       A  A       A  C
PAIRS = [(0, 1000), (1, 1001), (2, 1003), (3, 1004), (4, 1005), (5, 1006), (6, 1008), (7, 1009)]
SNVS = ((1003, "C"), (1009, "A"))


def test_mk_snvs_dict():
    r_items = tuple(sorted(mk_snvs_dict(Q_SEQ, PAIRS, REF_SEQ, 994, 994, 1000).items()))
    assert SNVS == r_items


def test_mk_snvs_hash():
    r_items = tuple(sorted(mk_snvs_hash(Q_SEQ, PAIRS, REF_SEQ, 994, 994, 1000).items()))
    assert SNVS == r_items
