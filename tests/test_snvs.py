from strkit_rust_ext import find_pair_by_ref_pos, get_snvs_simple


REF_SEQ = "ACACACATGGCCATAC"
Q_SEQ = "ATCCCAAA"

#         A  A       T  T       C  G       C  C       C  C       A  A       A  A       A  C
PAIRS = [(0, 1000), (1, 1001), (2, 1003), (3, 1004), (4, 1005), (5, 1006), (6, 1008), (7, 1009)]
SNVS = ((1003, "C"), (1009, "A"))


def test_pair_bin_search():
    assert find_pair_by_ref_pos(PAIRS, 1004) == (3, True)
    assert find_pair_by_ref_pos(PAIRS, 1007) == (6, False)


def test_get_snvs_simple():
    r_items = tuple(sorted(get_snvs_simple(Q_SEQ, PAIRS, REF_SEQ, 994, 994, 1000).items()))
    assert SNVS == r_items
