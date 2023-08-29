from strkit_rust_ext import get_snvs_dbsnp, get_snvs_meticulous, get_snvs_simple, get_read_snvs


REF_SEQ = "ACACACATGGCCATAC"
Q_SEQ = "ATCCCAAA"

#         A  A       T  T       C  G       C  C       C  C       A  A       A  A       A  C
PAIRS = [(0, 1000), (1, 1001), (2, 1003), (3, 1004), (4, 1005), (5, 1006), (6, 1008), (7, 1009)]
SNVS = ((1003, "C"), (1009, "A"))

SNV_CATALOG = [
    (1003, "snp1003", "G", ["C"]),
    (1009, "snp1009", "C", ["A"]),
]


def test_get_snvs_dbsnp():
    r_items = tuple(sorted(get_snvs_dbsnp(SNV_CATALOG, Q_SEQ, PAIRS, 994, 1000).items()))
    assert SNVS == r_items


def test_get_snvs_meticulous():
    r_items = tuple(sorted(get_snvs_meticulous(Q_SEQ, PAIRS, REF_SEQ, 994, 994, 1000, 0, 5, 10, 0.0).items()))
    assert SNVS == r_items


def test_get_snvs_simple():
    r_items = tuple(sorted(get_snvs_simple(Q_SEQ, PAIRS, REF_SEQ, 994, 994, 1000, 10, 0.0).items()))
    assert SNVS == r_items

def test_get_read_snvs():
    r_items = tuple(sorted(get_read_snvs(Q_SEQ, PAIRS, REF_SEQ, 994, 994, 1000, 0, 5, 20, 10, 0.0).items()))
    assert SNVS == r_items
