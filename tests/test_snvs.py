from strkit_rust_ext import get_snvs_dbsnp, get_snvs_meticulous, get_snvs_simple, get_read_snvs
from .common import REF_SEQ, Q_SEQ, PAIRS, ALIGN_COORDS_Q, ALIGN_COORDS_R, SNVS, SNV_CATALOG


def test_get_snvs_dbsnp():
    r_items = tuple(sorted(get_snvs_dbsnp(SNV_CATALOG, Q_SEQ, PAIRS, 994, 1000).items()))
    assert SNVS == r_items


def test_get_snvs_meticulous():
    r_items = tuple(sorted(get_snvs_meticulous(Q_SEQ, REF_SEQ, ALIGN_COORDS_Q, ALIGN_COORDS_R, 994, 994, 1000, 0, 5, 10, 0.0).items()))
    assert SNVS == r_items


def test_get_snvs_simple():
    r_items = tuple(sorted(get_snvs_simple(Q_SEQ, REF_SEQ, ALIGN_COORDS_Q, ALIGN_COORDS_R, 994, 994, 1000, 20, 10, 0.0).items()))
    assert SNVS == r_items

def test_get_read_snvs():
    r_items = tuple(sorted(get_read_snvs(Q_SEQ, REF_SEQ, ALIGN_COORDS_Q, ALIGN_COORDS_R, 994, 994, 1000, 0, 5, 20, 10, 0.0).items()))
    assert SNVS == r_items
