from strkit_rust_ext import get_aligned_pair_matches
from .common import CIGAR_OPS
from .cigar import get_aligned_pairs_from_cigar


def test_cigar():
    pairs1 = list(get_aligned_pairs_from_cigar(CIGAR_OPS, 0, 100000, True))
    pairs2 = list(zip(*get_aligned_pair_matches(CIGAR_OPS, 0, 100000)))
    assert pairs1 == pairs2
