import numpy as np
from strkit_rust_ext import calculate_seq_with_wildcards, normalize_contig


def test_calculate_seq_with_wildcards():
    assert calculate_seq_with_wildcards(
        "ACCA",
        np.fromiter((5, 3, 3, 5), dtype=np.uint8),
        3,
    ) == "AXXA"


def test_normalize_contig():
    assert normalize_contig("chr5", True) == "chr5"
    assert normalize_contig("5", True) == "chr5"
    assert normalize_contig("X", True) == "chrX"
    assert normalize_contig("chr5", False) == "5"
    assert normalize_contig("chrX", False) == "X"
