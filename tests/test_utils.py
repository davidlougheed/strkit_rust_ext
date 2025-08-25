import numpy as np
from strkit_rust_ext import calculate_seq_with_wildcards


def test_calculate_seq_with_wildcards():
    assert calculate_seq_with_wildcards(
        "ACCA",
        np.fromiter((5, 3, 3, 5), dtype=np.uint8),
        3,
    ) == "AXXA"
