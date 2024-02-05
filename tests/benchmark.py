from datetime import datetime
from strkit_rust_ext import (
    shannon_entropy, get_snvs_simple, get_snvs_meticulous, get_read_snvs, get_aligned_pair_matches
)
from .common import REF_SEQ, Q_SEQ, PAIRS, ALIGN_COORDS_Q, ALIGN_COORDS_R, CIGAR_OPS
from .cigar import get_aligned_pairs_from_cigar



strs = [
    b"ATGCGCGATAGAGCTAGTCGATGCC",
    b"ATGCTGATCGATCGGCGCGATATAC",
    b"AAAAAAAATTTTCCCTCTCTGGGAA",
    b"ATATTTTATATTTATTTATATATAT",
]


def main():
    dt = datetime.now()
    for _ in range(100000):
        for s in strs:
            shannon_entropy(s)
    print(f"shannon took {datetime.now() - dt}")

    dt = datetime.now()
    for _ in range(5000000):
        get_snvs_simple(Q_SEQ, REF_SEQ, ALIGN_COORDS_Q, ALIGN_COORDS_R, 994, 994, 1000, 10, 0.0)
    print(f"get_snvs_simple took {datetime.now() - dt}")

    dt = datetime.now()
    for _ in range(5000000):
        get_snvs_meticulous(Q_SEQ, REF_SEQ, ALIGN_COORDS_Q, ALIGN_COORDS_R, 994, 994, 1000, 0, 5, 10, 0.0)
    print(f"get_snvs_meticulous took {datetime.now() - dt}")

    dt = datetime.now()
    for _ in range(5000000):
        get_read_snvs(Q_SEQ, REF_SEQ, ALIGN_COORDS_Q, ALIGN_COORDS_R, 994, 994, 1000, 0, 5, 20, 10, 0.0)
    print(f"get_read_snvs took {datetime.now() - dt}")

    dt = datetime.now()
    for _ in range(10000):
        list(get_aligned_pairs_from_cigar(CIGAR_OPS, 0, 100000, True))
    print(f"get_aligned_pairs_from_cigar (py) took {datetime.now() - dt}")

    dt = datetime.now()
    for _ in range(10000):
        get_aligned_pair_matches(CIGAR_OPS, 0, 100000)
    print(f"get_aligned_pair_matches (rs) took {datetime.now() - dt}")


if __name__ == "__main__":
    main()
