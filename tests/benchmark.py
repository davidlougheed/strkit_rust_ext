from datetime import datetime
from strkit_rust_ext import (
    shannon_entropy, get_read_snvs, get_aligned_pair_matches, consensus_seq, 
    process_read_snvs_for_locus_and_calculate_useful_snvs
)
from .common import REF_SEQ, Q_SEQ, PAIRS, ALIGN_COORDS_Q, ALIGN_COORDS_R, CIGAR_OPS
from .cigar import get_aligned_pairs_from_cigar



strs = [
    b"ATGCGCGATAGAGCTAGTCGATGCC",
    b"ATGCTGATCGATCGGCGCGATATAC",
    b"AAAAAAAATTTTCCCTCTCTGGGAA",
    b"ATATTTTATATTTATTTATATATAT",
]
strs_st = tuple(map(lambda s: s.decode("ascii"), strs))

def main():
    dt = datetime.now()
    for _ in range(500000):
        for s in strs:
            shannon_entropy(s)
    print(f"shannon took {datetime.now() - dt}")

    dt = datetime.now()
    for _ in range(2000000):
        get_read_snvs(Q_SEQ, REF_SEQ, ALIGN_COORDS_Q, ALIGN_COORDS_R, 994, 994, 1000, 0, 5, 20, 10, 0.0)
    print(f"get_read_snvs took {datetime.now() - dt}")

    dt = datetime.now()
    for _ in range(50000):
        list(get_aligned_pairs_from_cigar(CIGAR_OPS, 0, 100000, True))
    print(f"get_aligned_pairs_from_cigar (py) took {datetime.now() - dt}")

    dt = datetime.now()
    for _ in range(50000):
        get_aligned_pair_matches(CIGAR_OPS, 0, 100000)
    print(f"get_aligned_pair_matches (rs) took {datetime.now() - dt}")

    dt = datetime.now()
    print(consensus_seq(strs_st))
    for _ in range(10000):
        consensus_seq(strs_st)
    print(f"consensus_seq took {datetime.now() - dt}")

    # dt = datetime.now()
    # for _ in range(1000):
    #     process_read_snvs_for_locus_and_calculate_useful_snvs(

    #     )


if __name__ == "__main__":
    main()
