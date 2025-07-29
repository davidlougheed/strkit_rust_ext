import logging
from datetime import datetime
from strkit_rust_ext import (
    shannon_entropy, get_aligned_pair_matches, consensus_seq, get_repeat_count
)
from .common import REF_SEQ, Q_QUALS, Q_SEQ, PAIRS, ALIGN_COORDS_Q, ALIGN_COORDS_R, CIGAR_OPS
from .cigar import get_aligned_pairs_from_cigar

logger = logging.getLogger(__name__)



strs = [
    b"ATGCGCGATAGAGCTAGTCGATGCC",
    b"ATGCTGATCGATCGGCGCGATATAC",
    b"AAAAAAAATTTTCCCTCTCTGGGAA",
    b"ATATTTTATATTTATTTATATATAT",
]
strs_st = tuple(map(lambda s: s.decode("ascii"), strs))

short_repeats = [
    "GTTTTTGTTTGTTTGTTTGTTTGTTTTTT",
    "GTTTTTTGTTTGTTTGTTTGTTTGTTTTTT",
    "GTTTTTTGTTTGTTTGTTTGTTTGTTTTTT",
    "GTTTTTTGTTTGTTTGTTTGTTTGTTTTTT",
    "GTTTTTTGTTTGTTTGTTTGTTTGTTTTTT",
    "GTTTTTTGTTTGTTTGTTTGTTTGTTTTTT",
    "GTTTTTTGTTTGTTTGTTTGTTTGTTTTTT",
    "GTTTTTTGTTTGTTTGTTTGTTTGTTTTTT",
    "GTTTTTTGTTTGTTTGTTTGTTTGTTTTTT",
    "GTTTTTTGTTTGTTTGTTTGTTTGTTTTTT",
    "GTTTTTTGTTTGTTTGTTTGTTTGTTTTTT",
    "GTTTTTTGTTTGTTTGTTTGTTTGTTTTTT",
    "GTTTTTTGTTTGTTTGTTTGTTTGTTTTTT",
    "GTTTTTTGTTTGTTTGTTTGTTTGTTTTTT",
    "GTTTTTTGTTTGTTTGTTTGTTTGTTTTTT",
    "GTTTTTTGTTTGTTTGTTTGTTTGTTTTTT",
    "GTTTTTTGTTTGTTTGTTTGTTTGTTTTTT",
    "GTTTTTTGTTTGTTTGTTTGTTTGTTTTTT",
    "GTTTTTTGTTTGTTTGTTTGTTTGTTTTTT",
    "GTTTTTTGTTTGTTTGTTTGTTTGTTTTTT",
    "GTTTTTTGTTTGTTTGTTTGTTTGTTTTTT",
    "GTTTTTTGTTTGTTTGTTTGTTTGTTTTTT",
    "GTTTTTTGTTTGTTTGTTTGTTTGTTTTTT",
    "GTTTTTTGTTTGTTTGTTTGTTTGTTTTTT",
    "GTTTTTTGTTTGTTTGTTTGTTTGTTTTTT",
    "GTTTTTTGTTTGTTTGTTTGTTTGTTTTTT",
    "GTTTTTTGTTTGTTTGTTTGTTTGTTTTTT",
    "GTTTTTTGTTTGTTTGTTTGTTTGTTTTTT",
    "GTTTTTTGTTTGTTTGTTTGTTTGTTTTTT",
    "GTTTTTTGTTTGTTTGTTTGTTTGTTTTTT",
    "GTTTTTTGTTTGTTTGTTTGTTTGTTTTTT",
    "GTTTTTTGTTTGTTTGTTTGTTTGTTTTTT",
    "GTTTTTTGTTTGTTTGTTTGTTTGTTTTTT",
    "GTTTTTTGTTTGTTTTGTTTGTTTGTTTTTT",
    "GTTTTTTTGTTTGTTTGTTTGTTTGTTTTTT",
    "GTTTTTTTGTTTGTTTGTTTGTTTGTTTTTT",
    "GTTTTTTGTTTGTTTGTTTGTTTGTTTTTTT",
    "GTTTTTTGTTTGTTTGTTTGTTTGTTTTTTT",
    "GTTTTTTGTTTTGTTTGTTTGTTTGTTTTTTT",
]

LONG_REPEATS = (("AAC" * 200, "AAC" * 220, "AAC" * 320, "AAC" * 335))

def main():
    dt = datetime.now()
    for _ in range(500000):
        for s in strs:
            shannon_entropy(s)
    print(f"shannon took {datetime.now() - dt}")

    # dt = datetime.now()
    # for _ in range(2000000):
    #     # query_sequence: str,
    #     # query_quals: list[int],
    #     # ref_seq: str,
    #     # query_coords: list[int],
    #     # ref_coords: list[int],
    #     # ref_coord_start: int,
    #     # tr_start_pos: int,
    #     # tr_end_pos: int,
    #     # contiguous_threshold: int,
    #     # max_snv_group_size: int,
    #     # too_many_snvs_threshold: int,
    #     # entropy_flank_size: int,
    #     # entropy_threshold: float,
    #     get_read_snvs(Q_SEQ, Q_QUALS, REF_SEQ, ALIGN_COORDS_Q, ALIGN_COORDS_R, 994, 994, 1000, 0, 5, 20, 10, 0.0)
    # print(f"get_read_snvs took {datetime.now() - dt}")

    dt = datetime.now()
    for _ in range(50000):
        list(get_aligned_pairs_from_cigar(CIGAR_OPS, 0, 100000, True))
    print(f"get_aligned_pairs_from_cigar (py) took {datetime.now() - dt}")

    dt = datetime.now()
    for _ in range(50000):
        get_aligned_pair_matches(CIGAR_OPS, 0, 100000)
    print(f"get_aligned_pair_matches (rs) took {datetime.now() - dt}")

    dt = datetime.now()
    print(consensus_seq(strs_st, logger, 100))
    for _ in range(10000):
        consensus_seq(strs_st, logger, 100)
    print(f"10000 iters of consensus_seq with short non-repeat sequences took {datetime.now() - dt}")

    dt = datetime.now()
    print(consensus_seq(short_repeats, logger, 100))
    for _ in range(10000):
        consensus_seq(short_repeats, logger, 100)
    print(f"10000 iters of consensus_seq with short repeats took {datetime.now() - dt}")

    n_iters = 100
    dt = datetime.now()
    for _ in range(n_iters):
        consensus_seq(LONG_REPEATS, logger, 100000)
    print(f"{n_iters} iters of consensus_seq with long repeats took {datetime.now() - dt}")

    dt = datetime.now()
    for _ in range(n_iters):
        consensus_seq(LONG_REPEATS, logger, 0)
    print(f"{n_iters} iters of consensus_seq (best_rep) with long repeats took {datetime.now() - dt}")

    dt = datetime.now()
    for _ in range(10000):
        rep = "AAC" * 56
        get_repeat_count(56, rep + "A", "GGC", "CGG", "AAC", 100, 3, 1)
        get_repeat_count(57, "AGC" + rep, "GGC", "CGG", "AAC", 100, 3, 1)
    print(f"10000 iters of get_repeat_count with slight imperfections took {datetime.now() - dt}")

    # dt = datetime.now()
    # for _ in range(1000):
    #     process_read_snvs_for_locus_and_calculate_useful_snvs(

    #     )


if __name__ == "__main__":
    main()
