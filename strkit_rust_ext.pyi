import numpy.typing
from logging import Logger
from typing import Optional, Sequence, Union

# consensus

def best_representatives(seqs: Sequence[str]) -> set[str]: ...
def consensus_seq(seqs: Sequence[str]) -> Optional[str]: ...

# locus

def get_pairs_and_tr_read_coords(
    cigar: list[tuple[int, int]],
    segment_start: int,
    left_flank_coord: int,
    left_coord: int,
    right_coord: int,
    right_flank_coord: int,
    motif: str,
    motif_size: int,
    query_seq: str,
) -> tuple[Optional[tuple[numpy.typing.NDArray, numpy.typing.NDArray]], int, int, int, int]: ...

def process_read_snvs_for_locus_and_calculate_useful_snvs(
    left_coord_adj: int,
    right_coord_adj: int,
    left_most_coord: int,
    ref_cache: str,
    read_dict_extra: dict[str, dict],
    read_q_coords: dict[str, numpy.typing.NDArray],
    read_r_coords: dict[str, numpy.typing.NDArray],
    candidate_snvs_dict: dict[int, dict[str, Union[str, tuple[str, ...]]]],
    min_allele_reads: int,
    significant_clip_snv_take_in: int,
    only_known_snvs: bool,
    logger: Logger,
    locus_log_str: str,
) -> list[tuple[int, int]]: ...

# snvs

def shannon_entropy(
    seq: bytes,
) -> float: ...

def get_read_snvs(
    query_sequence: str,
    ref_seq: str,
    query_coords: list[int],
    ref_coords: list[int],
    ref_coord_start: int,
    tr_start_pos: int,
    tr_end_pos: int,
    contiguous_threshold: int,
    max_snv_group_size: int,
    too_many_snvs_threshold: int,
    entropy_flank_size: int,
    entropy_threshold: float,
) -> dict[int, str]: ...

def get_aligned_pair_matches(
    cigar: list[tuple[int, int]],
    query_start: int,
    ref_start: int,
) -> tuple[list[int], list[int]]: ...
