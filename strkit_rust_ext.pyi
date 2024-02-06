from typing import Optional

def shannon_entropy(
    seq: bytes,
) -> float: ...

def get_snvs_dbsnp(
    candidate_snvs_dict_items_flat: list[tuple[int, str, str, list[str]]],
    query_sequence: str,
    pairs: list[tuple[int, int]],
    tr_start_pos: int,
    tr_end_pos: int,
) -> dict[int, str]: ...

def get_snvs_meticulous(
    query_sequence: str,
    ref_seq: str,
    query_coords: list[int],
    ref_coords: list[int],
    ref_coord_start: int,
    tr_start_pos: int,
    tr_end_pos: int,
    contiguous_threshold: int,
    max_snv_group_size: int,
    entropy_flank_size: int,
    entropy_threshold: float,
) -> dict[int, str]: ...

def get_snvs_simple(
    query_sequence: str,
    ref_seq: str,
    query_coords: list[int],
    ref_coords: list[int],
    ref_coord_start: int,
    tr_start_pos: int,
    tr_end_pos: int,
    too_many_snvs_threshold: int,
    entropy_flank_size: int,
    entropy_threshold: float,
) -> dict[int, str]: ...

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

def consensus_seq(
    seqs: tuple[str, ...],
) -> Optional[str]: ...
