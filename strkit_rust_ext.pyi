def get_snvs_simple(
    query_sequence: str,
    pairs: list[tuple[int, int]],
    ref_seq: str,
    ref_coord_start: int,
    tr_start_pos: int,
    tr_end_pos: int,
    entropy_flank_size: int,
    entropy_threshold: float,
) -> dict[int, str]: ...
