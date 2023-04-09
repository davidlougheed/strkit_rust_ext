def find_pair_by_ref_pos(
    pairs: list[tuple[int, int]], 
    target: int,
) -> tuple[int, bool]: ...


def get_snvs_simple(
    query_sequence: str,
    pairs: list[tuple[int, int]],
    ref_seq: str,
    ref_coord_start: int,
    tr_start_pos: int,
    tr_end_pos: int,
) -> dict[int, str]: ...
