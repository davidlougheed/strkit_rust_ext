import numpy
from logging import Logger
from numpy.typing import NDArray
from typing import Literal, Optional, Sequence, Union

# aligned_coords

class STRkitAlignedCoords:
    query_coords: list[int]
    ref_coords: list[int]

    def __init__(self, query_coords: NDArray[numpy.uint64], ref_coords: NDArray[numpy.uint64]): ...

    def query_coord_at_idx(self, idx: int) -> int: ...


# consensus

def consensus_seq(seqs: Sequence[str], logger: Logger, max_mdn_poa_length: int) -> Optional[tuple[str, Literal["single", "poa", "best_rep"]]]: ...

# locus

class STRkitLocus:
    t_idx: int
    locus_id: str

    contig: str

    left_coord: int
    left_flank_coord: int  # left_coord - flank_size
    right_coord: int
    right_flank_coord: int  # right_coord + flank_size
    ref_size: int

    motif: str
    motif_size: int

    n_alleles: int

    def __init__(
        self,
        t_idx: int,
        locus_id: str,
        contig: str,
        left_coord: int,
        right_coord: int,
        motif: str,
        n_alleles: int,
        flank_size: int,
    ): ...

    def log_str(self) -> str: ...

    def to_dict(self) -> dict: ...

    def with_ref_data(
        self,
        left_coord_adj: int,
        right_coord_adj: int,
        ref_contig: str,
        ref_cn: int,
        ref_seq: str,
        ref_left_flank_seq: str,
        ref_right_flank_seq: str,
        ref_total_seq: str,
        ref_time: float,
    ) -> "STRkitLocusWithRefData":
        ...


class STRkitLocusWithRefData:
    locus_def: STRkitLocus

    left_coord_adj: int
    right_coord_adj: int

    ref_contig: str
    ref_cn: int
    ref_seq: str
    ref_left_flank_seq: str
    ref_right_flank_seq: str
    ref_total_seq: str
    ref_time: float


def get_read_coords_from_matched_pairs(
    locus_with_ref_data: STRkitLocusWithRefData,
    query_seq: str,
    aligned_coords: STRkitAlignedCoords,
) -> tuple[int, int, int, int]: ...

def get_pairs_and_tr_read_coords(
    cigar: NDArray[numpy.uint32],
    segment_start: int,
    locus_with_ref_data: STRkitLocusWithRefData,
    query_seq: str,
) -> tuple[Optional[STRkitAlignedCoords], int, int, int, int]: ...

def process_read_snvs_for_locus_and_calculate_useful_snvs(
    left_coord_adj: int,
    right_coord_adj: int,
    # ---
    left_most_coord: int,
    ref_cache: str,
    # ---
    read_dict_extra: dict[str, dict],
    read_aligned_coords: dict[str, STRkitAlignedCoords],
    candidate_snvs: CandidateSNVs,
    min_allele_reads: int,
    significant_clip_snv_take_in: int,
    only_known_snvs: bool,
    # ---
    logger: Logger,
    locus_log_str: str,
) -> list[tuple[int, int]]: ...

# snvs

class CandidateSNVs:
    def get(self, pos: int) -> Optional[dict]: ...


class STRkitVCFReader:
    def __init__(self, path: str): ...

    def get_candidate_snvs(
        self,
        snv_vcf_contigs: tuple[str, ...],
        snv_vcf_file_format: Literal["chr", "num", "acc", ""],
        contig: str,
        left_most_coord: int,
        right_most_coord: int,
    ): ...


def shannon_entropy(
    seq: bytes,
) -> float: ...

def get_aligned_pair_matches(
    cigar: NDArray[numpy.uint32],
    query_start: int,
    ref_start: int,
) -> STRkitAlignedCoords: ...


# reads

class STRkitAlignedSegment:
    name: str
    length: int
    start: int
    end: int
    is_reverse: bool
    query_sequence: str
    query_qualities: NDArray[numpy.uint8]
    raw_cigar: NDArray[numpy.uint32]
    hp: Optional[int]
    ps: Optional[int]


class STRkitLocusBlockSegments:
    left_most_coord: int
    right_most_coord: int

    def get_segments_for_locus(
        self,
        locus: STRkitLocus,
    ) -> tuple[NDArray, int, NDArray[numpy.ulonglong], dict[str, int], int, int]: ...


class STRkitBAMReader:
    references: list[str]

    def __init__(
        self,
        path: str,
        ref_path: str,
        max_reads: int,
        skip_supp: bool,
        skip_sec: bool,
        use_hp: bool,
        logger: Logger,
        debug_logs: bool,
    ): ...

    def get_overlapping_segments_and_related_data_for_block(
        self,
        contig: str,
        left_coord: int,
        right_coord: int,
        log_str: str,
    ) -> STRkitLocusBlockSegments: ...


# repeats

def get_repeat_count(
    start_count: int,
    tr_seq: str,
    flank_left_seq: str,
    flank_right_seq: str,
    motif: str,
    max_iters: int,
    local_search_range: int,
    step_size: int,
    use_shortcuts: bool,
) -> tuple[tuple[int, int], int, int]: ...

# utils

def find_coord_idx_by_ref_pos(aligned_coords: STRkitAlignedCoords, target: int, start_left: int) -> tuple[int, bool]: ...

def calculate_seq_with_wildcards(qs: str, quals: NDArray[numpy.uint8] | None, base_wildcard_threshold: int) -> str: ...
