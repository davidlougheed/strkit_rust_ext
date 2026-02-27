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

    def find_coord_idx_by_ref_pos(self, target: int, start_left: int) -> tuple[int, bool]: ...


# consensus

def consensus_seq(seqs: Sequence[str], logger: Logger, max_mdn_poa_length: int) -> Optional[tuple[str, Literal["single", "poa", "best_rep"]]]: ...

# exceptions

class LowMeanBaseQual(Exception):
    mean_base_qual: int

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

    flank_size: int

    annotations: list[str]

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
        annotations: list[str],
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

    @property
    def left_flank_coord(self) -> int: ...

    @property
    def right_flank_coord(self) -> int: ...

    def log_str(self) -> str: ...


class STRkitLocusBlock:
    contig: str

    def __init__(self, loci: list[STRkitLocus], left: int, right: int): ...
    def __len__(self) -> int: ...
    def __iter__(self): ...


def get_read_coords_from_matched_pairs(
    locus_with_ref_data: STRkitLocusWithRefData,
    segment: STRkitAlignedSegment,
    aligned_coords: STRkitAlignedCoords,
) -> tuple[int, int, int, int]: ...

def get_pairs_and_tr_read_coords(
    locus_with_ref_data: STRkitLocusWithRefData,
    segment: STRkitAlignedSegment,
) -> tuple[Optional[STRkitAlignedCoords], int, int, int, int]: ...

def process_read_snvs_for_locus_and_calculate_useful_snvs(
    block_segments: STRkitLocusBlockSegments,
    locus_with_ref_data: STRkitLocusWithRefData,
    # ---
    left_most_coord: int,
    ref_cache: str,
    # ---
    read_dict_extra: dict[str, dict],
    read_aligned_coords: dict[str, STRkitSegmentAlignmentDataForLocus],
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
    def __init__(self, path: str, is_sample_vcf: bool, sample_id: Optional[str]): ...

    def get_candidate_snvs(self, locus_block: STRkitLocusBlock) -> CandidateSNVs: ...


def shannon_entropy(
    seq: bytes,
) -> float: ...

def get_aligned_pair_matches(
    cigar: NDArray[numpy.uint32],
    query_start: int,
    ref_start: int,
) -> STRkitAlignedCoords: ...


# reads

class STRkitSegmentAlignmentDataForLocus:
    left_flank_end: int
    realigned: bool

    def __init__(
        self,
        aligned_coords: STRkitAlignedCoords,
        left_flank_start: int,
        left_flank_end: int,
        right_flank_start: int,
        right_flank_end: int,
        realigned: bool,
    ): ...


    def query_coord_at_idx(self, idx: int) -> int: ...
    def find_coord_idx_by_ref_pos(self, target: int, start_left: int) -> tuple[int, bool]: ...

class STRkitAlignedSegmentSequenceDataForLocus:
    flank_left_seq_wc: str
    flank_right_seq_wc: str
    tr_seq: str
    tr_seq_wc: str
    tr_len_with_flank: int

    def get_motif_size_kmers(self) -> list[str]: ...
    def get_est_copy_num(self) -> int: ...
    def calc_adj_score(self, read_cn_score: int) -> float: ...

class STRkitAlignedSegment:
    name: str
    length: int
    start: int
    end: int
    is_reverse: bool
    query_sequence: str
    query_qualities: NDArray[numpy.uint8]
    hp: Optional[int]
    ps: Optional[int]
    ps_remapped: Optional[int]

    def soft_clip_overlaps_locus(self, locus: STRkitLocus) -> bool: ...

    def get_sequence_data_for_locus(
        self,
        locus: STRkitLocus,
        segment_alignment_data_for_locus: STRkitSegmentAlignmentDataForLocus,
        min_avg_phred: float,
        base_wildcard_threshold: int,
    ) -> STRkitAlignedSegmentSequenceDataForLocus: ...

    def get_vcf_anchor_for_locus(
        self,
        locus_with_ref_data: STRkitLocusWithRefData,
        segment_alignment_data_for_locus: STRkitSegmentAlignmentDataForLocus,
        vcf_anchor_size: int,
    ) -> str: ...


class STRkitLocusSegmentsIter:
    def __iter__(self) -> "STRkitLocusSegmentsIter": ...
    def __next__(self) -> STRkitAlignedSegment: ...


class STRkitLocusSegments:
    n_segments: int
    sorted_read_lengths: NDArray[numpy.ulonglong]
    left_most_coord: int
    right_most_coord: int

    def __iter__(self) -> STRkitLocusSegmentsIter: ...

    def get_chimeric_read_status(self, rn: str) -> int: ...


class STRkitLocusBlockSegments:
    left_most_coord: int
    right_most_coord: int

    def get_segments_for_locus(
        self,
        locus: STRkitLocus,
    ) -> STRkitLocusSegments: ...

    def set_segment_ps_remapped(self, rn: str, ps_remapped: int): ...


class STRkitBAMReader:
    def __init__(
        self,
        path: str,
        ref_path: str,
        max_reads: int,
        skip_supp: bool,
        skip_sec: bool,
        use_hp: bool,
        significant_clip_threshold: int,
        logger: Logger,
        debug_logs: bool,
    ): ...

    def get_overlapping_segments_and_related_data_for_block(
        self,
        locus_block: STRkitLocusBlock,
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

def calculate_seq_with_wildcards(qs: str, quals: NDArray[numpy.uint8] | None, base_wildcard_threshold: int) -> str: ...

def normalize_contig(contig: str, has_chr: bool) -> str: ...
