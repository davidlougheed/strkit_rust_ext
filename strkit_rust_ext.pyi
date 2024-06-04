import numpy
from logging import Logger
from numpy.typing import NDArray
from typing import Literal, Optional, Sequence, Union

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
) -> tuple[Optional[tuple[NDArray[numpy.uint64], NDArray[numpy.uint64]]], int, int, int, int]: ...

def process_read_snvs_for_locus_and_calculate_useful_snvs(
    left_coord_adj: int,
    right_coord_adj: int,
    left_most_coord: int,
    ref_cache: str,
    read_dict_extra: dict[str, dict],
    read_q_coords: dict[str, NDArray[numpy.uint64]],
    read_r_coords: dict[str, NDArray[numpy.uint64]],
    candidate_snvs_dict: CandidateSNVs,
    min_allele_reads: int,
    significant_clip_snv_take_in: int,
    only_known_snvs: bool,
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


class STRkitBAMReader:
    references: list[str]

    def __init__(self, path: str, ref_path: str): ...

    def get_overlapping_segments_and_related_data(
        self,
        contig: str,
        left_coord: int,
        right_coord: int,
        max_reads: int,
        logger: Logger,
        locus_log_str: str,
    ) -> tuple[NDArray, int, NDArray[numpy.ulonglong], dict[str, int], int, int]: ...
