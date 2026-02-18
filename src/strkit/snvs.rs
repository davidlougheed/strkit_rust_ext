use bytecount;
use pyo3::exceptions::PyException;
use pyo3::exceptions::PyValueError;
use pyo3::intern;
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyDict, PyString};
use regex::Regex;
use rust_htslib::bcf;
use rust_htslib::bcf::Read;
use std::cmp;
use std::collections::{HashMap, HashSet};
use std::sync::Mutex;

use crate::aligned_coords::STRkitAlignedCoords;
use crate::locus::STRkitLocusBlock;
use crate::reads::{STRkitAlignedSegment, STRkitLocusBlockSegments, STRkitSegmentAlignmentDataForLocus};

static SNV_OUT_OF_RANGE_CHAR: char = '-';
static SNV_GAP_CHAR: char = '_';


pub struct UsefulSNVsParams {
    pub contiguous_threshold: usize,
    pub max_snv_group_size: usize,
    pub too_many_snvs_threshold: usize,
    pub entropy_flank_size: usize,
    pub entropy_threshold: f32,
}


#[derive(Clone)] // TODO: rust: zero-copy when possible
#[pyclass(frozen)]
pub struct CandidateSNV {
    id: String,
    ref_base: char,
    alts: Vec<char>,
    // This struct must be built with STRkitVCFReader.get_candidate_snvs(...) via CandidateSNVs
}

#[pyclass(frozen)]
pub struct CandidateSNVs {
    pub snvs: HashMap<usize, CandidateSNV>,
}

#[pymethods]
impl CandidateSNVs {
    fn get<'py>(&self, py: Python<'py>, pos: usize) -> Option<Bound<'py, PyDict>> {
        self.snvs.get(&pos).map(move |c| {
            let res = PyDict::new(py);
            res.set_item("id", &c.id).unwrap();
            res.set_item("ref_base", c.ref_base).unwrap();
            res.set_item("alts", &c.alts).unwrap();
            res
        })
    }
}


fn _human_chrom_to_refseq_accession<'x>(contig: &str, snv_vcf_contigs: &'x Vec<String>) -> Option<&'x str> {
    let mut c = contig;
    c = c.strip_prefix("chr").unwrap_or(c);
    match c {
        "X" => { c = "23"; }
        "Y" => { c = "24"; }
        "M" => { c = "12920"; }
        _ => {}
    };

    let nc_fmt: String = format!("NC_{:06}", c);
    c = nc_fmt.as_str();

    let mut ret: Option<&str> = None;

    snv_vcf_contigs.iter().for_each(|vcf_contig| {
        if vcf_contig.starts_with(c) {
            ret = Some(vcf_contig);  // ret = vcf_contig;  // .to_string();
        }
    });

    ret
}

enum VCFContigFormat {
    Chr, // UCSC style (chr1, chr2, ...)
    Num, // Numeral contig pattern (1, 2, ..., X, Y, M)
    Acc, // Accession pattern (NC_)
    Unk, // Unknown
}

const NUMERAL_CONTIG_PATTERN: &str = r"^(\d{1,2}|X|Y)$";
const ACCESSION_PATTERN: &str = r"^NC_\d+";

fn _get_vcf_contig_format(snv_vcf_contigs: &Vec<String>) -> VCFContigFormat {
    if snv_vcf_contigs.is_empty() || snv_vcf_contigs[0].starts_with("chr") {
        VCFContigFormat::Chr
    } else if Regex::new(NUMERAL_CONTIG_PATTERN).unwrap().is_match(&snv_vcf_contigs[0]) {
        VCFContigFormat::Num
    } else if Regex::new(ACCESSION_PATTERN).unwrap().is_match(&snv_vcf_contigs[0]) {
        VCFContigFormat::Acc
    } else {
        VCFContigFormat::Unk
    }
}


#[pyclass]
pub struct STRkitVCFReader {
    reader: Mutex<bcf::IndexedReader>,
    contig_names: Vec<String>,
    contig_format: VCFContigFormat,
    // If sample_idx is Some(...), then this is a sample VCF with sample at index <...>.
    // Otherwise, this is not a sample VCF but merely a catalog.
    sample_idx: Option<usize>,
}

#[pymethods]
impl STRkitVCFReader {
    #[new]
    fn py_new(path: &str, is_sample_vcf: bool, sample_id: Option<String>) -> PyResult<Self> {
        let r = bcf::IndexedReader::from_path(path);

        if let Ok(rdr) = r {
            // fetch SNV VCF contigs
            let header = rdr.header();
            let contig_names: Vec<String> = header.header_records().into_iter().filter_map(|h: bcf::HeaderRecord| {
                match h {
                    bcf::HeaderRecord::Contig { key: _, values } => {
                        if let Some(id) = values.get("ID") {
                            Some(id.clone())
                        } else {
                            None
                        }
                    },
                    _ => None
                }
            }).collect();

            let mut sample_idx: Option<usize> = None;

            if is_sample_vcf {
                sample_idx = if let Some(sample_id_str) = sample_id {
                    header.sample_id(sample_id_str.as_bytes())
                } else {
                    Some(0)
                };

                if sample_idx.is_none() {
                    // Still None --> couldn't find sample in header
                    return Err(PyException::new_err("Could not find sample in sample VCF"));
                }
            }

            let contig_format = _get_vcf_contig_format(&contig_names);

            let reader = Mutex::new(rdr);
            Ok(STRkitVCFReader { reader, contig_names, contig_format, sample_idx })
        } else {
            Err(PyErr::new::<PyValueError, _>(format!("Could not load VCF from path: {}", path)))
        }
    }

    fn get_candidate_snvs<'py>(
        &mut self,
        py: Python<'py>,
        locus_block: &STRkitLocusBlock,
    ) -> PyResult<Bound<'py, CandidateSNVs>> {
        let mut reader = self.reader.lock().unwrap();

        let header = reader.header();

        let mut candidate_snvs = HashMap::<usize, CandidateSNV>::new();

        let lb_contig: &str = &locus_block.loci[0].contig;
        let snv_contig: &str = if header.name2rid(lb_contig.as_bytes()).is_err() {
            match self.contig_format {
                VCFContigFormat::Num => lb_contig.strip_prefix("chr").unwrap_or(lb_contig),
                VCFContigFormat::Acc => _human_chrom_to_refseq_accession(lb_contig, &self.contig_names).unwrap(),
                _ => lb_contig, // Otherwise, leave as-is (will cause an error below with name2rid)
            }
        } else {
            lb_contig // Exists, leave as-is
        };

        match header.name2rid(snv_contig.as_bytes()) {
            Ok(contig_rid) => {
                reader.fetch(contig_rid, locus_block.left, Some(locus_block.right + 1)).unwrap();

                let mut record = reader.empty_record();

                while let Some(r) = reader.read(&mut record) {
                    match r {
                        Ok(_) => {
                            let alleles = record.alleles();

                            let snv_ref = alleles[0];
                            if snv_ref.len() == 1 {
                                let snv_alts = alleles[1..]
                                    .iter()
                                    .filter(|a| a.len() == 1)
                                    .map(|aa| aa[0] as char)
                                    .collect::<Vec<char>>();

                                if !snv_alts.is_empty() {
                                    candidate_snvs.insert(record.pos() as usize, CandidateSNV {
                                        id: String::from_utf8(record.id()).unwrap(),
                                        ref_base: snv_ref[0] as char,
                                        alts: snv_alts,
                                    });
                                }
                            }
                        }
                        Err(_) => {
                            return Err(PyErr::new::<PyException, _>("Reading VCF record failed"));
                        }
                    }
                }

                Bound::new(py, CandidateSNVs { snvs: candidate_snvs })
            },
            Err(_) => {
                Err(PyErr::new::<PyException, _>(format!("Could not find contig in VCF: {}", lb_contig)))
            }
        }
    }
}


// We check entropy against a threshold in order to make sure the SNVs we find are somewhat
// useful, i.e., surrounded by a nice mixture of different bases rather than some possibly
// mis-mappable base inside a homopolymer (I used to get results like AAAAAAAAAACAAAAAAAA
// where one of those As could be a gap instead... really useless stuff that made it
// through the filters.)

fn _byte_entropy_f32(data: &[u8], data_len: f32, byte: u8) -> f32 {
    let mut p = bytecount::count(data, byte) as f32;
    if p > 0.0f32 {
        p /= data_len;
        p *= p.log2();
    }
    p
}

/// Calculates the Shannon entropy of a slice of bytes representing DNA bases. This slice should only consist of
/// ATGCatgc in ASCII form and no other bytes.
fn _shannon_entropy_dna(data: &[u8]) -> f32 {
    let data_len = data.len() as f32;

    // Calculate entropy based on [atgcATGC] content:
    let entropy: f32 = 0.0
        - _byte_entropy_f32(data, data_len, 65u8) // A
        - _byte_entropy_f32(data, data_len, 67u8) // C
        - _byte_entropy_f32(data, data_len, 71u8) // G
        - _byte_entropy_f32(data, data_len, 84u8) // T
        - _byte_entropy_f32(data, data_len, 97u8) // a
        - _byte_entropy_f32(data, data_len, 99u8) // c
        - _byte_entropy_f32(data, data_len, 103u8) // g
        - _byte_entropy_f32(data, data_len, 116u8); // t

    entropy
}

#[pyfunction]
pub fn shannon_entropy(data: &Bound<'_, PyBytes>) -> f32 {
    _shannon_entropy_dna(data.as_bytes())
}

pub fn get_snvs_meticulous(
    query_sequence: &[u8],
    query_quals: &[u8],
    ref_seq: &[u8],
    query_coords: &[u64],
    ref_coords: &[u64],
    ref_coord_start: usize,
    tr_start_pos: usize,
    tr_end_pos: usize,
    useful_snvs_params: &UsefulSNVsParams,
) -> HashMap<usize, (char, u8)> {
    let qry_seq_len = query_sequence.len();

    let mut lhs_contiguous: usize = 0;
    let mut rhs_contiguous: usize = 0;
    let mut last_rp: Option<usize> = None;

    let mut snv_group = Vec::<(usize, (char, u8))>::new();

    let mut snvs = HashMap::new();

    let contiguous_threshold = useful_snvs_params.contiguous_threshold;
    let max_snv_group_size = useful_snvs_params.max_snv_group_size;
    let entropy_flank_size = useful_snvs_params.entropy_flank_size;
    let entropy_threshold = useful_snvs_params.entropy_threshold;

    for i in 0..(ref_coords.len()) {
        let ref_pos = ref_coords[i] as usize;

        if tr_start_pos <= ref_pos && ref_pos < tr_end_pos {
            continue;
        }

        let read_pos = query_coords[i] as usize;

        let read_base = query_sequence[read_pos];
        let ref_base = ref_seq[ref_pos - ref_coord_start];

        let contiguous_at_base = match last_rp {
            Some(lr) => useful_snvs_params.contiguous_threshold == 0 || ref_pos - lr == 1,
            None => true,
        };

        if read_base == ref_base && contiguous_at_base {
            let snv_group_len = snv_group.len();

            if snv_group_len > 0 {
                rhs_contiguous += 1;
            } else {
                lhs_contiguous += 1;
            }

            if lhs_contiguous >= contiguous_threshold && rhs_contiguous >= contiguous_threshold {
                if snv_group_len <= max_snv_group_size {
                    for &(snv_pos, snv_a) in snv_group.iter() {
                        snvs.insert(snv_pos, snv_a);
                    }
                }

                // Otherwise, it might be a little mismapped area or a longer deletion vs reference, so ignore it.
                lhs_contiguous = 0;
                rhs_contiguous = 0;
                snv_group.clear();
            }

            last_rp = Some(ref_pos);
            continue;
        }

        if !contiguous_at_base {  // Non-contiguous jump; insertion in query
            lhs_contiguous = 0;
            last_rp = Some(ref_pos);
            continue;
        }

        // Otherwise, contiguous

        if read_base != ref_base {
            // If our entropy is ok, add this to the SNV group
            let seq = &query_sequence[
                read_pos - cmp::min(entropy_flank_size, read_pos)..cmp::min(read_pos + entropy_flank_size, qry_seq_len)
            ];
            if _shannon_entropy_dna(seq) >= entropy_threshold {
                snv_group.push((ref_pos, (read_base as char, query_quals[read_pos])));
            }

            // Don't reset either contiguous variable; instead, take this as part of a SNP group
            last_rp = Some(ref_pos);
        }
    }

    // Special case: if we have stuff in the SNV group with no contiguous requirements,
    // add it to the SNV HashMap.
    let sgl = snv_group.len();
    if contiguous_threshold == 0 && sgl > 0 && sgl <= max_snv_group_size {
        for &(snv_pos, snv_a) in snv_group.iter() {
            snvs.insert(snv_pos, snv_a);
        }
    }

    snvs
}

pub fn get_snvs_simple(
    query_sequence: &[u8],
    query_quals: &[u8],
    ref_seq: &[u8],
    query_coords: &[u64],
    ref_coords: &[u64],
    ref_coord_start: usize,
    tr_start_pos: usize,
    tr_end_pos: usize,
    useful_snvs_params: &UsefulSNVsParams,
) -> HashMap<usize, (char, u8)> {
    let qry_seq_len = query_sequence.len();

    let mut n_snvs = 0;
    let mut res = HashMap::new();

    let too_many_snvs_threshold = useful_snvs_params.too_many_snvs_threshold;
    let entropy_flank_size = useful_snvs_params.entropy_flank_size;
    let entropy_threshold = useful_snvs_params.entropy_threshold;

    for i in 0..query_coords.len() {
        let ref_pos = ref_coords[i] as usize;

        if tr_start_pos <= ref_pos && ref_pos < tr_end_pos {
            continue;
        }

        let read_pos = query_coords[i] as usize;
        let qry_byte_at_pos = query_sequence[read_pos];

        if qry_byte_at_pos == ref_seq[ref_pos - ref_coord_start] {
            continue;
        }

        let seq = &query_sequence[
            read_pos - cmp::min(entropy_flank_size, read_pos)..cmp::min(read_pos + entropy_flank_size, qry_seq_len)
        ];

        if _shannon_entropy_dna(seq) >= entropy_threshold {
            n_snvs += 1;

            // Below it's '>=' but we can skip one set_item by using '>'
            if too_many_snvs_threshold > 0 && n_snvs > too_many_snvs_threshold {
                return res;
            }

            res.insert(ref_pos, (qry_byte_at_pos as char, query_quals[read_pos]));
        }
    }

    res
}

pub trait GetReadSNVs {
    fn get_read_snvs(
        &self,
        ref_seq: &str,
        query_coords: &[u64],
        ref_coords: &[u64],
        ref_coord_start: usize,
        tr_start_pos: usize,
        tr_end_pos: usize,
        useful_snvs_params: &UsefulSNVsParams,
    ) -> HashMap<usize, (char, u8)>;
}

impl GetReadSNVs for STRkitAlignedSegment {
    fn get_read_snvs(
        &self,
        ref_seq: &str,
        query_coords: &[u64],
        ref_coords: &[u64],
        ref_coord_start: usize,
        tr_start_pos: usize,
        tr_end_pos: usize,
        useful_snvs_params: &UsefulSNVsParams,
    ) -> HashMap<usize, (char, u8)> {
        // Given a list of tuples of aligned (read pos, ref pos) pairs, this function finds non-reference SNVs which are
        // surrounded by a stretch of aligned bases of a specified size on either side.
        // Returns a hash map of <position, base>

        let qry_seq_bytes = self.query_sequence.as_bytes();
        let query_quals = self.query_qualities_slice();
        let ref_seq_bytes = ref_seq.as_bytes();

        let snvs = get_snvs_simple(
            qry_seq_bytes,
            query_quals,
            ref_seq_bytes,
            query_coords,
            ref_coords,
            ref_coord_start,
            tr_start_pos,
            tr_end_pos,
            useful_snvs_params,
        );

        if snvs.len() >= useful_snvs_params.too_many_snvs_threshold {  // TOO MANY, some kind of mismapping going on?
            get_snvs_meticulous(
                qry_seq_bytes,
                query_quals,
                ref_seq_bytes,
                query_coords,
                ref_coords,
                ref_coord_start,
                tr_start_pos,
                tr_end_pos,
                useful_snvs_params,
            )
        } else {
            snvs
        }
    }
}

fn find_base_at_pos(
    query_sequence_bytes: &[u8],
    aligned_coords: &STRkitAlignedCoords,
    t: usize,
    start_left: usize,
) -> (char, usize, bool) {
    let (idx, found) = aligned_coords.find_coord_idx_by_ref_pos(t, start_left);

    if found {
        // Even if not in SNV set, it is not guaranteed to be a reference base, since
        // it's possible it was surrounded by too much other variation during the original
        // SNV getter algorithm.
        let qc = query_sequence_bytes[aligned_coords.query_coords[idx] as usize] as char;
        (qc, idx, found)
    } else {
        // Nothing found, so must have been a gap
        (SNV_GAP_CHAR, idx, found)
    }
}

pub fn calculate_useful_snvs(
    py: Python<'_>,
    block_segments: &STRkitLocusBlockSegments,
    read_dict_extra: Bound<'_, PyDict>,
    read_locus_alignment_data: &Bound<'_, PyDict>,
    read_snvs: HashMap<String, HashMap<usize, (char, u8)>>,
    locus_snvs: HashSet<usize>,
    min_allele_reads: usize,
) -> Result<Vec<(usize, usize)>, PyErr> {
    // Mutates read_dict_extra - adds snv_bases keys to read entries

    let n_reads = read_dict_extra.len();

    let mut sorted_snvs: Vec<usize> = locus_snvs.into_iter().collect();
    sorted_snvs.sort();

    let mut snv_counters: HashMap<usize, HashMap<char, usize>> =
        sorted_snvs
            .iter()
            .map(|&s| (s, HashMap::new()))
            .collect();

    for rn in read_dict_extra.keys().into_iter().map(|x| x.cast_into::<PyString>().unwrap()) {
        let rn_str = rn.to_str()?;

        let Some(snvs) = read_snvs.get(rn_str) else {
            continue;
        };

        let segment = block_segments.get_segment_by_name(rn_str)
            .expect("Block segments should include all from read_dict_extra");

        let qs = segment.query_sequence.as_bytes();
        let fqqs = segment.query_qualities_slice();

        let segment_alignment_data_for_locus =
            read_locus_alignment_data
                .get_item(&rn)?
                .unwrap()
                .cast_into::<STRkitSegmentAlignmentDataForLocus>()?
                .borrow();
        let aligned_coords = &segment_alignment_data_for_locus.aligned_coords;

        let segment_start = segment.start as usize;
        let segment_end = segment.end as usize;

        let mut last_pair_idx: usize = 0;
        let snv_list: Vec<(char, u8)> = sorted_snvs.iter().map(|&snv_pos| {
            let mut base: char = SNV_OUT_OF_RANGE_CHAR;
            let mut qual: u8 = 0;

            if segment_start <= snv_pos && snv_pos <= segment_end {
                if let Some(&bb) = snvs.get(&snv_pos) {
                    base = bb.0;
                    qual = bb.1;
                } else {
                    // Binary search for base from correct pair
                    //  - We go in order, so we don't have to search left of the last pair index we tried.
                    let (sb, idx, found) = find_base_at_pos(qs, aligned_coords, snv_pos, last_pair_idx);
                    base = sb;
                    // 0 === indeterminate quality for out-of-range/gaps:
                    qual = if found { fqqs[idx] } else { 0 };

                    //  - Side effect: mutate last_pair_idx
                    last_pair_idx = idx;
                }
            }

            // Otherwise, leave as out-of-range (used immediately below)

            //  - Side effect: mutate snv_counters[snv_pos]
            if base != SNV_OUT_OF_RANGE_CHAR && base != SNV_GAP_CHAR && qual >= 20 {
                // Only count SNV bases where we've actually found the base for the read and it's of a high enough
                // quality to consider.
                let counter = snv_counters.get_mut(&snv_pos).unwrap();
                let count = counter.entry(base).or_insert_with(|| 0);
                *count += 1;
            }

            (base, qual)
        }).collect();

        let read_dict_extra_for_read = read_dict_extra.get_item(&rn)?.unwrap().downcast_into::<PyDict>()?;

        // TODO: set snv_bases as tuple
        read_dict_extra_for_read.set_item(intern!(py, "snv_bases"), snv_list)?;
    }

    // Enough reads to try for SNV based separation

    // require 2 alleles for the SNV, both with at least 20% of the reads, in order to differentiate alleles.
    // require over ~60% of the reads to have the SNV present and at a sufficient quality; otherwise it becomes too
    // easy to occasionally get cases with disjoint sets of SNVs.

    let reads_with_snv_allele_proportion = 0.2;  // TODO: parametrize
    let reads_with_snv_locus_proportion = 0.6;  // TODO: parametrize

    let allele_read_threshold = cmp::max(
        (n_reads as f32 * reads_with_snv_allele_proportion).round() as usize, min_allele_reads);
    let total_read_threshold = cmp::max((n_reads as f32 * reads_with_snv_locus_proportion).round() as usize, 5);
    // snv_counters is guaranteed by the previous inner loop to not have SNV_OUT_OF_RANGE_CHAR or SNV_GAP_CHAR

    let res: Vec<(usize, usize)> = sorted_snvs.iter().enumerate().filter_map(|(s_idx, s_pos)| {
        let snv_counter = snv_counters.get(s_pos).unwrap();

        let n_alleles_meeting_threshold: usize = snv_counter
            .values()
            .map(|&v| (v >= allele_read_threshold) as usize)
            .sum();

        let n_reads_with_this_snv_called: usize = snv_counter.values().sum();

        (
            n_alleles_meeting_threshold >= 2
            && n_reads_with_this_snv_called >= total_read_threshold
        ).then(|| (s_idx, sorted_snvs[s_idx]))
    }).collect();

    Ok(res)
}
