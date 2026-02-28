import pickle
from io import BytesIO
from strkit_rust_ext import STRkitLocus

LOCUS0_ARGS = (1, "locus1", "chr1", 10000, 12000, "CAG", 2, 70, ["annot1"])


def _assert_locus0_props_methods(locus: STRkitLocus):
    # test properties
    assert locus.t_idx == 1
    assert locus.locus_id == "locus1"
    assert locus.contig == "chr1"
    assert locus.left_coord == 10000
    assert locus.left_flank_coord == 10000 - 70
    assert locus.right_coord == 12000
    assert locus.right_flank_coord == 12000 + 70
    assert locus.ref_size == 2000
    assert locus.motif == "CAG"
    assert locus.motif_size == 3
    assert locus.n_alleles == 2
    assert locus.annotations == ["annot1"]

    # test methods
    log_str = "locus 1 (id=locus1): chr1:10000-12000 [CAG]"
    assert locus.log_str() == log_str
    assert locus.to_dict() == {
        "locus_index": 1,
        "locus_id": "locus1",
        "contig": "chr1",
        "start": 10000,
        "end": 12000,
        "motif": "CAG",
        "annotations": ["annot1"],
    }
    assert repr(locus) == (
        f"<STRkitLocus t_idx=1 locus_id=locus1 contig=chr1 left_coord=10000 left_flank_coord=9930 right_coord=12000 "
        f"right_flank_coord=12070 ref_size=2000 motif=CAG motif_size=3 n_alleles=2 flank_size=70 _log_str={log_str}>"
    )


def test_locus_construction():
    locus = STRkitLocus(*LOCUS0_ARGS)

    # test we can hash locus
    locus_hash = hash(locus)

    # test we can get a repr of the locus - value will be validated in _assert_locus0_props_methods
    locus_repr = repr(locus)

    # test properties and methods
    _assert_locus0_props_methods(locus)

    # test pickling
    p = pickle.dumps(locus)
    locus_unpickled = pickle.loads(p)
    _assert_locus0_props_methods(locus_unpickled)

    # test reprs are equal
    assert repr(locus_unpickled) == locus_repr

    # test hashes are equal
    assert hash(locus_unpickled) == locus_hash


def test_locus_with_ref_data_construction():
    locus = STRkitLocus(*LOCUS0_ARGS)
    ref_cags = "CAG" * 667
    lfs = "TCAGT" * 14
    rfs = "TGACT" * 14

    lwrd = locus.with_ref_data(9999, 12001, "chr1", 667, ref_cags, lfs, rfs, lfs + ref_cags + rfs, 0.1)

    assert lwrd.left_coord_adj == 9999
    assert lwrd.right_coord_adj == 12001
    assert lwrd.ref_cn == 667
    assert lwrd.ref_seq == ref_cags
    assert lwrd.ref_left_flank_seq == lfs
    assert lwrd.ref_right_flank_seq == rfs
    assert lwrd.ref_total_seq == lfs + ref_cags + rfs
    assert lwrd.ref_time == 0.1

    assert lwrd.locus_def.left_flank_coord == 10000 - 70
    _assert_locus0_props_methods(lwrd.locus_def)
