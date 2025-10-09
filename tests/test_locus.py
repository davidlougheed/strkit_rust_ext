import pickle
from io import BytesIO
from strkit_rust_ext import STRkitLocus

LOCUS0_ARGS = (0, "locus0", "chr1", 10000, 12000, "CAG", 2, 70)


def _assert_locus0_props_methods(locus: STRkitLocus):
    # test properties
    assert locus.t_idx == 0
    assert locus.locus_id == "locus0"
    assert locus.contig == "chr1"
    assert locus.left_coord == 10000
    assert locus.left_flank_coord == 10000 - 70
    assert locus.right_coord == 12000
    assert locus.right_flank_coord == 12000 + 70
    assert locus.ref_size == 2000
    assert locus.motif == "CAG"
    assert locus.motif_size == 3
    assert locus.n_alleles == 2

    # test methods
    assert locus.log_str() == "locus 0 (id=locus0): chr1:10000-12000 [CAG]"
    assert locus.to_dict() == {
        "locus_index": 0,
        "locus_id": "locus0",
        "contig": "chr1",
        "start": 10000,
        "end": 12000,
        "motif": "CAG",
    }


def test_locus_construction():
    locus = STRkitLocus(*LOCUS0_ARGS)

    # test properties and methods
    _assert_locus0_props_methods(locus)

    # test pickling
    p = pickle.dumps(locus)
    locus_unpickled = pickle.loads(p)
    _assert_locus0_props_methods(locus_unpickled)
