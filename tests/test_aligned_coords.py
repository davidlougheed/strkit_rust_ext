import numpy as np
import pickle
from strkit_rust_ext import STRkitAlignedCoords


def test_aligned_coords():
    qc = np.array([0, 1, 2, 3, 4, 5], dtype=np.uint64)
    rc = np.array([100, 101, 103, 104, 105, 107], dtype=np.uint64)

    ac = STRkitAlignedCoords(qc, rc)

    assert ac.query_coord_at_idx(0) == 0
    assert ac.query_coord_at_idx(1) == 1
    assert ac.query_coord_at_idx(2) == 2
    assert ac.query_coord_at_idx(3) == 3
    assert ac.query_coord_at_idx(4) == 4
    assert ac.query_coord_at_idx(5) == 5

    assert ac.find_coord_idx_by_ref_pos(100, 0) == (0, True)
    assert ac.find_coord_idx_by_ref_pos(101, 0) == (1, True)
    assert ac.find_coord_idx_by_ref_pos(102, 0) == (1, False)
    assert ac.find_coord_idx_by_ref_pos(103, 0) == (2, True)
    assert ac.find_coord_idx_by_ref_pos(104, 0) == (3, True)
    assert ac.find_coord_idx_by_ref_pos(105, 0) == (4, True)
    assert ac.find_coord_idx_by_ref_pos(106, 0) == (4, False)
    assert ac.find_coord_idx_by_ref_pos(107, 0) == (5, True)
    assert ac.find_coord_idx_by_ref_pos(108, 0) == (6, False)  # off the end

    p = pickle.dumps(ac)
    ac2 = pickle.loads(p)

    assert ac2.query_coord_at_idx(0) == 0
    assert ac2.query_coord_at_idx(1) == 1
    assert ac2.query_coord_at_idx(2) == 2
    assert ac2.query_coord_at_idx(3) == 3
    assert ac2.query_coord_at_idx(4) == 4
    assert ac2.query_coord_at_idx(5) == 5
