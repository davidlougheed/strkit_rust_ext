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

    p = pickle.dumps(ac)
    ac2 = pickle.loads(p)

    assert ac2.query_coord_at_idx(0) == 0
    assert ac2.query_coord_at_idx(1) == 1
    assert ac2.query_coord_at_idx(2) == 2
    assert ac2.query_coord_at_idx(3) == 3
    assert ac2.query_coord_at_idx(4) == 4
    assert ac2.query_coord_at_idx(5) == 5
