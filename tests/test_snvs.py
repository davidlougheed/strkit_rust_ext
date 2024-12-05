import sys
import time
from strkit_rust_ext import shannon_entropy, get_read_snvs
from .common import REF_SEQ, Q_SEQ, Q_QUALS, ALIGN_COORDS_Q, ALIGN_COORDS_R, SNVS


def test_shannon_entropy():
    assert shannon_entropy(b"AAAA") == 0.0

    iters = 500000
    t = time.perf_counter()
    for i in range(iters):
        shannon_entropy(b"ATGCATGCATGCAAAAATTTTTAATATATGCGCCCCCCATGCATGCATGCAAAAATTTTTAATATATGCGCCCCCCATGCATGCATGCAAAAATTTTTAATATATGCGCCCCCC")
    print(f"{iters} iters took {time.perf_counter() - t}", file=sys.stderr, flush=True)


def test_get_read_snvs():
    r_items = tuple(
        sorted(
            get_read_snvs(
                Q_SEQ, Q_QUALS, REF_SEQ, ALIGN_COORDS_Q, ALIGN_COORDS_R, 994, 994, 1000, 0, 5, 20, 10, 0.0
            ).items()
        )
    )
    assert SNVS == r_items
