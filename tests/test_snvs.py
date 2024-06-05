from strkit_rust_ext import get_read_snvs
from .common import REF_SEQ, Q_SEQ, Q_QUALS, ALIGN_COORDS_Q, ALIGN_COORDS_R, SNVS


def test_get_read_snvs():
    r_items = tuple(
        sorted(
            get_read_snvs(
                Q_SEQ, Q_QUALS, REF_SEQ, ALIGN_COORDS_Q, ALIGN_COORDS_R, 994, 994, 1000, 0, 5, 20, 10, 0.0
            ).items()
        )
    )
    assert SNVS == r_items
