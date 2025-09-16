import pytest
from strkit_rust_ext import get_repeat_count


@pytest.mark.parametrize("grc_args,result", [
    (
        (5, "CACACACACA", "ATGC", "ATGC", "CA", 50, 3, 1, True),
        ((5, 36), 1, 0),
    ),
    (
        (5, "CACACACA", "ATGC", "ATGC", "CA", 50, 3, 1, True),
        ((4, 32), 5, -1),
    ),
    (
        (4, "CACACCCC", "ATGC", "ATGC", "CA", 50, 3, 1, True),
        ((4, 14), 5, 0),
    ),
    # ----
    (
        (5, "CACACACACA", "ATGC", "ATGC", "CA", 50, 3, 1, False),
        ((5, 36), 1, 0),
    ),
    (
        (5, "CACACACA", "ATGC", "ATGC", "CA", 50, 3, 1, False),
        ((4, 32), 9, -1),
    ),
    (
        (4, "CACACCCC", "ATGC", "ATGC", "CA", 50, 3, 1, False),
        ((4, 14), 9, 0),
    ),
])
def test_get_repeat_count(grc_args, result):
    assert get_repeat_count(*grc_args) == result
