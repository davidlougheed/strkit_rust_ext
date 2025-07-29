from strkit_rust_ext import get_repeat_count


def test_get_repeat_count():
    assert get_repeat_count(5, "CACACACACA", "ATGC", "ATGC", "CA", 50, 3, 1) == ((5, 36), 1, 0)
    assert get_repeat_count(5, "CACACACA", "ATGC", "ATGC", "CA", 50, 3, 1) == ((4, 32), 5, -1)
    assert get_repeat_count(4, "CACACCCC", "ATGC", "ATGC", "CA", 50, 3, 1) == ((4, 14), 5, 0)
