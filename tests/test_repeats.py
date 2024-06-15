from strkit_rust_ext import get_repeat_count


def test_get_repeat_count():
    assert get_repeat_count(5, "CACACACACA", "ATGC", "ATGC", "CA", 50, 3, 1) == ((5, 36), 5, 0.0)
