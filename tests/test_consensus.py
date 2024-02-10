from strkit_rust_ext import best_representatives, consensus_seq


def test_best_representatives():
    assert best_representatives(()) == set()
    assert best_representatives(("AA", "AB", "AA")) == {"AA"}
    assert best_representatives(("AAACAAA", "AACAAA", "AAACAA")) == {"AAACAAA"}


def test_consensus():
    assert consensus_seq(()) is None
    assert consensus_seq(("AA", "AB", "AA")) == "AA"
    assert consensus_seq(("AAACAAA", "AACAAA", "AAACAA")) == "AAACAAA"
    assert consensus_seq(("A", "B", "C")) is None
    assert consensus_seq(("AB", "AAA", "AB", "AB")) == "AB"
