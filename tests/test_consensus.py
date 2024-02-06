from strkit_rust_ext import consensus_seq


def test_consensus():
    assert consensus_seq(("AA", "AB", "AA")) == "AA"
    assert consensus_seq(("AAACAAA", "AACAAA", "AAACAA")) == "AAACAAA"
    assert consensus_seq(("A", "B", "C")) is None
