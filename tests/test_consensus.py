import logging
import pytest
from strkit_rust_ext import consensus_seq

logger = logging.getLogger(__name__)


test_params = [
    ((), None),
    (("", "", "", "", ""), ("", "single")),
    (("", "", "", "", "A"), ("", "best_rep")),
    (("AA", "AB", "AA"), ("AA", "best_rep")),
    (("", "AA", "AA"), ("AA", "best_rep")),
    (("", "AA", "AA", "AB"), ("AA", "best_rep")),
    (("", "AA", "AA", "AB", "AC"), ("AA", "poa")),
    (("A", "A", "A", "A"), ("A", "single")),
    (("", "A", "A", "A", "A"), ("A", "best_rep")),
    (("", "A", "A", "A", "B"), ("A", "best_rep")),
    (("A", "A", "A", "AAA", "AAB"), ("A", "best_rep")),  # should decide on best_rep since median non-blank <= 1
    (("", "A", "A", "AAA", "AAB"), ("AA", "poa")),
    (("AAACAAA", "AACAAA", "AAACAA"), ("AAACAAA", "poa")),
]


@pytest.mark.parametrize(("seqs", "res"), test_params)
def test_consensus(seqs, res):
    assert consensus_seq(seqs, logger, 100) == res
