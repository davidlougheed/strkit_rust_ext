import logging
import pytest
from strkit_rust_ext import consensus_seq

logger = logging.getLogger(__name__)


def test_best_representatives():
    assert best_representatives(()) == set()
    assert best_representative(()) == None

    assert best_representatives(("AA", "AB", "AA")) == {"AA"}
    assert best_representative(("AA", "AB", "AA")) == "AA"

    assert best_representatives(("AAACAAA", "AACAAA", "AAACAA")) == {"AAACAAA"}
    assert best_representative(("AAACAAA", "AACAAA", "AAACAA")) == "AAACAAA"


test_params = [
    ((), None),
    (("", "", "", "", ""), ("", "single")),
    (("AA", "AB", "AA"), ("AA", "poa")),
    (("", "AA", "AA"), ("AA", )),
    (("", "A", "A", "A", "A"), ("A", "single")),
    (("", "A", "A", "A", "B"), ("A", "best_rep")),
]


@pytest.mark.parametrize(("seqs", "res"))
def test_consensus(seqs, res):
    assert consensus_seq(seqs, logger, 100) == res
