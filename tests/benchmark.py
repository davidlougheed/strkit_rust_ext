from datetime import datetime
from strkit_rust_ext import shannon_entropy, get_read_snvs
from .common import REF_SEQ, Q_SEQ, PAIRS, SNVS, SNV_CATALOG



strs = [
    b"ATGCGCGATAGAGCTAGTCGATGCC",
    b"ATGCTGATCGATCGGCGCGATATAC",
    b"AAAAAAAATTTTCCCTCTCTGGGAA",
    b"ATATTTTATATTTATTTATATATAT",
]


def main():
    dt = datetime.now()
    for _ in range(100000):
        for s in strs:
            shannon_entropy(s)
    print(f"shannon took {datetime.now() - dt}")

    dt = datetime.now()
    for _ in range(1000000):
        get_read_snvs(Q_SEQ, PAIRS, REF_SEQ, 994, 994, 1000, 0, 5, 20, 10, 0.0)
    print(f"SNVs took {datetime.now() - dt}")


if __name__ == "__main__":
    main()
