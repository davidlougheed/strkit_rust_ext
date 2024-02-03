from datetime import datetime
from strkit_rust_ext import shannon_entropy, get_snvs_simple, get_snvs_meticulous, get_read_snvs
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
        get_snvs_simple(Q_SEQ, PAIRS, REF_SEQ, 994, 994, 1000, 10, 0.0)
    print(f"get_snvs_simple took {datetime.now() - dt}")

    dt = datetime.now()
    for _ in range(1000000):
        get_snvs_meticulous(Q_SEQ, PAIRS, REF_SEQ, 994, 994, 1000, 0, 5, 10, 0.0)
    print(f"get_snvs_meticulous took {datetime.now() - dt}")

    dt = datetime.now()
    for _ in range(1000000):
        get_read_snvs(Q_SEQ, PAIRS, REF_SEQ, 994, 994, 1000, 0, 5, 20, 10, 0.0)
    print(f"get_read_snvs took {datetime.now() - dt}")


if __name__ == "__main__":
    main()
