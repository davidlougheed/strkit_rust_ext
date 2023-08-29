from datetime import datetime
from strkit_rust_ext import shannon_entropy


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


if __name__ == "__main__":
    main()
