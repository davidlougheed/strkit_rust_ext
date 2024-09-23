import numpy as np

REF_SEQ = "ACACACATGGCCATAC"
Q_SEQ = "ATCCCAAA"
Q_QUALS = np.fromiter([90] * len(Q_SEQ), dtype=np.uint8)

#         A  A       T  T       C  G       C  C       C  C       A  A       A  A       A  C
PAIRS = [(0, 1000), (1, 1001), (2, 1003), (3, 1004), (4, 1005), (5, 1006), (6, 1008), (7, 1009)]
SNVS = (1003, ("C", 90)), (1009, ("A", 90))

ALIGN_COORDS_Q = np.fromiter(map(lambda x: x[0], PAIRS), dtype=np.uint64)
ALIGN_COORDS_R = np.fromiter(map(lambda x: x[1], PAIRS), dtype=np.uint64)

SNV_CATALOG = [
    (1003, "snp1003", "G", ["C"]),
    (1009, "snp1009", "C", ["A"]),
]

CIGAR_OPS = np.array([[5, 50], [0, 560], [1, 50], [0, 2400], [1, 1]], dtype=np.uint32)
