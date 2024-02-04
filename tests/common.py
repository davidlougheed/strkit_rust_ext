REF_SEQ = "ACACACATGGCCATAC"
Q_SEQ = "ATCCCAAA"

#         A  A       T  T       C  G       C  C       C  C       A  A       A  A       A  C
PAIRS = [(0, 1000), (1, 1001), (2, 1003), (3, 1004), (4, 1005), (5, 1006), (6, 1008), (7, 1009)]
SNVS = ((1003, "C"), (1009, "A"))

ALIGN_COORDS_Q = list(map(lambda x: x[0], PAIRS))
ALIGN_COORDS_R = list(map(lambda x: x[1], PAIRS))

SNV_CATALOG = [
    (1003, "snp1003", "G", ["C"]),
    (1009, "snp1009", "C", ["A"]),
]

CIGAR_OPS = [(5, 50), (0, 560), (1, 50), (0, 2400), (1, 1)]
