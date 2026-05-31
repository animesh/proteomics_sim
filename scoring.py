import math

def dot_product(a, b):
    keys = set(a) | set(b)
    num = sum(a.get(k, 0) * b.get(k, 0) for k in keys)
    da = math.sqrt(sum(vv * vv for vv in a.values()))
    db = math.sqrt(sum(vv * vv for vv in b.values()))
    return 0 if da == 0 or db == 0 else num / (da * db)


def score_against_library(vector, library_vectors):
    return {pid: dot_product(vector, lib_vector) for pid, lib_vector in library_vectors.items()}
