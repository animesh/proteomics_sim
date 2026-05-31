import math
from peptides import MASS

def fragments(seq):
    b=[]; y=[]
    for i in range(1, len(seq)):
        b.append(sum(MASS[a] for a in seq[:i]) + 1)
        y.append(sum(MASS[a] for a in seq[i:]) + 19)
    return b + y


def fragment_efficiency(ce, optimal_ce, width=5):
    """Simple Gaussian efficiency model for collision energy matching."""
    diff = float(ce) - float(optimal_ce)
    return float(math.exp(-0.5 * (diff / width) ** 2))
