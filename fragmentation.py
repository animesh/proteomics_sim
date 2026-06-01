
import math
from peptides import MASS

def fragments(seq):
    out=[]
    for i in range(1,len(seq)):
        out.append(sum(MASS[a] for a in seq[:i])+1)
        out.append(sum(MASS[a] for a in seq[i:])+19)
    return out

def fragment_efficiency(ce,optimal,width=5):
    return math.exp(-0.5*((ce-optimal)/width)**2)
