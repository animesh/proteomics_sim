
import numpy as np

def assign_chromatography(peptides, rt_max, seed=42):
    rng=np.random.default_rng(seed)
    s=sorted(peptides,key=lambda x:x[3])
    rts=np.linspace(0.05*rt_max,0.95*rt_max,len(s))
    rtmap={p[0]:rt for p,rt in zip(s,rts)}
    abund={p[0]:float(10**rng.uniform(3,7)) for p in s}
    return [(*p,float(rtmap[p[0]]),abund[p[0]]) for p in peptides]
