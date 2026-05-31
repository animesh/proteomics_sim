
import numpy as np
import pandas as pd
from chromatography import peak
from fragmentation import fragment_efficiency

def acquire(peptides, time_axis, window=2.0, interference_alpha=0.05):
    out=[]
    bins=np.arange(400,1300+window,window)
    mz=peptides.mz.values
    dt=time_axis[1]-time_axis[0] if len(time_axis)>1 else 1.0
    traces={}
    for r in peptides.itertuples():
        traces[r.peptide_id]=peak(r.rt,r.abundance,time_axis)*fragment_efficiency(30,r.optimal_ce,5)
    rng=np.random.default_rng(42)
    for low in bins[:-1]:
        high=low+window
        idx=np.where((mz>=low)&(mz<high))[0]
        if len(idx)==0: continue
        signal=0.0
        n_prec=len(idx)
        interference=1.0+interference_alpha*np.log1p(max(n_prec-1,0))
        for i in idx:
            r=peptides.iloc[i]
            signal += traces[r.peptide_id].sum()*dt
        noisy=signal*(1+rng.normal(0,0.02*interference))
        out.append((low,high,noisy,n_prec,interference))
    return pd.DataFrame(out,columns=["low","high","signal","n_precursors","interference"])
