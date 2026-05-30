import pandas as pd, numpy as np
try:
    from .chromatography import peak
    from .fragmentation import fragment_efficiency
except ImportError:
    from chromatography import peak
    from fragmentation import fragment_efficiency
def acquire(peptides,time_axis,window=2.0):
    out=[]
    bins=np.arange(400,1300+window,window)
    mz=peptides.mz.values
    for t in time_axis:
      for low in bins[:-1]:
        high=low+window
        idx=np.where((mz>=low)&(mz<high))[0]
        if len(idx)==0: continue
        s=0.0
        for i in idx:
          r=peptides.iloc[i]
          s+=peak(r.rt,r.abundance,t)*fragment_efficiency(30,r.optimal_ce,5)
        out.append((low,high,s))
    return pd.DataFrame(out,columns=["low","high","signal"])
