import pandas as pd
try:
    from .chromatography import peak
    from .fragmentation import fragment_efficiency
except ImportError:
    from chromatography import peak
    from fragmentation import fragment_efficiency
def acquire(peptides,time_axis,topn=10,dynamic_exclusion=5):
    exclusion={}; rows=[]
    for t in time_axis:
        active=[]
        for r in peptides.itertuples():
            inten=peak(r.rt,r.abundance,t)
            if inten>100: active.append((r,inten))
        active.sort(key=lambda x:x[1],reverse=True)
        sel=[]
        for r,inten in active:
            if r.peptide_id in exclusion and t<exclusion[r.peptide_id]: continue
            exclusion[r.peptide_id]=t+dynamic_exclusion
            sel.append((r,inten))
            if len(sel)>=topn: break
        for r,inten in sel:
            rows.append((r.peptide_id,inten*fragment_efficiency(27+3*r.charge,r.optimal_ce)))
    return pd.DataFrame(rows,columns=["peptide_id","ms2"]).groupby("peptide_id",as_index=False).sum()
