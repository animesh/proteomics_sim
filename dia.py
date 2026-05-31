def acquire(peptides,window=20):
    bins={}
    for p in peptides:
        low=int(p[3]//window)*window
        bins.setdefault(low,[]).append(p)
    return bins
