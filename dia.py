
from chromatography import peak

def acquire(peptides,time_axis,window=20.0):
    windows = {}

    for p in peptides:
        low = int(p[3] // window) * window
        trace = []

        for t in time_axis:
            inten = peak(p[4], p[5], t)
            if inten > 1:
                trace.append((t, float(inten)))

        if not trace:
            continue

        windows.setdefault(low, []).append({
            'peptide': p,
            'trace': trace,
            'area': sum(i for _, i in trace)
        })

    return windows
