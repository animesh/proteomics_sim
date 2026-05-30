
import numpy as np
import pandas as pd

try:
    from .chromatography import peak
    from .fragmentation import fragment_efficiency
except ImportError:
    from chromatography import peak
    from fragmentation import fragment_efficiency

def acquire(peptides, time_axis, window=2.0):

    out = []

    bins = np.arange(400, 1300 + window, window)
    mz = peptides.mz.values

    dt = time_axis[1] - time_axis[0] if len(time_axis) > 1 else 1.0

    traces = {}

    for r in peptides.itertuples():
        traces[r.peptide_id] = (
            peak(r.rt, r.abundance, time_axis) *
            fragment_efficiency(30, r.optimal_ce, 5)
        )

    for low in bins[:-1]:

        high = low + window

        idx = np.where((mz >= low) & (mz < high))[0]

        if len(idx) == 0:
            continue

        signal = 0.0

        for i in idx:
            r = peptides.iloc[i]
            signal += traces[r.peptide_id].sum() * dt

        out.append(
            (
                low,
                high,
                signal,
                len(idx)
            )
        )

    return pd.DataFrame(
        out,
        columns=[
            "low",
            "high",
            "signal",
            "n_precursors"
        ]
    )
