import numpy as np,pandas as pd
try:
    from .digest import generate_peptides
    from .chromatography import peak
    from .dda import acquire as dda_acquire
    from .ndia import acquire as dia_acquire
    from .quantification import cv_table
except ImportError:
    from digest import generate_peptides
    from chromatography import peak
    from dda import acquire as dda_acquire
    from ndia import acquire as dia_acquire
    from quantification import cv_table
def trapz(y, x=None, dx=1.0):
    y = np.asarray(y)
    if x is None:
        return np.sum((y[:-1] + y[1:]) * 0.5) * dx
    x = np.asarray(x)
    return np.sum((y[:-1] + y[1:]) * 0.5 * np.diff(x))
class ProteomicsStudy:
    def __init__(self, seed=42):
        self.seed = seed
    def run(self, n_peptides=5000, n_replicates=3, gradient_min=22, dia_window=2.0):
        gradient_length_sec = gradient_min * 60
        pep = generate_peptides(n=n_peptides, gradient_length_sec=gradient_length_sec, seed=self.seed)
        t = np.arange(0, gradient_length_sec, 0.6)
        ms1_all = []
        dda_all = []
        rng = np.random.default_rng(self.seed)
        for rep in range(n_replicates):
            rows = []
            for r in pep.itertuples():
                rt = (r.rt +
                    rng.normal(0, 0.4))
                abund = (r.abundance *
                    rng.lognormal(0, 0.15))
                trace = peak(rt, abund, t)
                area = trapz(trace, x=t)
                rows.append((r.peptide_id, r.protein_id, area, rep))
            ms1 = pd.DataFrame(rows, columns=["peptide_id", "protein_id", "ms1_area", "rep"])
            dda = dda_acquire(pep, t)
            dda["rep"] = rep
            ms1_all.append(ms1)
            dda_all.append(dda)
        ms1_all = pd.concat(ms1_all, ignore_index=True)
        dda_all = pd.concat(dda_all, ignore_index=True)
        dia = dia_acquire(pep, t, dia_window)
        return {"ms1": ms1_all, "dda": dda_all, "dia": dia, "cv": cv_table(ms1_all, "ms1_area")}
