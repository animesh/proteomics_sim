
import numpy as np
import pandas as pd

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

    def run(self, n_peptides=5000, n_replicates=3,
            gradient_min=22, dia_window=2.0):

        gradient_length_sec = gradient_min * 60
        print(f"Starting simulation: n_peptides={n_peptides}, n_replicates={n_replicates}, gradient_min={gradient_min}, dia_window={dia_window}")
        print(f"Gradient length: {gradient_length_sec} sec, time sampling 0.6 sec")

        pep = generate_peptides(
            n=n_peptides,
            gradient_length_sec=gradient_length_sec,
            seed=self.seed
        )

        print(f"Generated peptide library: {len(pep)} peptides")

        t = np.arange(0, gradient_length_sec, 0.6)

        ms1_all = []
        dda_all = []
        dia_all = []

        rng = np.random.default_rng(self.seed)

        for rep in range(n_replicates):
            print(f"\n=== Replicate {rep+1}/{n_replicates} ===")
            print("Applying RT drift and abundance variability")

            pep_rep = pep.copy()

            pep_rep["rt"] += rng.normal(0, 0.4, len(pep_rep))
            pep_rep["abundance"] *= rng.lognormal(0, 0.15, len(pep_rep))

            print("Simulating MS1 chromatograms and peak integration")
            rows = []

            for r in pep_rep.itertuples():
                trace = peak(r.rt, r.abundance, t)
                area = trapz(trace, x=t)

                rows.append(
                    (r.peptide_id, r.protein_id, area, rep)
                )

            ms1 = pd.DataFrame(
                rows,
                columns=[
                    "peptide_id",
                    "protein_id",
                    "ms1_area",
                    "rep"
                ]
            )

            print("Running DDA acquisition")
            dda = dda_acquire(pep_rep, t)
            dda["rep"] = rep

            print("Running DIA window acquisition")
            dia = dia_acquire(pep_rep, t, dia_window)
            dia["rep"] = rep

            ms1_all.append(ms1)
            dda_all.append(dda)
            dia_all.append(dia)

        ms1_all = pd.concat(ms1_all, ignore_index=True)
        dda_all = pd.concat(dda_all, ignore_index=True)
        dia_all = pd.concat(dia_all, ignore_index=True)
        print(f"Finished simulation: total MS1 rows={len(ms1_all)}, DDA rows={len(dda_all)}, DIA rows={len(dia_all)}")

        return {
            "ms1": ms1_all,
            "dda": dda_all,
            "dia": dia_all,
            "cv": cv_table(ms1_all, "ms1_area"),
            "completeness": float(ms1_all["ms1_area"].notna().mean()*100),
            "identified_peptides": int(dda_all["peptide_id"].nunique()) if len(dda_all)>0 else 0,
            "cv_dia": cv_table(
                dia_all,
                "signal",
                groupby=["low", "high"]
            )
        }
