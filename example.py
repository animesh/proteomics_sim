from simulator import ProteomicsStudy
import numpy as np
# -----------------------------
# Run simulation
# -----------------------------
N_PEPTIDES = 1000
N_REPLICATES = 30
GRADIENT_MIN = 22
DIA_WINDOW = 2.0
print("Proteomics Simulation")
print(f"Peptides : {N_PEPTIDES}")
print(f"Replicates : {N_REPLICATES}")
print(f"Gradient : {GRADIENT_MIN} min")
print(f"DIA window : {DIA_WINDOW} Th")
results = ProteomicsStudy(seed=42).run(n_peptides=N_PEPTIDES, n_replicates=N_REPLICATES, gradient_min=GRADIENT_MIN, dia_window=DIA_WINDOW,)
# -----------------------------
# Extract tables
# -----------------------------
ms1 = results["ms1"]
dda = results["dda"]
dia = results["dia"]
cv = results["cv"]
# -----------------------------
# Basic statistics
# -----------------------------
print("PEPTIDE POPULATION")
print(f"Total peptides : {len(cv)}")
print(f"Total proteins : {ms1['protein_id'].nunique()}")
# -----------------------------
# Detection statistics
# -----------------------------
valid_cv = cv["cv"].dropna()
detected = len(valid_cv)
print("DETECTION")
print(f"Detected peptides : {detected}")
print(f"Missing peptides : {len(cv) - detected}")
print(f"Detection rate : {100*detected/len(cv):.1f}%")
# -----------------------------
# Quantification precision
# -----------------------------
print("MS1 QUANTIFICATION")
print(f"Median CV : {valid_cv.median():.2f}%")
print(f"Mean CV : {valid_cv.mean():.2f}%")
print(f"Best CV : {valid_cv.min():.2f}%")
print(f"Worst CV : {valid_cv.max():.2f}%")
print(f"CV < 10% : {(valid_cv < 10).mean()*100:.1f}%")
print(f"CV < 20% : {(valid_cv < 20).mean()*100:.1f}%")
print(f"CV < 30% : {(valid_cv < 30).mean()*100:.1f}%")
# -----------------------------
# DDA summary
# -----------------------------
print("DDA ACQUISITION")
dda_ids = dda["peptide_id"].nunique()
print(f"Identified peptides : {dda_ids}")
print(f"Coverage : {100*dda_ids/len(cv):.1f}%")
dda_counts = (dda.groupby("peptide_id").size().sort_values(ascending=False))
print("DDA SELECTION STATISTICS")
print("Median selections :", dda_counts.median())
print("Mean selections   :", dda_counts.mean())
print("Max selections    :", dda_counts.max())
# -----------------------------
# DIA summary
# -----------------------------
print("nDIA ACQUISITION")
print(f"DIA windows : {len(dia)}")
signal_windows = (dia["signal"] > 0).sum()
print(f"Active windows : {signal_windows}")
# -----------------------------
# Intensity distribution
# -----------------------------
print("INTENSITY")
print(f"Median MS1 area : "
    f"{ms1['ms1_area'].median():.3e}")
print(f"Max MS1 area : "
    f"{ms1['ms1_area'].max():.3e}")
print(f"Min non-zero area : "
    f"{ms1.loc[ms1['ms1_area']>0,'ms1_area'].min():.3e}")
# -----------------------------
# Top peptides
# -----------------------------
print("TOP 10 PEPTIDES BY MEAN INTENSITY")
top = (cv
    .sort_values("mean", ascending=False)
    .head(10))
print(top[["mean", "cv"]])
print("END OF REPORT")
