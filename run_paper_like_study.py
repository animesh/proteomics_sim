from simulator import ProteomicsStudy
s=ProteomicsStudy()
r=s.run(n_peptides=5000,n_replicates=30,dia_window=2.0)
print("MS1 median CV:", r["cv"]["cv"].median())
print("MS1 mean CV:", r["cv"]["cv"].mean())
print("DDA unique peptides:", r["dda"]["peptide_id"].nunique())
print("DDA coverage (%):", 100 * r["dda"]["peptide_id"].nunique() / r["ms1"]["peptide_id"].nunique())
print("DIA median precursors/window:", r["dia"]["n_precursors"].median())
print("DIA max precursors/window:", r["dia"]["n_precursors"].max())
if "cv_dia" in r:
    print("DIA median CV:", r["cv_dia"]["cv"].median())
