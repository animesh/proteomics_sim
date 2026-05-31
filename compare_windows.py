
from simulator import ProteomicsStudy
import pandas as pd
rows=[]
for w in [4.0,2.0,1.2]:
    r=ProteomicsStudy(seed=42).run(n_peptides=5000,n_replicates=3,dia_window=w)
    rows.append({
        "window":w,
        "ms1_cv":r["cv"]["cv"].median(),
        "dia_cv":r["cv_dia"]["cv"].median(),
        "median_precursors":r["dia"]["n_precursors"].median(),
        "max_precursors":r["dia"]["n_precursors"].max(),
        "median_interference":r["dia"]["interference"].median(),
        "completeness":r.get("completeness"),
        "identified_peptides":r.get("identified_peptides")
    })
print(pd.DataFrame(rows))
