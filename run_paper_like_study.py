from simulator import ProteomicsStudy
s=ProteomicsStudy()
r=s.run(n_peptides=5000,n_replicates=30,dia_window=2.0)
print(r["cv"].describe())
