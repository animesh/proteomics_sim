from simulator import Simulator

s = Simulator()
r = s.run(n_peptides=5000, window=20, gradient_min=22, topn=10)

print("Total peptides:", len(r["peptides"]))
print("Chromatogram peptides:", len(r["chromatogram"]["peptides"]))
print("DDA selected peptides:", len(r["dda"]))
print("DIA windows:", len(r["dia"]))
print("Library entries:", len(r["library"]))
