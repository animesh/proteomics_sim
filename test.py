from simulator import Simulator

if __name__ == "__main__":
    N_PEPTIDES = 100
    GRADIENT_MIN = 22
    DIA_WINDOW = 20.0
    TOPN = 10

    print("Proteomics Simulation")
    print(f"Peptides : {N_PEPTIDES}")
    print(f"Gradient : {GRADIENT_MIN} min")
    print(f"DIA window : {DIA_WINDOW} Th")

    results = Simulator().run(
        n_peptides=N_PEPTIDES,
        window=DIA_WINDOW,
        gradient_min=GRADIENT_MIN,
        topn=TOPN
    )

    print("Total peptides:", len(results["peptides"]))
    print("Chromatogram peptides:", len(results["chromatogram"]["peptides"]))
    print("DDA selected peptides:", len(results["dda"]))
    print("DIA windows:", len(results["dia"]))
    print("Library entries:", len(results["library"]))
