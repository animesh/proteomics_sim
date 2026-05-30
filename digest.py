import numpy as np, pandas as pd
def generate_peptides(n=5000, gradient_length_sec=1320, seed=42):
    rng = np.random.default_rng(seed)
    return pd.DataFrame({"peptide_id": range(n), "protein_id":
                rng.integers(0, 800, n), "mz":
                rng.uniform(400, 1300, n), "charge":
                rng.choice([2, 3, 4], n, p=[0.6, 0.3, 0.1]), "rt":
                rng.uniform(0.05 * gradient_length_sec, 0.95 * gradient_length_sec, n), "abundance":
                10 ** rng.uniform(3, 7, n), "length":
                rng.integers(7, 35, n), "optimal_ce":
                rng.normal(30, 4, n)})
