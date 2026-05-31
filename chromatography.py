
import numpy as np

def peak(rt, abundance, t, width=8):
    sigma = width / 2.355  # width is FWHM in seconds
    return abundance * np.exp(-((t - rt) ** 2) / (2 * sigma ** 2))


def assign_chromatography(peptides, rt_min=0.0, rt_max=1320.0, intensity=1.0, seed=None):
    """Assign retention times by increasing precursor m/z and random peptide abundances.

    Returns a new peptide list where each peptide tuple is extended with
    (rt, abundance) at the end.

    MS1 abundance is assigned as a uniform random value between 0 and 1
    (optionally scaled by `intensity`).
    """
    if not peptides:
        return []

    sorted_peptides = sorted(peptides, key=lambda p: p[3])
    n = len(sorted_peptides)

    if n == 1:
        rts = np.array([0.5 * (rt_min + rt_max)])
    else:
        rts = np.linspace(rt_min, rt_max, n)

    rt_by_id = {p[0]: float(rt) for p, rt in zip(sorted_peptides, rts)}
    rng = np.random.default_rng(seed)
    abundances = rng.random(n) * float(intensity)
    abundance_by_id = {p[0]: round(float(abundance), 4) for p, abundance in zip(sorted_peptides, abundances)}

    return [
        (*p, rt_by_id[p[0]], abundance_by_id[p[0]])
        for p in peptides
    ]
