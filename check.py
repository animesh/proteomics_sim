import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
from digest import generate_peptides as digest_generate_peptides
from peptides import generate_peptides as random_generate_peptides
from chromatography import peak
from chrom_assign import assign_chromatography
from fragmentation import fragment_efficiency
from dda import acquire as dda_acquire
from ndia import acquire as dia_acquire
from quantification import cv_table
from simulator import Simulator
def assert_true(condition, message):
    if not condition:
        raise AssertionError(message)
def test_digest():
    peptides = digest_generate_peptides(n=50, seed=1)
    assert_true(len(peptides) == 50, "digest.generate_peptides returned incorrect row count")
    expected_columns = {"peptide_id", "protein_id", "mz", "charge", "rt", "abundance", "length", "optimal_ce"}
    assert_true(set(peptides.columns) == expected_columns, "digest.generate_peptides returned wrong columns")
    assert_true(peptides["mz"].between(400, 1300).all(), "mz values out of expected range")
    assert_true(peptides["charge"].isin([2, 3, 4]).all(), "unexpected charge states")
def test_chromatography():
    t = np.linspace(0, 100, 11)
    values = peak(rt=50, abundance=1e6, t=t, width=8)
    assert_true(values.shape == (11,), "chromatography.peak returned wrong shape")
    assert_true(np.all(values >= 0), "chromatography.peak returned negative intensities")
    assert_true(values.max() == values[5], "chromatography.peak peak center not at expected index")
def test_fragmentation():
    eff_same = fragment_efficiency(30, 30)
    eff_off = fragment_efficiency(35, 30)
    assert_true(0 < eff_same <= 1, "fragment_efficiency should be between 0 and 1")
    assert_true(eff_off < eff_same, "fragment_efficiency should decrease away from optimal CE")
def test_dda():
    peptides = random_generate_peptides(n=30, seed=2)
    peptides = assign_chromatography(peptides, rt_max=90.0)
    time_axis = np.linspace(30, 90, 7)
    selections = dda_acquire(peptides, time_axis, topn=5)
    assert_true(isinstance(selections, list), "dda.acquire returned wrong type")
    assert_true(len(selections) > 0, "dda.acquire returned no selections")
    assert_true(
        all(
            isinstance(obs, tuple)
            and len(obs) == 3
            and isinstance(obs[1], tuple)
            and len(obs[1]) >= 6
            for obs in selections
        ),
        "dda.acquire returned invalid selection tuples"
    )
def test_ndia():
    peptides = digest_generate_peptides(n=20, seed=3)
    time_axis = np.linspace(30, 90, 5)
    dia_table = dia_acquire(peptides, time_axis, window=5.0)
    assert_true({"low", "high", "signal"}.issubset(dia_table.columns), "ndia.acquire returned wrong columns")
    assert_true((dia_table["signal"] >= 0).all(), "ndia.acquire produced negative signal values")
def test_quantification():
    df = pd.DataFrame({"peptide_id": [0, 0, 1, 1], "ms1_area": [100.0, 120.0, 50.0, 60.0],})
    cv = cv_table(df, "ms1_area")
    assert_true("cv" in cv.columns, "quantification.cv_table did not compute cv")
    assert_true(cv.loc[0, "cv"] > 0, "quantification.cv_table computed an invalid cv")
def test_simulator():
    study = Simulator()
    results = study.run(n_peptides=100, window=20, gradient_min=10, topn=5)
    expected = {"library", "peptides", "chromatogram", "ms1", "ms2", "dda", "dia", "metrics"}
    assert_true(set(results.keys()) == expected, "Simulator.run returned unexpected keys")
    assert_true("peptides" in results, "Simulator.run output missing peptides")
    assert_true("chromatogram" in results and "peptides" in results["chromatogram"], "Simulator.run output missing chromatogram peptides")
    assert_true("ms1" in results and "peptides" in results["ms1"], "Simulator.run output missing ms1 peptides")
    assert_true("dda" in results, "Simulator.run output missing dda")
    assert_true("dia" in results, "Simulator.run output missing dia")
    assert_true("library" in results, "Simulator.run output missing library")
    assert_true("ms2" in results and "dda" in results["ms2"] and "dia" in results["ms2"], "Simulator.run output missing ms2 structure")
    assert_true("metrics" in results, "Simulator.run output missing metrics")
def main():
    tests = [test_digest, test_chromatography, test_fragmentation, test_dda, test_ndia, test_quantification, test_simulator,]
    for test in tests:
        print(f"Running {test.__name__}...", end=" ")
        test()
        print("PASSED")
    print("\nAll tests passed.\n")
if __name__ == "__main__":
    main()
