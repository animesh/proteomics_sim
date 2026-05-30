import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
from digest import generate_peptides
from chromatography import peak
from fragmentation import fragment_efficiency
from dda import acquire as dda_acquire
from ndia import acquire as dia_acquire
from quantification import cv_table
from simulator import ProteomicsStudy
def assert_true(condition, message):
    if not condition:
        raise AssertionError(message)
def test_digest():
    peptides = generate_peptides(n=50, seed=1)
    assert_true(len(peptides) == 50, "digest.generate_peptides returned incorrect row count")
    expected_columns = {"peptide_id", "protein_id", "mz", "charge", "rt", "abundance", "length", "optimal_ce",}
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
    peptides = generate_peptides(n=30, seed=2)
    time_axis = np.linspace(30, 90, 7)
    dda_table = dda_acquire(peptides, time_axis, topn=5, dynamic_exclusion=5)
    assert_true("peptide_id" in dda_table.columns and "ms2" in dda_table.columns, "dda.acquire returned wrong columns")
    assert_true((dda_table["ms2"] >= 0).all(), "dda.acquire produced negative MS2 intensities")
def test_ndia():
    peptides = generate_peptides(n=20, seed=3)
    time_axis = np.linspace(30, 90, 5)
    dia_table = dia_acquire(peptides, time_axis, window=5.0)
    assert_true(set(dia_table.columns) == {"low", "high", "signal"}, "ndia.acquire returned wrong columns")
    assert_true((dia_table["signal"] >= 0).all(), "ndia.acquire produced negative signal values")
def test_quantification():
    df = pd.DataFrame({"peptide_id": [0, 0, 1, 1], "ms1_area": [100.0, 120.0, 50.0, 60.0],})
    cv = cv_table(df, "ms1_area")
    assert_true("cv" in cv.columns, "quantification.cv_table did not compute cv")
    assert_true(cv.loc[0, "cv"] > 0, "quantification.cv_table computed an invalid cv")
def test_simulator():
    study = ProteomicsStudy(seed=42)
    results = study.run(n_peptides=100, n_replicates=2, gradient_min=10, dia_window=2.0)
    assert_true(set(results.keys()) == {"ms1", "dda", "dia", "cv"}, "ProteomicsStudy.run returned unexpected keys")
    assert_true("peptide_id" in results["ms1"].columns, "ProteomicsStudy.run ms1 output missing peptide_id")
    assert_true("ms2" in results["dda"].columns, "ProteomicsStudy.run dda output missing ms2")
    assert_true("signal" in results["dia"].columns, "ProteomicsStudy.run dia output missing signal")
def main():
    tests = [test_digest, test_chromatography, test_fragmentation, test_dda, test_ndia, test_quantification, test_simulator,]
    for test in tests:
        print(f"Running {test.__name__}...", end=" ")
        test()
        print("PASSED")
    print("\nAll tests passed.\n")
if __name__ == "__main__":
    main()
