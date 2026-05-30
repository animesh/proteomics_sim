# proteomics_sim
The package provides the proteomics LC-MS/MS acquisition (DDA/DIA) quantification simulation framework on which those components can be added for exploring the transition between Data-Dependent Acquisition (DDA) and narrow-window Data-Independent Acquisition (nDIA) in modern proteomics. The goal of this project is to provide a transparent, modifiable, and educational environment for studying how acquisition strategy, chromatography, fragmentation, and quantification interact in shotgun proteomics experiments.

Repository:

[proteomics_sim GitHub Repository](https://github.com/animesh/proteomics_sim)

---

# Motivation

Recent advances in high-speed mass spectrometry, particularly Orbitrap Astral technology, have enabled DIA acquisition with extremely narrow precursor isolation windows (1.2-2 Th). These acquisition strategies produce increasingly DDA-like MS/MS spectra while retaining systematic precursor coverage.

This repository was created to investigate questions such as:

* How narrow must DIA windows become before DIA behaves similarly to DDA?
* How do MS1 and MS2 quantification differ?
* How does dynamic exclusion affect identification coverage?
* How do chromatographic peak widths influence quantitative precision?
* How does collision-energy optimization affect DDA and DIA differently?
* How do completeness, coefficient of variation (CV), and ratio accuracy change under different acquisition schemes?
* What aspects of DDA and DIA are fundamentally different and what aspects are converging?

---

# Scientific Background

This project is inspired by:

**Naomi O'Sullivan,  Florian P Bayer,  Carolin Mogler,  Bernhard Kuster.**

*High-Speed Mass Spectrometers diminish the difference between Data-Dependent and Data-Independent Acquisition Proteomics.*

The study demonstrated that ultra-fast DIA acquisition using narrow precursor windows can achieve comprehensive proteome coverage while reducing spectral complexity and increasing similarity to DDA spectra.

Relevant links:

* https://www.biorxiv.org/content/10.64898/2026.05.26.727836v1 (https://doi.org/10.64898/2026.05.26.727836)
* Data Availability section says that it should be there on PRIDE but i could not find it https://www.ebi.ac.uk/pride/archive?keyword=High-Speed%20Mass%20Spectrometers%20diminish%20the%20difference%20between%20Data%20Dependent%20and%20Data-Independent%20Acquisition%20Proteomics%20

The simulator attempts to reproduce the qualitative mechanisms discussed in the manuscript:

* DDA Top-N precursor selection
* Dynamic exclusion
* Narrow DIA windows
* Peptide-specific fragmentation efficiency
* Fixed versus adaptive collision energy
* Chromatographic variability
* Run-to-run variability
* MS1 quantification
* MS2 quantification
* Quantitative precision and completeness

The simulator is inspired by the paper but does not reproduce instrument firmware, Orbitrap Astral hardware behavior, DIA-NN processing, or identification algorithms used in the publication.

---

# Current Status

⚠️ Research prototype

Current implementation includes:

* Synthetic peptide generation
* Protein assignment
* Chromatographic peak simulation
* Retention-time drift
* Ionization variability
* DDA acquisition
* nDIA acquisition
* Dynamic exclusion
* Collision-energy modeling
* MS1 quantification
* MS2 quantification
* CV estimation
* Data completeness estimation

Current implementation does not include:

* Real peptide identification
* FDR estimation
* DIA deconvolution
* Protein inference
* PTM localization
* Match-between-runs
* Spectral libraries
* Isotope envelopes
* Orbitrap transient simulation
* Astral detector simulation

---

# Repository Structure

```text
proteomics_sim/
├── __init__.py
├── digest.py
├── chromatography.py
├── fragmentation.py
├── dda.py
├── ndia.py
├── quantification.py
├── simulator.py
├── run_paper_like_study.py
├── example.py
└── check.py
```

## digest.py

Synthetic proteome generation.

Generates:

* proteins
* peptides
* charge states
* abundances
* precursor m/z values
* retention times

---

## chromatography.py

Chromatographic simulation.

Provides:

* Gaussian peak model
* peak integration
* retention-time behavior
* peak width interpreted as chromatographic full width at half maximum (FWHM)

---

## fragmentation.py

Fragmentation efficiency model.

Provides:

* peptide-specific optimal collision energies
* fragmentation efficiency penalties

---

## dda.py

Data-dependent acquisition.

Provides:

* Top-N selection
* dynamic exclusion
* charge-dependent collision energy

---

## ndia.py

Narrow-window DIA acquisition.

Provides:

* fixed-width windows
* systematic precursor coverage
* fixed collision energy
* per-window precursor density as `n_precursors`

---

## quantification.py

Quantitative analysis.

Provides:

* coefficient of variation calculations
* peptide-level summaries

---

## simulator.py

Main orchestration layer.

Combines:

* peptide generation
* chromatography
* acquisition
* quantification

into a complete simulated experiment.

---

## Example plots

The repository includes `example.py`, which runs a small simulation and saves plot files into a `plots/` folder.

The saved filenames are generated from the simulation parameters, for example:

* `plots/summary_n1000_rep10_g22_w2.0.png`
* `plots/chromatogram_stacks_n1000_rep10_g22_w2.0.png`

The summary plot shows:

* MS1 coefficient-of-variation distribution
* DDA selected peptide counts per peptide
* average DIA window signal
* DIA window CV distribution when available

The stacked chromatogram plot shows:

* left panel: per-replicate DDA chromatograms with black tick marks indicating picked DDA acquisition times
* right panel: a per-replicate DIA window heatmap showing the number of precursors in each m/z window

Each replicate is labeled by its replicate number, so the stacked view is interpreted as replicate index rather than raw intensity units.

---

# Installation

Clone repository:

```bash
git clone https://github.com/animesh/proteomics_sim.git

cd proteomics_sim
```

Create environment:

```bash
python -m venv venv
source venv/bin/activate
```

Windows:

```powershell
venv\Scripts\activate
```

Install dependencies:

```bash
pip install numpy pandas matplotlib
```

Optional plotting support is provided by `matplotlib` for `example.py`.

---

# Default Parameters

The current implementation uses the following defaults.

## Peptide Population

```python
n_peptides = 5000
n_proteins = 800
```

m/z range:

```python
400 <= mz <= 1300
```

Charge states:

```python
[2,3,4]
```

Probabilities:

```python
[0.6,0.3,0.1]
```

Peptide length:

```python
7-35 amino acids
```

---

# Usage

Run the simulator directly from the package folder:

```powershell
cd \Download\proteomics_sim
python run_paper_like_study.py
```

Or invoke the main entrypoint interactively:

```powershell
cd \Download\proteomics_sim
python -c "from simulator import ProteomicsStudy; r=ProteomicsStudy(seed=42).run(); print(r['cv'].head())"
```

If you want to import using the package name, run Python from the parent folder:

```powershell
cd \Download
python -c "from proteomics_sim import ProteomicsStudy; r=ProteomicsStudy(seed=42).run(); print(r['cv'].head())"
```

Note: this repository does not currently include packaging metadata like `pyproject.toml` or `setup.py`, so editable install is not available by default.

# Testing

A lightweight smoke test is included in `check.py`.

```powershell
cd \Download\proteomics_sim
python check.py
```

A full simulation example is available in `example.py`.

```powershell
cd \Download\proteomics_sim
python example.py
```

The smoke test validates core functionality for:

* `digest.generate_peptides`
* `chromatography.peak`
* `fragmentation.fragment_efficiency`
* `dda.acquire`
* `ndia.acquire`
* `quantification.cv_table`
* `ProteomicsStudy.run`

The example script exercises the full reporting workflow, generates plots, and prints summary statistics for MS1, DDA, and DIA results.

The simulator now also computes DIA window CVs in `results["cv_dia"]`.

If all checks pass, the check script prints `All tests passed.` and exits normally.

Abundance distribution:

```python
10**Uniform(3,7)
```

Equivalent range:

```text
1,000 - 10,000,000 arbitrary units
```

---

## Chromatography

Gradient:

```python
22 minutes
```

Time sampling:

```python
0.6 sec
```

Peak shape:

```python
Gaussian
```

Peak width:

```python
8 sec
```

Retention-time drift:

```python
Normal(0,0.4 sec)
```

---

## Ionization Variability

Run-to-run ionization variability:

```python
LogNormal(
    mean=0,
    sigma=0.15
)
```

Purpose:

* emulate electrospray variation
* generate realistic CV distributions

---

## DDA Defaults

TopN:

```python
10
```

Isolation window:

```python
1.2 Th
```

Dynamic exclusion:

```python
5 sec
```

Collision energy:

```python
CE = 27 + 3 * charge
```

Examples:

| Charge | CE |
| ------ | -- |
| 2+     | 33 |
| 3+     | 36 |
| 4+     | 39 |

---

## nDIA Defaults

Window width:

```python
2.0 Th
```

Collision energy:

```python
30
```

Recommended comparisons:

```python
1.2 Th
2.0 Th
4.0 Th
```

These correspond to the principal window sizes discussed in the narrow-window DIA publication.

---

## Fragmentation Efficiency

Current model:

```python
efficiency =
exp(
-(CE-optimal_CE)^2 /
(2*sigma^2)
)
```

Default:

```python
sigma = 4
```

Purpose:

* simulate peptide-specific fragmentation preferences
* generate DDA-only and DIA-only identifications

---

## Quantification

Current implementation:

MS1:

```python
chromatographic peak integration
```

DDA MS2:

```python
summed fragment intensity
```

DIA MS2:

```python
summed window intensity
```

No isotope envelopes are currently modeled.

---

# Quick Start

```python
from proteomics_sim import ProteomicsStudy

study = ProteomicsStudy(seed=42)

results = study.run(
    n_peptides=5000,
    n_replicates=30,
    gradient_min=22,
    dia_window=2.0
)

print(results["cv"].describe())
```

---

# Understanding run_paper_like_study.py

The script provides a minimal example intended to mimic the conceptual comparisons reported in the narrow-window DIA study.

Current source:

```python
from proteomics_sim import ProteomicsStudy

s = ProteomicsStudy()

r = s.run(
    n_peptides=5000,
    n_replicates=30,
    dia_window=2.0
)

print(r["cv"].describe())
```

## What the Script Simulates

### Step 1

Generate a synthetic proteome.

```text
5000 peptides
800 proteins
```

### Step 2

Assign each peptide:

* abundance
* precursor m/z
* retention time
* charge state
* optimal collision energy

### Step 3

Generate 30 independent LC-MS runs.

Each replicate receives:

* retention-time drift
* ionization variability

### Step 4

Simulate DDA acquisition.

```text
TopN = 10
Isolation = 1.2 Th
Dynamic exclusion = 5 sec
```

### Step 5

Simulate nDIA acquisition.

```text
Window width = 2.0 Th
Full precursor coverage
```

### Step 6

Perform quantification.

Compute:

* MS1 area
* DDA MS2 signal
* DIA MS2 signal
* CV

### Step 7

Return summary tables.

---

# Returned Objects

## results["ms1"]

Columns:

```text
peptide_id
protein_id
ms1_area
rep
```

---

## results["dda"]

Columns:

```text
peptide_id
ms2
rep
```

---

## results["dia"]

Columns:

```text
low
high
signal
```

---

## results["cv"]

Columns:

```text
mean
std
cv
```

---

# Reproducing Paper-Like Experiments

## DIA Window Comparison

```python
for w in [1.2,2.0,4.0]:

    study.run(
        dia_window=w
    )
```

Question:

How rapidly does DIA become DDA-like as windows narrow?

---

## Gradient Comparison

```python
for g in [5,15,30,60]:

    study.run(
        gradient_min=g
    )
```

Question:

How does chromatographic sampling affect quantitative precision?

---

## Replicate Precision

```python
study.run(
    n_replicates=30
)
```

Question:

What fraction of peptides achieve CV < 20%?

---

## Complexity Scaling

```python
for n in [1000,5000,10000]:

    study.run(
        n_peptides=n
    )
```

Question:

How does acquisition strategy behave as proteome complexity increases?

---

# Typical Use Cases

## Teaching

Demonstrate:

* DDA acquisition
* DIA acquisition
* chromatographic sampling
* MS1 versus MS2 quantification

### Method Development

Explore:

* alternative DIA window schemes
* collision-energy strategies
* scan scheduling approaches

### Benchmarking

Evaluate:

* completeness
* precision
* robustness

---

# Planned Features

* tryptic digestion from FASTA
* isotope envelopes
* b-ion generation
* y-ion generation
* fragment chromatograms
* DIA deconvolution
* target-decoy FDR
* protein inference
* MaxLFQ rollup
* phosphoproteomics
* PTM localization
* match-between-runs
* DIA-NN-style processing
* Orbitrap Astral scheduling models

Optional:

```bash
pip install pyarrow numba jupyter scipy matplotlib
```

## Research Planning

Perform in silico experiments before instrument acquisition.

## Software Development

Generate synthetic datasets for:

* search engines
* quantification tools
* visualization software

---

# Scientific Caveats

This repository is a mechanistic simulator. Major simplifications include:

* Gaussian peaks
* simplified fragmentation
* no isotope distributions
* no detector physics
* no identification engine
* no FDR estimation
* no protein grouping

---

# Contributing

Contributions are welcome.

Areas of particular interest:

* realistic fragmentation models
* phosphoproteomics support
* DIA deconvolution
* protein inference
* benchmarking datasets
* visualization modules

