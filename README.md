# 4775-final-project
RNA Secondary Structure Prediction

This project implements and compares two RNA secondary structure prediction algorithms: the Nussinov algorithm and the Zuker algorithm. The pipeline takes RNA sequences from a CSV file, runs both algorithms, and evaluates their agreement using overlap scoring.

Overview

	•	Nussinov Algorithm: Predicts RNA secondary structures by maximizing base pairs using dynamic programming.
	•	Zuker Algorithm: Predicts RNA structures by minimizing free energy, using the ViennaRNA library.
	•	Comparison: The predicted structures are compared based on overlap (fraction of shared base pairs).

## Requirements

To ensure all notebooks and scripts in this repository work as expected, please install the following Python packages:

- `pandas`: For data manipulation and analysis.
- `matplotlib`: For creating visualizations and plots.
- `seaborn`: For enhanced data visualization.
- `numpy`: For numerical computations.
- `biopython`: For handling biological sequences and related utilities.
- `RNA` (from ViennaRNA): For RNA secondary structure prediction and stochastic traceback.

### Installation

You can install these packages using `pip`:

```bash
pip install pandas matplotlib seaborn numpy biopython