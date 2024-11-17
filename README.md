# 4775-final-project
RNA Secondary Structure Prediction

This project implements and compares two RNA secondary structure prediction algorithms: the Nussinov algorithm and the Zuker algorithm. The pipeline takes RNA sequences from a CSV file, runs both algorithms, and evaluates their agreement using overlap scoring.

Overview

	•	Nussinov Algorithm: Predicts RNA secondary structures by maximizing base pairs using dynamic programming.
	•	Zuker Algorithm: Predicts RNA structures by minimizing free energy, using the ViennaRNA library.
	•	Comparison: The predicted structures are compared based on overlap (fraction of shared base pairs).