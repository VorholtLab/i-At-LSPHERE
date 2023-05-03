*At*-LSPHERE genome-scale metabolic model simulation scripts 
========================

These scripts will simulate competitive outcomes between previously-generated genome-scale models, and will compare these outcomes to experimental data. This guide is based on a recommended folder structure for storing models, but can be modified in each script.

# Local quickstart

Software requirements:
  * [MATLAB](https://www.mathworks.com/products/matlab.html) R2021a or higher
  * [COBRA Toolbox](https://opencobra.github.io/cobratoolbox/stable/) v2.24.3 or higher
  * [IBM CPLEX Solver](https://www.ibm.com/products/ilog-cplex-optimization-studio/cplex-optimizer) v12.10

## Compute competitive outcomes and compare to experimental data:

**Main script:**
* competitiveOutcomesPairs.m

**Key inputs:**
  * Curated models (in 'Models/Final/')
  * Medium composition ('Medium/minMedCSourceScreen.mat')

**Procedure:**
1. Open MATLAB and the 'competitiveOutcomesPairs.m' script. This script computes competitive outcomes between strain pairs and community compositions, and compares them to experimental outcomes if desired.

**Key outputs:**
  * Pairwise and community competitive outcomes and associated metabolic flux information