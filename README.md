*i-At*-LSPHERE collection of genome-scale metabolic models
========================

This repository contains the collection of genome-scale metabolic models reported in Schäfer, Pacheco, *et al.*, as well as the scripts required to generate the models and simulate community outcomes.

# Model collection

Located in the 'Models/Final/' directory, contains genome-scale models for 224 bacterial members of the *Arabidopsis thaliana* phyllosphere microbiome. The models are provided both in .mat format (as outputs of the model generation scripts and for downstream simulation) and in .sbml format for use with COBRApy. The genome files used for generating draft metabolic models are located in 'Models/Genomes.'

# Model generation and simulation scripts

Located in the 'Scripts/' directory, contains scripts for generating models and for carrying out FBA simulations. Please refer to the individual README files in each directory for instructions.

# Carbon source screen results and medium composition

Located in the 'Medium/' directory, contains results of the *in vitro* carbon source screen used to curate the genome-scale models, as well as the minimal medium composition for FBA simulations.