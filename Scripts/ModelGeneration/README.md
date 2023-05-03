*At*-LSPHERE genome-scale metabolic model generation pipeline
========================

This collection of scripts will output a set of curated metabolic models based on organism genomes and experimental information. It is divided into four subsections: (1) generation of draft models using CarveMe (Machado *et al.*, 2018), (2) initial gapfilling of the draft models using NICEgame (Vayena *et al.*, 2022), (3) Additional gapfilling of the models to resolve false positive and negative reactions, and (4) final model formatting and annotation, followed by verification using MEMOTE (Lieven *et al.*, 2020). This guide is based on a recommended folder structure for storing models and reports.

# Local quickstart

Software requirements:
  * [MATLAB](https://www.mathworks.com/products/matlab.html) R2021a or higher
  * [CarveMe](https://carveme.readthedocs.io/en/latest/installation.html)
  * [Python](https://www.python.org) 3.6 or 3.7
  * [COBRA Toolbox](https://opencobra.github.io/cobratoolbox/stable/) v2.24.3 or higher
  * [IBM CPLEX Solver](https://www.ibm.com/products/ilog-cplex-optimization-studio/cplex-optimizer) v12.10
  * NICEgame (from this repository)
  * [MEMOTE](https://memote.readthedocs.io/en/latest/)

## Generate draft metabolic reconstructions using CarveMe:

**Procedure:**
1. Download all desired genomes (in this repo, these are in 'Models/Genomes/'):

2. Using a command line interface, navigate to the CarveMe installation directory and initialize the software:
  ```bash
  $ python3 /Applications/carveme-master/carveme/__init__.py
  ```

3. To generate models for all genomes in a directory, navigate to the directory in which the genomes are stored (i.e., 'Models/Genomes/') and run:
```bash
{
for infile in *.faa.zip; do
   outfile=$(echo $infile | awk -F'[.]' '{print $1}')
   carve $infile -o "../CarveMe/sbml_noGF/$outfile.xml
done
}
```
This will create one SBML draft model corresponding to each genome, and will store them in the 'sbml_noGF' directory.

 Alternatively, to generate models for individual genomes, navigate to desired directory and run:
```bash
{
carve --refseq GCF_XXXXXXXXX.1 -o ../CarveMe/sbml_noGF/GCF_XXXXXXXXX.xml
done
}
```

**Key outputs:**
  * One draft genome-scale model (in SBML format) for each input genome

## Generate gapfilled models using NICEgame:

**Main script:**
* Gapfilling/NICEgame/gapFillModelTFA.m

**Key inputs:**
  * Draft models (in 'FBA/Models/CarveMe/sbml_noGF/')
  * Carbon source screen data ('Medium/CSourceScreen_Jul2022.xlsx')

**Procedure:**
1. Unpack the matTFA toolbox located in NICEgame/matTFA-master/matTFA.zip

2. Open MATLAB and the 'gapFillModelTFA.m' script. This script generates genome-scale metabolic models from previously-generated CarveMe reconstructions and experimental data using the matTFA (Thermodynamic Flux Analysis, Salvy *et al.*, 2019) and NICEgame (Vayena *et al.*, 2022) pipelines.

     This script takes a CarveMe draft metabolic model of an organism and its corresponding experimental data (in .xlsx format representing growth/no growth on carbon sources) as its main inputs. It performs gapfilling using NICEgame and matTFA, which merge the corresponding draft model with a universal metabolite/reaction database and constrains reactions using thermodynamic information. NICEgame then finds candidate reactions that need to be added to the reconstructions to enable growth on each carbon source.

     The script then selects the best combination of gapfilled reactions to use by predicting the growth/no growth phenotype of each model on combinations of solutions. It then saves COBRA model files for downstream curation.

**Key outputs:**
  * List of candidate reactions for gapfilling (in 'FBA/Models/NICEgame/GapfillingResults/')
  * Gapfilled models (one per organism, in 'FBA/Models/NICEgame/Gapfilled/')

## Perform additional model curation to resolve false negative and positive growth:

**Main scripts:**
  * Gapfilling/getModelAccuracy.m
  * Gapfilling/troubleshootFalsePosNeg.m

**Key inputs:**
  * Gapfilled models (one per organism, in 'FBA/Models/NICEgame/Gapfilled/')
  * Carbon source screen data ('Medium/CSourceScreen_Jul2022.xlsx')

**Procedure:**
1. Run the 'getModelAccuracy.m' script, which will output a .mat file containing accuracy statistics of all models in the relevant directory.

2. Run the 'troubleshootFalsePosNeg.m' script, which will reference other models within the collection to correct for false negative and positive growth predictions. Here, the threshold for false positives and the method of correction can be adjusted.

**Key outputs:**
  * FP/FN-corrected gapfilled models (one per organism, in 'FBA/Models/NICEgame/Gapfilled/FPFNCorrected/')

## Perform final model formatting:

**Main scripts:**
  * Final/finalModelFormatting.m

**Key inputs:**
  * FP/FN-corrected gapfilled models (one per organism, in 'FBA/Models/NICEgame/Gapfilled/FPFNCorrected/')
  * Annotation databases (in 'FBA/Scripts/ModelGeneration/Final/databases/')

**Procedure:**
1. Run the 'finalModelFormatting.m' script, which will attempt to annotate all model metabolites, genes, reactions, and subsystems. It will output a .mat file containing the formatted model in COBRA format, as well as an SBML model in .xml.

**Key outputs:**
  * Annotated models in .mat format (one per organism, in 'FBA/Models/Final/')
  * Annotated models in SBML format (one per organism, in 'FBA/Models/Final/sbml')

## Verify models using MEMOTE:

**Key inputs:**
  * Annotated models in SBML format (in 'FBA/Models/Final/sbml')

**Procedure:**

1. Navigate to the directory containing the gapfilled models in SBML format and run MEMOTE via a command line interface to verify the models:

```bash
for i in *.xml; do
  memote report snapshot --filename "../../Reports/${i%.*}.html" "$i" || break
done
```

**Key outputs:**
  * MEMOTE quality scores for each model (in 'FBA/Models/Reports/')