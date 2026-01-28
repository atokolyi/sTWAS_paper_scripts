# sTWAS paper scripts
For the manuscript "Biobank-scale Bayesian TWAS reveals splicing-mediated mechanisms of complex disease", available on [medRxiv](https://www.medrxiv.org/content/10.1101/2025.09.19.25336165v2).

Note: Due to the specificities of the UK Biobank and All of Us research environments, this generalized code will need to be adapted to the specific platform/environment of the user.

## 1. System requirements
These scripts were ran on linux with python 3.13.1 and R 4.4.2 though should generalize to any OS.
For TWAS fine-mapping, the package [B-LORE](https://github.com/soedinglab/b-lore) is required to be installed from source.
Large amounts of memory may be required for biobank scale fine-mapping in the per-cohort B-LORE stage (i.e. 400-500GB for UKBB, stage 6), though the remainder of the code should run smoothly on a standard system.

The following packages are required throughout the scripts:
### python
```
pandas_plink
scipy
sklearn
numpy
pandas
statsmodels
```
### R (only necessary for All of Us PheCode creation, steps 4.1 & 4.2)
```
bigrquery
tidyverse
stringr
```

## 2. Installation guide
After installing the above required packages, clone the source code to run the scripts locally.
```
git clone https://github.com/atokolyi/sTWAS_paper_scripts/
```

## 3. Demo
As data from the UKBB and AoU is protected, I have included simulated data to assist in running the B-LORE adaptation for TWAS (steps 6 and 7). After performing the installation instructions above, run the following to perform the Bayesian fine-mapping and meta-analysis for one condition (disease1) on two biobanks (biobank1 and biobank2).
Here as in the paper I am using a ZMAX of 10 (maximum number of causal features per LD block), and a PROB_THRESH of 0.5 (causal probability threshold of the posterior probabilities, only used for script logging).
```
cd sTWAS_paper_scripts
python scripts/6_finemapping_BLORE_per_biobank.py disease1 biobank1
python scripts/6_finemapping_BLORE_per_biobank.py disease1 biobank2
python scripts/7_finemapping_BLORE_meta_analysis.py 10 0.5
```
The expected output is as follows, two predicted causal expression features including one splicing event (splice3) and one gene (gene4), with approximately similar posterior probabilities due to procedural stochasticity.
```
Finished B-LORE meta-analysis for (disease1)
         - 2/20 features with posterior probability > 50.0%:
                 - splice3:     PP=100.0%
                 - gene4:       PP=85.34%
Done!
```
The output of the script is stored in `demo/outputs` and expected output is stored in `demo/expected_output/`.
A summary file is saved at `demo/outputs/disease1_meta_z10_results.tsv` with posterior probabilities for all tested features (10 transcript splicing events and 10 genes in this simulated case).
This demo should run in only a few seconds on a standard desktop computer. 

## 4. Instructions for use
Run each step in `scripts` sequentially using your QTLs, RNA-seq cohort expression and genotype matricies (steps 1-3), and biobank EHR and genotype data (steps 3-7).