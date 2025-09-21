import pandas as pd
import numpy as np
import sys
import time
import pickle
import gc
import sys
import collections
import re

# Load B-LORE source functions for adaptation (https://github.com/soedinglab/b-lore)
sys.path.append('b-lore/blore')
from inference.logistic_regression import LogisticRegression
from inference.optimize_regularizer import OptimizeRegularizer
from meta import optimize_hyperparameters
from iotools.io_summary import WriteSummary

ft = time.time()

# Load UKBB covariates and PheCodes
cov = pd.read_csv("covs_wbs_sexmatch.tsv",sep="\t")
cov['id'] = cov['id'].astype(str)
cond = pd.read_pickle("phenome.pickle");
cond['IID'] = cond['IID'].astype(str)
cond = cond.set_index('IID').loc[cov['id']].reset_index()
cov.index = cov['id']
cond.index = cond['IID']

# Load the sTWAS results
block_keep = pd.read_csv("stwas_out_filt.tsv",sep="\t",header=0)
phens_to_blore = pd.read_csv("stwas_sig_phecodes.txt",sep="\t",header=0)
CONDS = sorted(list(phens_to_blore['x'])) 

# Parallelise by PheCode
COND=CONDS[int(sys.argv[1])-1]
print("Doing ",sys.argv[1],"/",str(len(CONDS)),sep="")
print(COND)
sys.stdout.flush()

# Subset to valid samples for PheCode and 
valid = cond[COND].notna()
y = cond[COND][valid]

# Load imputed residualized splices, gene expression, and covs
splice = pd.read_pickle("splice_gene_age_sex_df.pkl")
features = list(splice.columns)
splice.drop(splice[~valid].index,inplace=True)
splice = splice.to_numpy()
covs = splice[:, -2:]
splice = splice[:, :-2]
features = features[:-2]

feature_meta = pd.read_csv("ld_block_mapping.tsv",header=0,sep="\t")

# Limit to features imputed in both UKBB and AoU
with open('../../final_ukbb_aou/blore/aou_features.pkl', 'rb') as file: 
    aou_features = pickle.load(file) 

fkeep = [x in aou_features for x in features]
features = np.array(features)[fkeep]
splice = splice[:, fkeep]
feature_meta = feature_meta.set_index('feature', drop=False)
feature_meta = feature_meta.loc[features]
features = list(features)

# Limit to AoU and UKBB sig splice blocks
sig_splices = list(set(list(block_keep[block_keep['cond']==COND]['splice'])))
keep_blocks = list(set(list(feature_meta[feature_meta['feature'].isin(sig_splices)]['ld_block'])))
feature_meta = feature_meta[feature_meta['ld_block'].isin(keep_blocks)]
keep_i = [features.index(x) for x in list(feature_meta['feature'])]
splice = splice[:,keep_i]
features = list(np.array(features)[keep_i])

# Construct matricies for B-LORE
gt = [splice[:,[features.index(x) for x in feature_meta[feature_meta['ld_block']==y]['feature']]] for y in list(set(list(feature_meta['ld_block'])))]
del splice
gc.collect()

phenotype = y
phenotype = tuple(np.array(phenotype))
gt = tuple([x for x in gt])

# Using default optimization parameters
mureg = 0.0
sigreg = 0.01
tolerance = 0.01
cov=covs

print ("Optimizing regularizer ...")
t = time.time()
sigreg_optim = OptimizeRegularizer(gt, phenotype, mureg, sigreg, tolerance, cov)
sigreg_optim.update()
mureg = sigreg_optim.mureg
sigreg = sigreg_optim.sigmareg
niter = sigreg_optim.niter
print("\tTime: ",(time.time()-t)/60," mins")

covinfo = [["age","sex"]]

# Run the Bayesian multiple logistic regression
t = time.time()
print ("Calculating summary statistics")
logreg = LogisticRegression(gt, phenotype, mureg, sigreg, covariates = cov)
logreg.fit()
v0 = logreg.v0
vmin = logreg.vmin
precll = logreg.precll
iscov = logreg.iscov
print("\tLogreg time: ",time.time()-t)

outdir = "ukbb_out/"+COND
file_prefix = COND

# Locus names (LD blocks) and transcriptomic events for output
SNPINFO_FIELDS = ['rsid', 'bp_location', 'ref_allele', 'alt_allele']
class SnpInfo(collections.namedtuple('_SnpInfo', SNPINFO_FIELDS)):
    __slots__ = ()
loci = [[x for x in feature_meta[feature_meta['ld_block']==y]['feature']] for y in list(set(list(feature_meta['ld_block'])))]
snpinfo = list()
for locus in loci:
    snp_locus = list()
    for feature in locus:
        this_snp = SnpInfo(rsid = feature,
                       bp_location = "1",
                       alt_allele = "X",
                       ref_allele = "X")
        snp_locus.append(this_snp)
    snpinfo.append(snp_locus)
locnames = list(set(list(feature_meta['ld_block'])))

# Write out the summary statistics
print ("Saving results")
dupes = [[] for x in range(len(loci))]
freq  = [[1 for x in range(len(loci[y]))] for y in range(len(loci))]
summary = WriteSummary(outdir, file_prefix)
summary.set_statistics(snpinfo, covinfo, dupes, freq, v0, vmin, precll,
                       mureg, sigreg, iscov, locnames, niter)
summary.write()