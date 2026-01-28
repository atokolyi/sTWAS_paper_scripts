import time
ft = time.time()
import sys
COND=sys.argv[1] # i.e. "disease1"
COHORT=sys.argv[2] # i.e. "biobank1"
print("Running per-cohort B-LORE for (",COND,") in (",COHORT,")",sep="")
print("============================================================")

import pandas as pd
import numpy as np
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

cov = pd.read_csv("demo/inputs/"+COHORT+"/covs.tsv",sep="\t")
cov['id'] = cov['id'].astype(str)
cond = pd.read_csv("demo/inputs/"+COHORT+"/phenome.tsv",sep="\t")
cond['id'] = cond['id'].astype(str)
cond = cond.set_index('id').loc[cov['id']].reset_index()
cov.index = cov['id']
cond.index = cond['id']

valid = cond[COND].notna()
y = cond[COND][valid]

splice = pd.read_csv("demo/inputs/"+COHORT+"/expression_mat.tsv",sep="\t",index_col=0)
features = list(splice.columns)
splice.drop(splice[~valid].index,inplace=True)
splice = splice.to_numpy(dtype=np.float64)

covs = cov.to_numpy()[:,1:]
covs = np.asarray(covs,dtype=np.float64)

feature_meta = pd.read_csv("demo/inputs/"+COHORT+"/feature_meta.tsv",header=0,sep="\t")
feature_meta = feature_meta.set_index('feature', drop=False)
feature_meta = feature_meta.loc[features]
features = list(features)

block_keep = pd.read_csv("demo/inputs/"+COHORT+"/block_keep.tsv",sep="\t",header=0)

# Limit LD blocks with sTWAS significant features
sig_splices = list(set(list(block_keep[block_keep['cond']==COND]['splice'])))
keep_blocks = list(set(list(feature_meta[feature_meta['feature'].isin(sig_splices)]['ld_block'])))
feature_meta = feature_meta[feature_meta['ld_block'].isin(keep_blocks)]
keep_i = [features.index(x) for x in list(feature_meta['feature'])]
splice = splice[:,keep_i]
features = list(np.array(features)[keep_i])

gt = [splice[:,[features.index(x) for x in feature_meta[feature_meta['ld_block']==y]['feature']]] for y in list(set(list(feature_meta['ld_block'])))]

phenotype = y==1
phenotype = tuple(np.array(phenotype, dtype=np.float64))
gt = tuple([x for x in gt])

# Using default optimization parameters
mureg = 0.0
sigreg = 0.01
tolerance = 0.01

print ("Optimizing regularizer...")
t = time.time()
sigreg_optim = OptimizeRegularizer(gt, phenotype, mureg, sigreg, tolerance, covs)
sigreg_optim.update()
mureg = sigreg_optim.mureg
sigreg = sigreg_optim.sigmareg
niter = sigreg_optim.niter
print("\t - Done in ",round((time.time()-t),4)," secs",sep="")

covinfo = [["age","sex"]]

# Run the Bayesian multiple logistic regression
t = time.time()
print ("Calculating summary statistics...")
logreg = LogisticRegression(gt, phenotype, mureg, sigreg, covariates = covs)
logreg.fit()
v0 = logreg.v0
vmin = logreg.vmin
precll = logreg.precll
iscov = logreg.iscov
print("\t - Done in ",round((time.time()-t),4)," secs",sep="")

outdir = "demo/outputs/"+COHORT+"/"+COND
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
print("Saving results to: ",outdir)
dupes = [[] for x in range(len(loci))]
freq  = [[1 for x in range(len(loci[y]))] for y in range(len(loci))]
summary = WriteSummary(outdir, file_prefix)
summary.set_statistics(snpinfo, covinfo, dupes, freq, v0, vmin, precll,
                       mureg, sigreg, iscov, locnames, niter)
summary.write()
print("\t - Done!")

