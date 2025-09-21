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

# Run meta-analysis for PheCodes which converged in both cohorts
aou_done = set(os.listdir("aou_out"))
ukbb_done = set(os.listdir("ukbb_out"))
done = list(aou_done & ukbb_done)
done = sorted(done)
for COND in done:
    statdir = ["aou_out/"+COND,"ukbb_out/"+COND]
    locusnamesfile = "aou_out/"+COND+"/" + COND + "_out.locusnames"
    featurefiles=None
    zmax=10 # Up to 10 causal variants per LD block-PheCode pair
    muvar=False
    mparams=None
    outdir = "meta_out_z10/"+COND
    file_prefix = COND
    optimize_hyperparameters(statdir, locusnamesfile, featurefiles, 
                         outdir, file_prefix, zmax, muvar, mparams)
