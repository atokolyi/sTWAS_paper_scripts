import pandas as pd
import numpy as np
import sys
import time
import pickle
import gc
import sys
import collections
import re
import os

print("Running B-LORE meta-analysis",sep="")
ZMAX=int(sys.argv[1]) # i.e. 10
PROB_THRESH=float(sys.argv[2]) # i.e. 0.5

# Load B-LORE source functions for adaptation (https://github.com/soedinglab/b-lore)
sys.path.append('b-lore/blore')
from inference.logistic_regression import LogisticRegression
from inference.optimize_regularizer import OptimizeRegularizer
from meta import optimize_hyperparameters
from iotools.io_summary import WriteSummary

# Run meta-analysis for PheCodes which converged in both cohorts

# Find cohorts
COHORTS = [f for f in os.listdir("demo/outputs/") if not (f.startswith('.') or "_meta_" in f)]
print("Found ",len(COHORTS)," cohort(s): ",", ".join(COHORTS))

# Find diseases that converged in all cohorts
done = set([f for f in os.listdir("demo/outputs/"+COHORTS[0]+"/") if not f.startswith('.')])
for i in range((len(COHORTS)-1)):
    biobank_n_done = set([f for f in os.listdir("demo/outputs/"+COHORTS[i+1]+"/") if not f.startswith('.')])
    done = set(done & biobank_n_done)
done = sorted(list(done))
print("Found ",len(done)," disease(s): ",", ".join(done))

for COND in done:
    print("Running meta-analysis for (",COND,") with ZMAX=",ZMAX," and PROB_THRESH=",PROB_THRESH,sep="",end="\n\n")
    statdir = ["demo/outputs/biobank1/"+COND,"demo/outputs/biobank2/"+COND]
    locusnamesfile = "demo/outputs/biobank1/"+COND+"/" + COND + ".locusnames"
    featurefiles=None
    zmax=ZMAX # Up to 10 causal variants per LD block-PheCode pair
    muvar=False
    mparams=None
    outdir = "demo/outputs/"+COND+"_meta_z"+str(ZMAX)
    file_prefix = COND
    optimize_hyperparameters(statdir, locusnamesfile, featurefiles, 
                         outdir, file_prefix, zmax, muvar, mparams)

    DIR="demo/outputs/"+COND+"_meta_z"+str(ZMAX)+"/"+COND+"_res"
    BLOCKS=[f for f in os.listdir(DIR) if not f.startswith('.')]
    BLOCKS=[b.replace(".res","") for b in BLOCKS if ".res" in b]
    
    all_res = []
    for block in BLOCKS:
        res = pd.read_csv("demo/outputs/"+COND+"_meta_z"+str(ZMAX)+"/"+COND+"_res/"+block+".res",
                            skiprows=1, header=None, sep='\\s+')
        res = res.iloc[:, [0,4,5,6,7,8]]
        res.columns = ["feature","prob","vmin","pi","mu","sigma"]
        res['block'] = block
        all_res.append(res)
    all_res = pd.concat(all_res)
    all_res = all_res.sort_values(by="prob",ascending=False)
    
    all_res_pp = all_res[all_res['prob']>PROB_THRESH]

    print("\nFinished B-LORE meta-analysis for (",COND,")",sep="")
    if len(all_res_pp)>0:
        print("\t - ",len(all_res_pp),"/",len(all_res)," features with posterior probability > ",PROB_THRESH*100,"%:",sep="")
        for index, row in all_res_pp.iterrows():
            print("\t\t - ",row['feature'],":\tPP=",round(row['prob']*100,2),"%",sep="")
        
    else:
        print("\t - 0 features with posterior probability > ",PROB_THRESH*100,"%:",sep="")
    
    all_res.to_csv("demo/outputs/"+COND+"_meta_z"+str(ZMAX)+"_results.tsv",
                      index=False,header=True,sep="\t")
print("Done!")