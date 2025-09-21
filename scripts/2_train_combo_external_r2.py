import math
import pickle

import numpy as np
import pandas as pd
import statistics as stat

from sklearn.metrics import r2_score
from scipy.stats import pearsonr
from sklearn.linear_model import BayesianRidge
from pandas_plink import read_plink

# Read the pruned cis and trans sQTLs ("combo") for model inclusion (from 80% train subset sQTLs)
gws = pd.read_table("combo_full_subset_pruned.tsv")

# Load plink bed of these sQTL SNPs
bim,fam,bed = read_plink("combo_pruned")
gen = bed.compute()
gen = pd.DataFrame(data=gen,index=bim["snp"])

# Read the residualized matrix of splicing events
with open("splice_resid_mat.pickle","rb") as f:
    all_resid = pickle.load(f)

splices = list(set(list(gws["phenotype_id"])))
len(splices)

# Using uninformative priors for BRR (https://doi.org/10.1038/s41586-023-05844-9)
alpha_1 = 10e-5
alpha_2 = 10e-5
lambda_1= 10e-5
lambda_2 = 10e-5

all_ids = list(range(gen.shape[1]))
cov_train = pd.read_table("sqtl_cov_train80.tsv") # Randomized 80% train cohort subset
cov_test = pd.read_table("sqtl_cov_test20.tsv") # Randomized 20% test cohort subset
all_samples = list(gen.columns)
train_samples = list(cov_train['sample_id'])
test_samples = list(cov_test['sample_id'])
train_samples_pos = [all_samples.index(x) for x in train_samples]
test_samples_pos = [all_samples.index(x) for x in test_samples]

CHUNK = math.floor(len(train_samples_pos)/5)

with open("combo_external_r2_out.tsv","w") as f:
    f.write('splice\tn_vars\tmean_int_pr\tmean_int_r2\tmean_ext_pr\tmean_ext_prp\tmean_ext_r2\tmean_int_pr_1snp\tmean_int_r2_1snp\tmean_ext_pr_1snp\tmean_ext_r2_1snp\tn_vars_cis\tmean_int_pr_cis\tmean_int_r2_cis\tmean_ext_pr_cis\tmean_ext_r2_cis\tn_vars_gws\tmean_int_pr_gws\tmean_int_r2_gws\tmean_ext_pr_gws\tmean_ext_r2_gws\n')
    f.flush()
    for s in splices:
        gws_s = gws[gws["phenotype_id"]==s]
        gws_s = gws_s.set_index("variant_rsid")
        resid = list(all_resid[s])
        var_ids = list(gws_s.index)
        X = np.array(gen.loc[var_ids])
        row_mean = np.nanmean(X, axis=1)
        inds = np.where(np.isnan(X))
        X[inds] = np.take(row_mean, inds[0])
        y = np.array(resid)
        X = np.transpose(X)
        n_vars = X.shape[1]

        all_ext_pr = []
        all_ext_r2 = []
        
        for i in range(1,6):
            test_pos_internal = train_samples_pos[CHUNK*(i-1):CHUNK*i]
            train_pos = [x for x in train_samples_pos if x not in test_pos_internal] # Remove test pos ids from all_ids
        
            x_train = X[train_pos,:]
            y_train = y[train_pos]
        
            model = BayesianRidge(alpha_1=alpha_1, alpha_2=alpha_2, lambda_1=lambda_1, lambda_2=lambda_2)
            model.fit(x_train, y_train)
           
            x_ext = X[test_samples_pos,:]
            y_ext = y[test_samples_pos]
            
            y_pred_ext = model.predict(x_ext)
            all_ext_pr.append(pearsonr(y_ext,y_pred_ext)[0])
            all_ext_r2.append(r2_score(y_ext,y_pred_ext))

        all_ext_pr_mean = stat.mean(all_ext_pr)
        all_ext_r2_mean = stat.mean(all_ext_r2)
        
        write_text = "{}\t{}\t{}\t{}\n".format(s,n_vars,
    all_ext_pr_mean,all_ext_r2_mean)
        f.write(write_text)
        f.flush()

