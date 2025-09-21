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
all_samples = list(gen.columns)
train_samples = list(cov_train['sample_id'])
train_samples_pos = [all_samples.index(x) for x in train_samples]

CHUNK = math.floor(len(train_samples_pos)/5)

with open("models_internal_r2_out.tsv","w") as f:
    f.write('splice_id\tn_vars_combo\tmean_int_pr_combo\tmean_int_r2_combo\tn_vars_cis\tmean_int_pr_cis\tmean_int_r2_cis\tn_vars_trans\tmean_int_pr_trans\tmean_int_r2_trans\n')
    f.flush()
    for s in splices:
        
        #
        # Compute internal R2s for combo (cis & trans sQTLs)
        #
        gws_s = gws[gws["phenotype_id"]==s]
        gws_s = gws_s.set_index("variant_rsid")
        resid = list(all_resid[s])
        var_ids = list(gws_s.index)

        # Mean imputation for SNP missingness
        X = np.array(gen.loc[var_ids])
        row_mean = np.nanmean(X, axis=1)
        inds = np.where(np.isnan(X))
        X[inds] = np.take(row_mean, inds[0])
        
        y = np.array(resid)
        X = np.transpose(X)
        n_vars = X.shape[1]

        all_int_pr = []
        all_int_r2 = []
        
        for i in range(1,6):
            test_pos = train_samples_pos[CHUNK*(i-1):CHUNK*i]
            train_pos = [x for x in train_samples_pos if x not in test_pos] # Remove test pos ids from all_ids
        
            x_train, x_test = X[train_pos,:], X[test_pos,:]
            y_train, y_test = y[train_pos], y[test_pos]
        
            model = BayesianRidge(alpha_1=alpha_1, alpha_2=alpha_2, lambda_1=lambda_1, lambda_2=lambda_2)
            model.fit(x_train, y_train)
            y_pred = model.predict(x_test)
            
            all_int_pr.append(pearsonr(y_test,y_pred)[0])
            all_int_r2.append(r2_score(y_test,y_pred))

        # Compute mean of folds
        all_int_pr_mean = stat.mean(all_int_pr)
        all_int_r2_mean = stat.mean(all_int_r2)

        
        #
        # If theres at least one CIS snp, do CIS only, otherwise set avgs to NA
        #   
        gws_s = gws[gws["phenotype_id"]==s]
        gws_s = gws_s.set_index("variant_rsid")
        if "CIS" in list(gws_s["type"].value_counts().keys()):
            gws_s = gws_s[gws_s["type"]=="CIS"]
            resid = list(all_resid[s])
            var_ids = list(gws_s.index)
            X = np.array(gen.loc[var_ids])
            row_mean = np.nanmean(X, axis=1)
            inds = np.where(np.isnan(X))
            X[inds] = np.take(row_mean, inds[0])
            y = np.array(resid)
            X = np.transpose(X)
            n_vars_cis = X.shape[1]
    
            all_int_pr_cis = []
            all_int_r2_cis = []

            for i in range(1,6):
                test_pos = train_samples_pos[CHUNK*(i-1):CHUNK*i]
                train_pos = [x for x in train_samples_pos if x not in test_pos] # Remove test pos ids from all_ids
            
                x_train, x_test = X[train_pos,:], X[test_pos,:]
                y_train, y_test = y[train_pos], y[test_pos]
            
                model = BayesianRidge(alpha_1=alpha_1, alpha_2=alpha_2, lambda_1=lambda_1, lambda_2=lambda_2)
                model.fit(x_train, y_train)
                y_pred = model.predict(x_test)
                
                all_int_pr_cis.append(pearsonr(y_test,y_pred)[0])
                all_int_r2_cis.append(r2_score(y_test,y_pred))

            all_int_pr_cis_mean = stat.mean(all_int_pr_cis)
            all_int_r2_cis_mean = stat.mean(all_int_r2_cis)

        else:
            n_vars_cis = 0
            all_int_pr_cis_mean = "NA"
            all_int_r2_cis_mean = "NA"

        
        #
        # If theres at least one GWS snp, do GWS only, otherwise set avgs to NA
        # 
        gws_s = gws[gws["phenotype_id"]==s]
        gws_s = gws_s.set_index("variant_rsid")
        if "GWS" in list(gws_s["type"].value_counts().keys()):
            gws_s = gws_s[gws_s["type"]=="GWS"]
            resid = list(all_resid[s])
            var_ids = list(gws_s.index)
            X = np.array(gen.loc[var_ids])
            row_mean = np.nanmean(X, axis=1)
            inds = np.where(np.isnan(X))
            X[inds] = np.take(row_mean, inds[0])
            y = np.array(resid)
            X = np.transpose(X)
            n_vars_gws = X.shape[1]
    
            all_int_pr_gws = []
            all_int_r2_gws = []

            for i in range(1,6):
                test_pos = train_samples_pos[CHUNK*(i-1):CHUNK*i]
                train_pos = [x for x in train_samples_pos if x not in test_pos] # Remove test pos ids from all_ids
            
                x_train, x_test = X[train_pos,:], X[test_pos,:]
                y_train, y_test = y[train_pos], y[test_pos]
            
                model = BayesianRidge(alpha_1=alpha_1, alpha_2=alpha_2, lambda_1=lambda_1, lambda_2=lambda_2)
                model.fit(x_train, y_train)
                y_pred = model.predict(x_test)
                
                all_int_pr_gws.append(pearsonr(y_test,y_pred)[0])
                all_int_r2_gws.append(r2_score(y_test,y_pred))

            all_int_pr_gws_mean = stat.mean(all_int_pr_gws)
            all_int_r2_gws_mean = stat.mean(all_int_r2_gws)

        else:
            n_vars_gws = 0
            all_int_pr_gws_mean = "NA"
            all_int_r2_gws_mean = "NA"
        
        write_text = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(s,n_vars,
                                                                                                                   all_int_pr_mean,all_int_r2_mean,                                                                      n_vars_cis,all_int_pr_cis_mean,all_int_r2_cis_mean,
n_vars_gws,all_int_pr_gws_mean,all_int_r2_gws_mean)
        f.write(write_text)
        f.flush()
