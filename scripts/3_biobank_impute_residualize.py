import math
import pickle

import numpy as np
import pandas as pd
import statistics as stat
import statsmodels.formula.api as smf

from sklearn.metrics import r2_score
from scipy.stats import pearsonr
from sklearn.linear_model import BayesianRidge
from pandas_plink import read_plink

# Read the pruned cis and trans sQTLs ("combo") for model inclusion (from full cohort)
## For splicing events with witheld test subset R2>0.01
gws = pd.read_table("combined_full_subset_ukbb_pruned.tsv")

# Load full INTERVAL genotypes
bim,fam,bed = read_plink("combo_pruned")
gen = bed.compute()
gen = pd.DataFrame(data=gen,index=bim["snp"])

splices = list(set(list(gws["phenotype_id"])))

len(splices)

alpha_1 = 10e-5
alpha_2 = 10e-5
lambda_1= 10e-5
lambda_2 = 10e-5

all_ids = list(range(gen.shape[1]))
all_samples = list(gen.columns)

# Load UKBB genotypes [subset to combo sQTL SNPs]
bim,fam,bed = read_plink("ukbb_combo_pruned")
gen_ukbb = bed.compute()
gen_ukbb = pd.DataFrame(data=gen_ukbb,index=bim["snp"],columns=fam['iid'])

# Load UKBB sample PCs
pcs = pd.read_table("ukbb_PCs.tsv",index_col=0)
pcs.index = [str(x) for x in list(pcs.index)]


with open("ukbb_splice_resid.tsv","w") as fr:
    fr.write("Splice\t"+"\t".join(list(pcs.index))+"\n")
    fr.flush()

    # Train on full INTERVAL cohort, and predict on UKBB
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

        model = BayesianRidge(alpha_1=alpha_1, alpha_2=alpha_2, lambda_1=lambda_1, lambda_2=lambda_2)
        model.fit(X, y)

        X = np.array(gen_ukbb.loc[var_ids])
        row_mean = np.nanmean(X, axis=1)
        inds = np.where(np.isnan(X))
        X[inds] = np.take(row_mean, inds[0])
        X = np.transpose(X)
        y_pred_ukbb = model.predict(X)

        # Regress out 17 genotype PCs from prediction
        df = pd.DataFrame({
            "id": gen_ukbb.columns,
            "pred": y_pred_ukbb
        })
        df.index = df['id']
        df = df.loc[list(pcs.index)] # Limit the output to those with PCs/QC subset
        pcs['pred'] = df['pred']
        form = "pred ~ PC" + " + PC".join([str(x) for x in list(range(1,17))])
        lm = smf.ols(form, data=pcs, missing='drop').fit()
        df['resid'] = lm.resid
        df = df.drop(['id'],axis=1)
        fr.write(s+"\t"+"\t".join(["{:.4f}".format(x) for x in list(df['resid'])])+"\n")
        fr.flush()