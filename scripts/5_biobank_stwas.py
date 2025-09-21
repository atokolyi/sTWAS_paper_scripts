import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
import statsmodels.api as sm
import subprocess
import sys

# Load subset UKBB COVs and PheCodes (https://github.com/umich-cphds/createUKBphenome)
cov = pd.read_csv("ukbb/covs_wbs_sexmatch.tsv",sep="\t")
cond = pd.read_pickle("phenome.pickle");
cond = cond.set_index('IID').loc[cov['id']].reset_index()

# Parallelise per splicing event
SPLICE_N = int(sys.argv[1])+1
STOP_N = SPLICE_N+1
with open(f"/tmp/ukbb.twas.{SPLICE_N}", "w") as out:
    subprocess.run(
        ["sed", "-n", f"1p;{SPLICE_N}p;{STOP_N}q;", "ukbb_splice_resid.tsv"],
        stdout=out,
        check=True
    )
splice = pd.read_csv(f"/tmp/ukbb.twas.{SPLICE_N}", sep="\t")
subprocess.run(["rm", f"/tmp/ukbb.twas.{SPLICE_N}"], check=True)

# Align covariate and predicted splicing matricies
splice_name = splice['Splice'][0]
splice.drop('Splice',axis=1,inplace=True)
splice = splice.T
splice.columns = ['splice']
splice.index = [int(x) for x in splice.index]
splice = splice.loc[cov['id']].reset_index()

cov['resid'] = splice['splice']

conds = list(cond.columns[1:])
MAX = len(conds)

df = []

# Run sTWAS for all PheCodes for this splicing event
for c in range(0,MAX):
    cov['cond'] = cond.iloc[:, 1+c]
    cond_name = conds[c]
    formula = 'cond ~ resid+sex+age'
    model_glm = smf.glm(formula=formula, data=cov, family=sm.families.Binomial()).fit()
    res = model_glm.summary2().tables[1]
    out = {"splice": splice_name,
       "cond":cond_name,
       "glm_est":res['Coef.']['resid'],
       "glm_se":res['Std.Err.']['resid'],
      "glm_z":res['z']['resid'],
      "glm_p":res['P>|z|']['resid'],
      "glm_low":res['[0.025']['resid'],
      "glm_high":res['0.975]']['resid']}
    df.append(out)

df_final = pd.DataFrame.from_dict(df)
df_final.to_csv("stwas_out/"+str(SPLICE_N-1)+".tsv",sep="\t",index=False)
