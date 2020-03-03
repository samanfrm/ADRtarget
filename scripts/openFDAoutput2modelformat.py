#!/usr/bin/env python
# coding: utf-8
# Robert Ietswaart 
# GNU General Public License v3.0

### Citation:
#Robert Ietswaart<sup>\*,#</sup>, Seda Arat<sup>\*,#</sup>, Amanda X. Chen<sup>\*</sup>, 
#Saman Farahmand<sup>\*</sup>, Bumjun Kim, William DuMouchel, 
#Duncan Armstrong, Alexander Fekete, Jeffrey J. Sutherland<sup>#</sup>, Laszlo Urban<sup>#</sup>  
#*Machine learning guided association of adverse drug reactions with in vitro target-based 
#pharmacology* (2019), [BioRxiv; 750950](https://www.biorxiv.org/content/10.1101/750950v2).

import sys
import os
import pandas as pd
import numpy as np
import copy
from scipy.stats import binom
from statsmodels.stats.multitest import fdrcorrection

# ### Load the openFDA output
repo_path = str(sys.argv[1])
path = os.path.join(repo_path,'data')
filename = 'v1_compounds_FDA.csv'
pulled = pd.read_table(os.path.join(path,filename),sep='\t',index_col=0)
pulled = pulled.rename(columns={'#reports': 'N_reports','HGLT_per_report':'HLGT_per_report'})
pulled = pulled.drop('Unnamed: 0.1',axis=1)
pulled.loc[pulled['N_reports']=='error','N_reports']=-1
pulled['N_reports'] = pulled['N_reports'].astype(dtype='int32')


# ### Merge of multiple compound names
# In pulled, some (unique) compound ID numbers occur multiple times as there are multiple names (1 row per generic name). Additionally when OpenFDA is queried for a particular generic name, it might return the same ADR report as for another name of the same compound. Here you merge the result for each generic name into one row for each unique compound ID.

UNIQ_number=list(set(pulled['number']))
for i in range(len(UNIQ_number)):
    TEMP=list(copy.deepcopy(pulled['number'][pulled['number']==UNIQ_number[i]].index))
#     print(TEMP)
   
    if len(TEMP)>1:
        if max(pulled.loc[TEMP,'N_reports'])>0:#there are reports
#             print(max(pulled.loc[TEMP,'N_reports']))
            MAX_ID=pulled['N_reports'][TEMP].idxmax()#get ID of row with max number of reports, merge rest into that row
#             print(MAX_ID)
            mrep_id=pulled['reports_id'][MAX_ID].split('|||')
            
            for idx in TEMP:
                if idx != MAX_ID:
                    if pulled.loc[idx,'N_reports'] > 0:
                        reps=pulled['reports_id'][idx].split('|||')
                        pts=pulled['PTs_per_report'][idx].split('|||')
                        hlgts=pulled['HLGT_per_report'][idx].split('|||')
                        socs=pulled['SOC_per_report'][idx].split('|||')

                        for j in range(len(reps)):
                            if reps[j] not in mrep_id:
                                mrep_id.append(reps[j])
                                pulled.loc[MAX_ID,'PTs_per_report']=pulled.loc[MAX_ID,'PTs_per_report']+'|||'+pts[j]
                                pulled.loc[MAX_ID,'HLGT_per_report']=pulled.loc[MAX_ID,'HLGT_per_report']+'|||'+hlgts[j]
                                pulled.loc[MAX_ID,'SOC_per_report']=pulled.loc[MAX_ID,'SOC_per_report']+'|||'+socs[j]
                    
                    pulled=pulled.drop(idx)#delete the row
#                     print('del idx',idx)
            pulled.loc[MAX_ID,'N_reports']=len(mrep_id)
            pulled.loc[MAX_ID,'reports_id']='|||'.join(mrep_id)
            
        else:#none of them have any reports, keep only first row, drop rest
            pulled=pulled.drop(TEMP[1:])
#             print('del idx',TEMP[1:])

#drop the column name as we have merged the results for different generic names
pulled = pulled.drop(['name'], axis=1)


# ### Calculate occurrence proportions and their significance (1 or 0)
def P_binomial_ocr(ocr_counts, N_reps,N_t):
    p1=1.0/N_t
    pval = 1-binom.cdf(ocr_counts-1, N_reps, p1)#pval=prob(x >= ocr)=1-prob(x<= ocr-1)=1-cdf(ocr-1)
    return pval

def discretize_ocr(dict_pval):
    alpha_FDR=0.01
    pvals=list(dict_pval.values())
    BOOL, qvals = fdrcorrection(pvals, alpha=alpha_FDR)
    keys = list(dict_pval.keys())
    for i in range(len(keys)):
        if BOOL[i]: 
            dict_pval[keys[i]]= 1
        else:
            dict_pval[keys[i]]= 0
    return dict_pval

def get_ocr(onerow,level,level_list):
    if onerow[level] is not np.nan:
        N_terms=len(level_list)
        N_r=onerow['N_reports']
        ocr_raw = dict.fromkeys(level_list, 0)
        ocr_discretized = dict.fromkeys(level_list, 0)
        reps=onerow[level].split('|||')
        for r in reps:
            terms=set(r.split(';'))
            for t in terms:
                try:
                    ocr_raw[t]+=1
                except KeyError:
                    pass
        for akey in ocr_raw.keys():
            ocr_discretized[akey] = P_binomial_ocr(ocr_raw[akey],N_r,N_terms)#generate pvals
            ocr_raw[akey] = float(ocr_raw[akey])/N_r
        ocr_discretized=discretize_ocr(ocr_discretized)
    else:
        ocr_raw = {}
        ocr_discretized = {}
    return ocr_raw, ocr_discretized

def calc_occurences(pulled,term_list,level):
    res_ocr_raw = {}
    res_ocr_disc = {}
    for asoc in term_list:
        res_ocr_raw[asoc] = []
        res_ocr_disc[asoc] = []

    for index, arow in pulled.iloc[:].iterrows():
        raw_tmp, disc_tmp  = get_ocr(arow,level,term_list)
        for asoc in res_ocr_raw.keys():
            try:
                res_ocr_raw[asoc].append(raw_tmp[asoc])
            except KeyError:
                res_ocr_raw[asoc].append(np.NaN)
            try:
                res_ocr_disc[asoc].append(disc_tmp[asoc])
            except KeyError:
                res_ocr_disc[asoc].append(np.NaN)
    
    return res_ocr_raw , res_ocr_disc

# ### function to output results to csv
def ocr_to_csv(res_dict,pulled,path,filename):
    res_df = pd.DataFrame(res_dict)
    res_df.insert(loc=0,column='number', value=pulled['number'].tolist())
    res_df.insert(loc=1,column='N_reports', value=pulled['N_reports'].tolist())
    res_df.to_csv(os.path.join(path,filename))
    return res_df

# ### Write SOC level results to csv
#SOC = System Organ Classes
soc_acc = []
for i in pulled['SOC_per_report'].tolist():
    if i is not np.nan:
        i=i.replace('|||',';')
        soc_acc.extend(i.split(';'))
soc_list = set(soc_acc)
soc_list.remove('NaN')

res_ocr_raw, res_ocr_discretized = calc_occurences(pulled,soc_list,'SOC_per_report')

filename = 'v1_compounds_FDA_model_format_SOC_ocr_prob.csv'
res_df = ocr_to_csv(res_ocr_raw, pulled, path, filename)

filename = 'v1_compounds_FDA_model_format_SOC_ocr_bool.csv'
res_df = ocr_to_csv(res_ocr_discretized,pulled,path,filename)


# ### Write HLGT level results to csv
# HLGT = Higher level Group Terms
hlgt_acc = []
for i in pulled['HLGT_per_report'].tolist():
    if i is not np.nan:
        i=i.replace('|||',';')
        hlgt_acc.extend(i.split(';'))
hlgt_list = set(hlgt_acc)
hlgt_list.remove('NaN')

res_ocr_raw,res_ocr_discretized = calc_occurences(pulled,hlgt_list,'HLGT_per_report')

filename = 'v1_compounds_FDA_model_format_HLGT_ocr_prob.csv'
res_df = ocr_to_csv(res_ocr_raw,pulled,path,filename)

filename = 'v1_compounds_FDA_model_format_HLGT_ocr_bool.csv'
res_df = ocr_to_csv(res_ocr_discretized,pulled,path,filename)
