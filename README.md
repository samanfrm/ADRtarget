# Adverse drug reaction (ADR) - Target association prediction
Machine learning guided association of adverse drug reactions with in vitro off-target pharmacology 

## Overview
Adverse drug reactions (ADRs) are one of the leading causes of morbidity and mortality in health care. Understanding which drug targets are linked to ADRs can lead to the development of safer medicines. 

Here, we analyze in vitro secondary pharmacology of common (off) targets for 2134 marketed drugs. To associate these drugs with human ADRs, we utilized FDA Adverse Event Reports and developed Random Forest models that predict ADR occurrences from in vitro pharmacological profiles. 

By evaluating Gini importance scores of model features, we identify 250 target-ADR associations. Among these are established relations, such as the association of in vitro hERG binding with cardiac arrhythmias, which validate our machine learning approach. Evidence on bile acid metabolism supports our identification of associations between the Bile Salt Export Pump and renal, thyroid, lipid metabolism, skin, respiratory tract and central nervous system disorders. 

Unexpectedly, our model suggests PDE3 is associated with 44 ADRs, including congenital renal disorders. These associations provide a comprehensive resource to support drug development and human biology studies.

## System requirements: 

### All software dependencies and operating systems (including version numbers):
OS: MacOS 10.13 or higher  (Linux OS and Windows 10 work as well after necessary configurations)
Python version 3.7
R version 3.5.0
RStudio version 1.1.453

### Versions the software has been tested on:
Python v3.6 and 3.7
R versions 3.4.1 and 3.5.0
RStudio version 1.1.453

### Any required non-standard hardware:
computer with sufficient RAM: 8Gb or more

## Installation Guide

### Instructions:
-download/install python 3.7 from python.org
-in shell:
pip install virtualenv
cd ~/
virtualenv -p /usr/local/bin/python3 py3 #creates virtual environment named py3 using installed python3.7 version
source py3/bin/activate #activate virtual environment
pip install jupyter
pip install pandas
pip install --upgrade numpy
pip install --upgrade scipy
pip install requests
pip install statsmodels
pip install jellyfish

## Reproducibility analysis between 5 replicate runs
Identifying for which ADRs the model only predicts zeros: exclude those ADRs (below)

```python
version_model='v1'
level='HLGT'#'SOC'#
SEEDS=[112358,2662851,49,5332728]
ADR_excl=dict()
for i in range(5):
    filename='summ_metrics_ADRs.csv'
    if i == 0:
        path='/Users/horizon/Documents/HMS/Novartis2018Hackathon/PCS/Seda_modelling/3classes/'+ \
            version_model+'_'+level+'/'
        MetADR = pd.read_csv(path+filename,index_col=0)
    else:
        seed=SEEDS[i-1]
        path='/Users/horizon/Documents/HMS/Novartis2018Hackathon/PCS/Seda_modelling/v1_HLGT_diffSeed_'+str(seed)+'/'
        MetADR_OLD=copy.deepcopy(MetADR['ADR'])
        MetADR = pd.read_csv(path+filename,index_col=0)
        MetADR['ADR']=MetADR_OLD
#     MetADR
#     len(MetADR)
    MetADR_nonan=MetADR[~np.isnan(MetADR['precision'])]
    MetADR_nonan['Pos']=MetADR_nonan['TP']+MetADR_nonan['FN']

    ADR_excl[i]=MetADR[np.isnan(MetADR['precision'])]['ADR']
    ADR_excl[i]=list(ADR_excl[i])
    for j in range(len(ADR_excl[i])):
        ADR_excl[i][j]=re.sub('\(|\)|,|\'|-',' ',ADR_excl[i][j])
    print(len(ADR_excl[i]))
```

First use meddra sheet to get mapping from HLGT to SOC level for ordering in heatmap

```python
path='/Users/horizon/Documents/HMS/Novartis2018Hackathon/PCS/'
filename='meddr1a_full_12092006.xlsx'
meddra = pd.read_excel(path+filename)

HLGT2SOC=dict()
for index, arow in meddra.iterrows():
    hlgt = str(arow['HLGT_TXT']).lower()
    hlgt=re.sub('\(|\)|,|\'|-',' ',hlgt)
    try:
        HLGT2SOC[hlgt]
    except KeyError:
        HLGT2SOC[hlgt] = str(arow['SOC_TXT']).lower()
#HLGT2SOC
```

### Load all unfiltered target-ADRs (as above for example seed 49)

```python
version_model='v1'
level='HLGT'
path='/Users/horizon/Documents/HMS/Novartis2018Hackathon/PCS/'
SEEDS=[112358,2662851,49,5332728]

AT=dict()
filename=version_model+'_'+level+'_ADR_target_assoc_annotated.csv'#original seed (7332)
AT[0]= pd.read_csv(path+filename,index_col=0)
AT[0]=AT[0].rename(index=str, columns={"pos correlation (between target IC50 and ADR prob)": "pos_corr"})
for i in range(1,5):
    seed=SEEDS[i-1]
    filename=version_model+'_'+level+'_diffSeed'+str(seed)+'_ADR_target_assoc.csv'
    AT[i]= pd.read_csv(path+filename,index_col=0)

for i in range(5):
    SOC=[]
    for j in AT[i].index:
        SOC.append(HLGT2SOC[AT[i]['ADR'][j]])
    AT[i].insert(loc=1,column='SOC',value=SOC)
```

### Filtering out ADRs from AT
Filter out ADRs that have no predictive value (ADR_excl)

Also filter out the following SOC classes: 

-general disorders and administration site conditions

-injury, poisoning and procedural complications

-poisoning and procedural complications

-surgical and medical procedures

-neoplasms benign, malignant and unspecified (incl cysts and polyps)

-investigations

-social circumstances

```python
SOC_excl=['surgical and medical procedures',
        'general disorders and administration site conditions',
        'injury, poisoning and procedural complications',
        'investigations',
        'poisoning and procedural complications',
        'neoplasms benign, malignant and unspecified (incl cysts and polyps)',
        'social circumstances']
for i in range(5):
    AT[i]=AT[i][~AT[i]['ADR'].isin(ADR_excl[i])]
    AT[i]=AT[i][~AT[i]['SOC'].isin(SOC_excl)]
    AT[i]=AT[i].sort_values(by=['SOC','ADR','target'])
    
AT[0]=AT[0].drop(['pos_corr','no_corr'],axis=1)
AT[0]=AT[0].reindex(columns=['SOC','ADR','target','assay','0-3uM','3-30uM','>30uM'])
for i in range(1,5):
    AT[i]=AT[i].drop(['no_corr'],axis=1)
    AT[i]=AT[i].reindex(columns=['SOC','ADR','target','assay','0-3uM','3-30uM','>30uM'])
```
### Merge (union = outer) to see overlaps

```python
ATM = pd.merge(AT[0],AT[1], on=['SOC','ADR','target','assay'],how='outer')
ATM.columns = ['SOC','ADR','target','assay','0:0-3uM','0:3-30uM','0:>30uM','1:0-3uM','1:3-30uM','1:>30uM']
for i in range(2,5):
    ATM= pd.merge(ATM,AT[i], on=['SOC','ADR','target','assay'],how='outer')
    colnames={'0-3uM':str(i)+':0-3uM','3-30uM':str(i)+':3-30uM','>30uM':str(i)+':>30uM'}
    ATM=ATM.rename(columns=colnames)#index=str, 
    
path='/Users/horizon/Documents/HMS/Novartis2018Hackathon/PCS/'
ATM=ATM.sort_values(by=['SOC','ADR','target'])
print(len(ATM))
```

## Get summary statistics

```python
ct={1:[],2:[],3:[]}
cl_dict={1:'0-3uM',2:'3-30uM',3:'>30uM'}
for j in ATM.index:
    for cl in cl_dict.keys():
        count=0
        for i in range(5):
            if not np.isnan(ATM[str(i)+':'+cl_dict[cl]][j]):
                count+=1
        ct[cl].append(count)
ATM.insert(loc=4,column='N:'+cl_dict[1],value=ct[1])
ATM.insert(loc=5,column='N:'+cl_dict[2],value=ct[2])
ATM.insert(loc=6,column='N:'+cl_dict[3],value=ct[3])

#get averages over reps
ct={1:[],2:[],3:[]}
se={1:[],2:[],3:[]}
for j in ATM.index:
    for cl in cl_dict.keys():
        ct[cl].append( \
    np.mean(ATM.loc[j,['0:'+cl_dict[cl],'1:'+cl_dict[cl],'2:'+cl_dict[cl],'3:'+cl_dict[cl],'4:'+cl_dict[cl]]]))
        se[cl].append( \
    np.std(ATM.loc[j,['0:'+cl_dict[cl],'1:'+cl_dict[cl],'2:'+cl_dict[cl],'3:'+cl_dict[cl],'4:'+cl_dict[cl]]]) \
                      / np.sqrt(ATM.loc[j,'N:'+cl_dict[cl]]))    
ATM.insert(loc=7,column='mean:'+cl_dict[1],value=ct[1])
ATM.insert(loc=8,column='mean:'+cl_dict[2],value=ct[2])
ATM.insert(loc=9,column='mean:'+cl_dict[3],value=ct[3])
ATM.insert(loc=10,column='sem:'+cl_dict[1],value=se[1])
ATM.insert(loc=11,column='sem:'+cl_dict[2],value=se[2])
ATM.insert(loc=12,column='sem:'+cl_dict[3],value=se[3])

```

### Filter out non-significant associations

```python
from scipy.stats import ttest_ind, ttest_rel
from statsmodels.stats.multitest import fdrcorrection
alpha_FDR=0.1
pvals=[]
cl_dict={1:'0-3uM',2:'3-30uM',3:'>30uM'}

#filter out ADR-targets that were only found 1: not reproducible
temp=ATM[['N:'+cl_dict[1],'N:'+cl_dict[2],'N:'+cl_dict[3]]]>1
ATM=ATM[temp.any(axis=1)]

for j in ATM.index:
    count=0
    data_sample1=[]
    data_sample2=[]
    for cl in cl_dict.keys():
        if ATM.loc[j,'N:'+cl_dict[cl]]>1: 
            if count==0:
                data_sample1=list(ATM.loc[j,['0:'+cl_dict[cl],'1:'+cl_dict[cl],'2:'+cl_dict[cl],'3:'+cl_dict[cl],'4:'+cl_dict[cl]]])
                count+=1
            elif count==1:
                data_sample2=list(ATM.loc[j,['0:'+cl_dict[cl],'1:'+cl_dict[cl],'2:'+cl_dict[cl],'3:'+cl_dict[cl],'4:'+cl_dict[cl]]])
    stat,pval=ttest_ind(data_sample1,data_sample2,nan_policy='omit')#independent samples t test

    if np.isnan(stat):#nan data so drop and do not test for significance
        print(stat,pval,'drop',j)
        ATM=ATM.drop(j) 
    else:
        pvals.append(pval)


BOOL, qvals = fdrcorrection(pvals, alpha=1)
ATM.insert(loc=7,column='pval',value=pvals)
ATM.insert(loc=8,column='padj',value=qvals)
ATM=ATM[ATM['padj']<alpha_FDR]



ATM=ATM.sort_values(by=['SOC','ADR','target'])
path='/Users/horizon/Documents/HMS/Novartis2018Hackathon/PCS/'
filename=version_model+'_'+level+'_ADR_target_assoc_final_demo.csv'
ATM.to_csv(path+filename)
print(len(ATM))
```
