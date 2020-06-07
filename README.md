# Adverse drug reaction (ADR) - Target association prediction
Machine learning guided association of adverse drug reactions with in vitro off-target pharmacology 

## Overview
Adverse drug reactions (ADRs) are one of the leading causes of morbidity and mortality 
in health care. Understanding which drug targets are linked to ADRs can lead to the 
development of safer medicines. 

Here, we analyze in vitro secondary pharmacology of common (off) targets for 2134 marketed 
drugs. To associate these drugs with human ADRs, we utilized FDA Adverse Event Reports and 
developed random forest models that predict ADR occurrences from in vitro pharmacological 
profiles. 

By evaluating Gini importance scores of model features, we identify 221 target-ADR 
associations, which co-occur in PubMed abstracts to a greater extent than expected 
by chance. Among these are established relations, such as the association of in vitro 
hERG binding with cardiac arrhythmias, which further validate our machine learning approach. 

Evidence on bile acid metabolism supports our identification of associations between the 
Bile Salt Export Pump and renal, thyroid, lipid metabolism, respiratory tract and central 
nervous system disorders. Unexpectedly, our model suggests PDE3 is associated with 40 ADRs. 
These associations provide a comprehensive resource to support drug development and human 
biology studies.

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

- Download/install python 3.7 from python.org
- In shell:
```shell
pip install virtualenv
cd ~/
virtualenv -p /usr/local/bin/python3 py3 #creates virtual environment py3 using python v3.7
source py3/bin/activate #activate virtual environment
pip install jupyter
pip install pandas
pip install --upgrade numpy
pip install --upgrade scipy
pip install requests
pip install statsmodels
pip install jellyfish
```
R/Rstudio:
The model is generated using utiml version 0.1.4
```R
Install.packages (c (“randomForest”, “mldr”, “utiml”, “mltools”), dependencies = T) 
library (randomForest) # Load libraries for these packages
library (mldr)
library (utiml)
options (utiml.empty.prediction = T)
library (mltools)
```

Install the ADRtarget repository
```
git clone https://github.com/samanfrm/ADRtarget
```

## Typical install time: 
less than 1h


### OpenFDA Adverse Event Reports query:
Argument 1: start number of queried compounds from `data/compounds.csv`.  
Argument 2: end number of queried compounds from `data/compounds.csv`.  
Argument 3: open FDA api key (to ensure fast run time speed): can be requested 
at https://open.fda.gov/apis/authentication/ .    
Argument 4: path to this local git ADRtarget repository.  
Run in shell:
```shell
python openFDAapi.py 0 2134 [openFDA api key] [local path]
```
0 and 2134 indicate the begin and end number of compounds (see data/compounds.csv)
that are used for querying FDA. The script can be run in smaller chunks of number of 
compounds, for instance 0 300, 300 1000, 1000 1500 and 1500 2134.
Only if run in smaller chunks: merge the results:
```shell
python merge_fda.py [local path]
```

Expected output:
`data/v1_compounds_FDA.csv`

### Determine ADR occurrences and associations for all compounds 
Run in shell
```shell
python openFDAoutput2modelformat.py [local path]
```
Expected output:  
`data/v1_compounds_FDA_model_format_SOC_ocr_prob.csv`  
`data/v1_compounds_FDA_model_format_SOC_ocr_bool.csv`  
`data/v1_compounds_FDA_model_format_HLGT_ocr_prob.csv`  
`data/v1_compounds_FDA_model_format_HLGT_ocr_bool.csv`  

For a demo version that describes the processing steps, 
see `OpenFDAoutput2modelformat_demo.ipynb` 

### Run random Forest models:
Run in R (or RStudio):
```shell
R_code_modeling.R
```

Please make sure that SOC/ and HLGT/ folders are in the same folder with this
R code (R_code_modeling.R)
Follow the instructions in the R script (e.g. please specify/change the working 
directory/folder on your local computer in scripts where indicated).

The 2 files provided, baseMatrix and v1_compounds_FDA_model_format_*_ocr_bool.csv, 
where * is either SOC or HLGT, are the inputs for the R script.

Please note that the R code contains references to baseMatrix.Rdata.  This is a 
matrix containing compound identifiers on rows (cid1, cid2, cid3, ...) and assay 
labels on columns.  In the work reported in the paper, this was obtained from 
Novartis in vitro assay data (Supplementary Table 2).  The column 
headers in these data matrices contain assay IDs (A1, A2, A3, 
for assays 1, 2, 3, etc. discretized into ranges as described in methods (“A1_CL1”, 
“A1_CL2”, “A1_CL3”, “A2_CL1”, “A2_CL2”, “A2_CL3”, etc). The mapping of assays number 
to targets can be found in `mergedassaynumber-to-onewordname.csv` and Supplemental 
Table 2. For brevity, while the 
code shows analyses on the redefined BSEP ranges, we excluded these sample input and 
output files from the example data. The redefined BSEP data still can be found in 
Supplemental table 2.
Likewise, for brevity, we excluded the sample input and output for the analysis using 
a different random number generator seed. But the different seed results can be generated 
by copying the available input files in their respective seed folders and running the 
R code accordingly. Please see the R code for further instructions.

With the latest version of utiml package, their accuracy formula was changed and 
so the summary metrics results were not consistent. If you use the latest utiml 
version as well, please use the custom code to use the standard accuracy formula 
(defined in our methods) starting on line 303.

Expected output:
 
`data/Features_ADRs_predictions.csv`  
`data/summ_metrics.csv`  
`data/summ_metrics_V2.csv` (if the code starting at line 296 is run)  
`data/summ_metrics_ADRs.csv` (for all models except in SOC)

### Target-ADR associations
The ADR target association code is provided:
```shell
ADR_target.ipynb
```
 

### Systematic Pubmed api query for validation of associations
Obtain an NBCI api key (follow instructions on NCBI website:
https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/).  
Set api key and local paths in the scripts below and run in R (or Rstudio):
```shell
pubmed_query_mesh_both_all.R    
pubmed_query_mesh_both_positive.R
```
Stastistics on the results are described in:
```shell
Pubmed_search_validation_statistics_rev.ipynb  
Pubmed_search_validation_statistics_negpairs_from_pos_rev.ipynb
```

### Run times
On a "normal" computer, OpenFDA query: with API-key (see details above) for OpenFDA:
up to 2 hours. Without API-key the query speed is capped and takes longer.

Random Forest training and cross validation: between 1 hour to 24 hours (depending 
on SOC and HLGT with 1 core, i.e. no parallel computing).

Target-ADR associations: less than 1h.

pubmed queries: ~24h.

### Instructions for user provided data
A list of unique identifiers and generic compound names in csv/tsv format (as a replacement 
for `data/compounds.csv`) is required to run the openFDA query.  
For the Random forest model, discretize and one hot encoded the pharmacology data as 
described in the manuscript methods and the R_code_modeling.R script.  
For the target-ADR associations: provide the mapping file between assay ID and gene names.
 
### Citation
Robert Ietswaart<sup>\*,#</sup>, Seda Arat<sup>\*,#</sup>, Amanda X. Chen<sup>\*</sup>, 
Saman Farahmand<sup>\*</sup>, Bumjun Kim, William DuMouchel, 
Duncan Armstrong, Alexander Fekete, Jeffrey J. Sutherland<sup>#</sup>, Laszlo Urban<sup>#</sup>  
*Machine learning guided association of adverse drug reactions with in vitro target-based 
pharmacology* (2019), [BioRxiv; 750950](https://www.biorxiv.org/content/10.1101/750950v2).


