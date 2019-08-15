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
- Download/install python 3.7 from python.org
- In shell:
```shell
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
```
R/Rstudio:

