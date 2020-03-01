#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Saman Farahmand / Robert Ietswaart
# MIT license

### Citation:
#Robert Ietswaart<sup>\*,#</sup>, Seda Arat<sup>\*,#</sup>, Amanda X. Chen<sup>\*</sup>, 
#Saman Farahmand<sup>\*</sup>, Bumjun Kim, William DuMouchel, 
#Duncan Armstrong, Alexander Fekete, Jeffrey J. Sutherland<sup>#</sup>, Laszlo Urban<sup>#</sup>  
#*Machine learning guided association of adverse drug reactions with in vitro target-based 
#pharmacology* (2019), [BioRxiv; 750950](https://www.biorxiv.org/content/10.1101/750950v2).



import pandas as pd
import os

repo_path = str(sys.argv[1])
path = os.path.join(repo_path,'data')

v1_1=pd.read_csv(os.path.join(path,'Compound_FDA_V1_from0_to_299.csv'),sep='\t')
v1_2=pd.read_csv(os.path.join(path,'Compound_FDA_V1_from300_to_999.csv'),sep='\t')
v1_3=pd.read_csv(os.path.join(path,'Compound_FDA_V1_from1000_to_1499.csv'),sep='\t')
v1_4=pd.read_csv(os.path.join(path,'Compound_FDA_V1_from1500_to_2133.csv'),sep='\t')

v1_concat=[v1_1,v1_2,v1_3,v1_4]
v1_concat=pd.concat(v1_concat,ignore_index=True)
v1_concat.to_csv('v1_compounds_FDA.csv',sep='\t')

