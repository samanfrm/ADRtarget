#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import pandas as pd

v1_1=pd.read_csv('Compound_FDA_V1_from0_to_299.csv',sep='\t')
v1_2=pd.read_csv('Compound_FDA_V1_from300_to_999.csv',sep='\t')
v1_3=pd.read_csv('Compound_FDA_V1_from1000_to_1499.csv',sep='\t')
v1_4=pd.read_csv('Compound_FDA_V1_from1500_to_2133.csv',sep='\t')

v1_concat=[v1_1,v1_2,v1_3,v1_4]
v1_concat=pd.concat(v1_concat,ignore_index=True)
v1_concat.to_csv('v1_compounds_FDA_demo.csv',sep='\t')

