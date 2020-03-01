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


import requests
import pandas as pd
import jellyfish
import string
import sys
import os

def fda_process(url,gname,HLGT_dict,SOC_dict,all_reports):
            resp = requests.get(url=url)
            data = resp.json()
            for result in data['results']:
                    try:
                        if('primarysource' in result.keys()):
                            if('qualification' in result['primarysource'].keys()):
                                if(result['primarysource']['qualification']=='1'):
                                    for drug in result['patient']['drug']:
                                            if('drugcharacterization' in drug.keys()):
                                                 if(drug['drugcharacterization']=='1'):
                                                    exclude = set(string.punctuation)
                                                    gname = ''.join(ch for ch in gname if ch not in exclude)
                                                    gname=gname.strip().upper()
                                                    flag=0
                                                    if('activesubstance' in drug.keys()):
                                                        drugname=drug['activesubstance']['activesubstancename'].upper()
                                                        score=jellyfish.jaro_distance(gname,drugname)
                                                        if(score>= 0.8):
                                                            flag=1
                                                    if('medicinalproduct' in drug.keys() and flag==0):
                                                        drugname=drug['medicinalproduct'].upper()
                                                        score=jellyfish.jaro_distance(gname,drugname)
                                                        if(score>= 0.8):
                                                            flag=1
                                                    else:
                                                        continue

                                    if(flag==1):
                                        reportid=result['safetyreportid']
                                        PTs=[]
                                        HGLTs=[]
                                        SOCs=[]
                                        for reaction in result['patient']['reaction']:
                                           if('reactionmeddrapt' in reaction.keys()):
                                               PT=reaction['reactionmeddrapt'].lower()
                                               PTs.append(PT)
                                               mapped_hglt=map_adr_to_meddra(PT, HLGT_dict)
                                               HGLTs.append(mapped_hglt)
                                               mapped_soc=map_adr_to_meddra(PT, SOC_dict)
                                               SOCs.append(mapped_soc)
                                        all_reports[reportid]={'PTs':PTs,'HGLTs':HGLTs,'SOCs':SOCs}

                    except KeyError:
                        continue
            return all_reports

def create_output(reports):
    len_report=str(len(reports))
    reports_ids=[]
    len_PTs=[]
    PTs=[]
    len_HGLTs=[]
    HGLTs=[]
    len_SOCs=[]
    SOCs=[]
    for key,value in reports.items():
        reports_ids.append(key)
        if(len(value['PTs'])>0):
            len_PTs.append(str(len(value['PTs'])))
            PTs.append(';'.join(value['PTs']))
        else:
            len_PTs.append('0')
            PTs.append('NaN')

        if(len(value['HGLTs'])>0):
            len_HGLTs.append(str(len(value['HGLTs'])))
            HGLTs.append(';'.join(value['HGLTs']))
        else:
            len_HGLTs.append('0')
            HGLTs.append('NaN')

        if(len(value['SOCs'])>0):
            len_SOCs.append(str(len(value['SOCs'])))
            SOCs.append(';'.join(value['SOCs']))
        else:
            len_SOCs.append('0')
            SOCs.append('NaN')
    return len_report,'|||'.join(reports_ids),'|||'.join(len_PTs),'|||'.join(PTs),'|||'.join(len_HGLTs),'|||'.join(HGLTs),'|||'.join(len_SOCs),'|||'.join(SOCs)

def map_adr_to_meddra(anentry, meddradict):

    try:
        result=meddradict[str(anentry.lower())]
        return result
    except KeyError:
        return 'NaN'

    

ind_start = int(sys.argv[1])
ind_end = int(sys.argv[2])
key = str(sys.argv[3])
repo_path = str(sys.argv[4])
path = os.path.join(repo_path,'data')


filename = 'compounds.csv'### two column file: first with unique compound ids and second with generic gene names used for the openFDA queries
generics = pd.read_csv(os.path.join(path,filename),header=0,names=['number','generic'])

##Create dictionary for HLGT and SOC
filename = 'meddr1a_full_12092006.xlsx'
meddra = pd.read_excel(os.path.join(path,filename))

HLGT_dict = {}
for index, arow in meddra.iterrows():
    pt = arow['PT_TXT'].lower()
    try:
        HLGT_dict[pt]
    except KeyError:
        HLGT_dict[str(pt)] = str(arow['HLGT_TXT']).lower()

SOC_dict = {}
for index, arow in meddra.iterrows():
    pt = arow['PT_TXT'].lower()
    try:
        SOC_dict[pt]
    except KeyError:
        SOC_dict[str(pt)] = str(arow['SOC_TXT']).lower()


limit = 100
qualification = 1
drug_characterization = 1

compound_ADRs = pd.DataFrame(columns=['number','name','#reports','reports_id','PTs_per_report','HGLT_per_report','SOC_per_report'])


for j in range(ind_start,ind_end):
    comp=generics.loc[j,'generic']
    comp_id=generics.loc[j,'number']
    gnames=comp.split(';')
    print(str(j) + ' out of '+ str(ind_end))
    for gname in gnames:
        reports={}
        gname=gname.strip()
        gname2=gname.replace(' ','+')
        print(comp_id +' - '+ gname)
        url='https://api.fda.gov/drug/event.json?api_key='+key+'&search="'+gname2+'"+AND+primarysource.qualification:"1"&limit='+str(limit)
        resp = requests.get(url=url)
        data = resp.json()
        if('error' in data.keys()):
            compound_ADRs=compound_ADRs.append({'number':comp_id,'name':gname,'#reports':'error','reports_id':data['error']['message'],'PTs_per_report':'','HGLT_per_report':'','SOC_per_report':''},ignore_index=True)
        else:
            total=data['meta']['results']['total']
            count=0
            if(total>100):
                count=total//100
            skipped=0
            for i in range(1,count+1):
                url='https://api.fda.gov/drug/event.json?api_key='+key+'&search="'+gname2+'"+AND+primarysource.qualification:"1"&limit='+str(limit)+"&skip="+str(skipped)
                reports=fda_process(url,gname,HLGT_dict,SOC_dict,reports)
                skipped=i*100

            if(total>skipped): #for the last part
                url='https://api.fda.gov/drug/event.json?api_key='+key+'&search="'+gname2+'"+AND+primarysource.qualification:"1"&limit='+str(limit)+"&skip="+str(skipped)
                reports=fda_process(url,gname,HLGT_dict,SOC_dict,reports)

            if(len(reports)>0):
                len_report,reports_ids,len_PTs,PTs,len_HGLTs,HGLTs,len_SOCs,SOCs=create_output(reports)
                compound_ADRs=compound_ADRs.append({'number':comp_id,'name':gname,'#reports':len_report,'reports_id':reports_ids,'PTs_per_report':PTs,'HGLT_per_report':HGLTs,'SOC_per_report':SOCs},ignore_index=True)
            else:
                compound_ADRs=compound_ADRs.append({'number':comp_id,'name':gname,'#reports':'0','reports_id':'','PTs_per_report':'','HGLT_per_report':'','SOC_per_report':''},ignore_index=True)

if ind_start == 0 and ind_end == 2134:
    filename = 'v1_compounds_FDA.csv'
else:
    filename = 'Compound_FDA_V1_from'+str(ind_start)+'_to_'+str(ind_end-1)+'.csv'
    
compound_ADRs.to_csv(os.path.join(path,filename),sep='\t')



