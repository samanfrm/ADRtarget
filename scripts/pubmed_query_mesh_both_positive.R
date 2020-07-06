#Saman Farahmand - GNU General Public License v3.0

### Citation:
#Robert Ietswaart<sup>\*,#</sup>, Seda Arat<sup>\*,#</sup>, Amanda X. Chen<sup>\*</sup>, 
#Saman Farahmand<sup>\*</sup>, Bumjun Kim, William DuMouchel, 
#Duncan Armstrong, Alexander Fekete, Jeffrey J. Sutherland<sup>#</sup>, Laszlo Urban<sup>#</sup>  
#*Machine learning guided association of adverse drug reactions with in vitro target-based 
#pharmacology*, Ebiomedicine (2020) <https://doi.org/10.1016/j.ebiom.2020.102837>.

library(rentrez)
library(org.Hs.eg.db)
library(dplyr)
library(tidyr)
require(biomaRt)

nkey="XXXXXXXX" #replace with ncbi api key
dir_path="path_to_ADRtarget_directory" #set the working directory path point to 'PATH_to/ADRtarget/' main directory
setwd(dir_path)
positive=read.csv(file = "data/predicted_ADR_target_final.txt",header = T,sep = '\t')#replace with local paths
ADR2mesh=read.csv(file = "data/predicted_meddra2mesh.txt",header = T,sep = '\t')#replace with local paths
gene2mesh=read.csv(file = "data/predicted_gene2mesh.txt",header = T, sep = '\t')#replace with local paths



positive=positive %>% mutate(meddra_name=MedDRA.term.name,meddra_id=MedDRA.ID,gene_symbol=Entrez.Gene.Symbol.for.target,target_mesh=MeSH.term.for.target)
positive_ADR=positive %>% dplyr::select(c("meddra_name","meddra_id")) %>% distinct()
positive_ADR=left_join(positive_ADR,ADR2mesh,by=c("meddra_id"="meddra_code")) %>% dplyr::select(c("meddra_name","meddra_id","mesh_term")) %>% distinct()
positive_ADR=positive_ADR %>% na.omit()

results_ADR=data.frame("HLGT"='',"mesh_term"='',"mesh_N"='',"mesh_pmids"='')

for (i in c(1:nrow(positive_ADR))) {
  ADR=as.character(positive_ADR[i,"meddra_name"])
  ADR_mesh=as.character(positive_ADR[i,"mesh_term"])
  query_term=paste(ADR_mesh,' [MeSH Terms]')
  dbs = entrez_search(db="pubmed", term=query_term,retmax="1000000000",api_key=nkey)
  if(length(dbs$ids)){
    results_ADR=rbind(results_ADR,data.frame("HLGT"=ADR,"mesh_term"=ADR_mesh,
                                     "mesh_N"=as.factor(length(dbs$ids)),
                                     "mesh_pmids"=paste(dbs$ids,collapse = ";")))
    
  }
  else{
    results_ADR=rbind(results_ADR,data.frame("HLGT"=ADR,"mesh_term"=ADR_mesh,
                                     "mesh_N"=as.factor(0),
                                     "mesh_pmids"=''))
  }
  print(paste0("Writing the results for ADR ",i,"_",length(dbs$ids)))
}
print("Writing the final results...")
results_ADR=results_ADR[-1,]
################ Summarise mesh terms for each unique ADR/HLGT

unique_HLGT=unique(results_ADR$HLGT)
results_ADR_uniq=data.frame("HLGT"='',"HLGT_N"='',"HLGT_pmids"='')
for (i in c(1:length(unique_HLGT))){
  temp=results_ADR %>% filter(HLGT==unique_HLGT[i])
  ccl=''
  for (j in c(1:nrow(temp))){
    arr=unlist(strsplit(x = as.character(temp[j,"mesh_pmids"]),split = ';'))
    ccl=c(ccl,arr)
  }
  ccl=unique(ccl[ccl!=''])
  results_ADR_uniq=rbind(results_ADR_uniq,data.frame("HLGT"=unique_HLGT[i],"HLGT_N"=as.factor(length(ccl)),
                                   "HLGT_pmids"=paste(ccl,collapse = ';'))
  )
}

results_ADR_uniq=results_ADR_uniq[-1,]
##########################  get PMIDs for the mesh terms of genes

positive_gene=positive %>% dplyr::select(c("gene_symbol","target_mesh")) %>% distinct()

results_gene=data.frame("gene"='',"mesh_tr"='',"mesh_N"='',"mesh_pmids"='')

for (i in c(1:nrow(positive_gene))) {
  gene=as.character(positive_gene[i,"gene_symbol"])
  gene_mesh=as.character(positive_gene[i,"target_mesh"])
  query_term=paste(gene_mesh,' [MeSH Terms]')
  dbs = entrez_search(db="pubmed", term=query_term,retmax="1000000000",api_key=nkey)
  if(length(dbs$ids)){
    results_gene=rbind(results_gene,data.frame("gene"=gene,"mesh_tr"=gene_mesh,
                                             "mesh_N"=as.factor(length(dbs$ids)),
                                             "mesh_pmids"=paste(dbs$ids,collapse = ";")))
  }
  else{
    results_gene=rbind(results_gene,data.frame("gene"=gene,"mesh_tr"=gene_mesh,
                                             "mesh_N"=as.factor(0),
                                             "mesh_pmids"=''))
  }
  print(paste0("Writing the results for gene mesh ",i,"_",length(dbs$ids)))
}
print("Writing the final results...")
results_gene=results_gene[-1,]

#################### get genes and ADRs results together and Do intersection

#### Do random crossing ADRs and genes
full_ADR_target=tidyr::crossing(results_ADR_uniq,results_gene) 
colnames(full_ADR_target)

############## Do intersection


get_inters=function(x,y){
  x=unlist(strsplit(as.character(x),';'))
  y=unlist(strsplit(as.character(y),';'))
  intr=intersect(x,y)
  return (paste(intr,collapse = ';'))
}
get_intrs_num=function(x,y){
  x=unlist(strsplit(as.character(x),';'))
  y=unlist(strsplit(as.character(y),';'))
  intr=intersect(x,y)
  return (length(intr))
}

#positive_results$intersect_N=''
full_ADR_target$intersect_N=''

for (i in c(1:nrow(full_ADR_target))){
  full_ADR_target[i,"intersect_N"]=get_intrs_num(as.character(full_ADR_target[i,"HLGT_pmids"]),
                                           as.character(full_ADR_target[i,"mesh_pmids"]))
}
full_ADR_target_summ=full_ADR_target %>% dplyr::select(c("HLGT","gene","mesh_tr","HLGT_N","mesh_N","intersect_N"))

##tag positive with True
control=positive %>% dplyr::select(c("meddra_name","gene_symbol"))
control$flag=TRUE
## join
full_ADR_target_summ=left_join(full_ADR_target_summ,control,by=c("HLGT"="meddra_name","gene"="gene_symbol"))
full_ADR_target_summ$flag[is.na(full_ADR_target_summ$flag)]=FALSE


write.table(x = full_ADR_target_summ,file = "data/Full_result_ADRmesh_Genemesh_pubmed_NEW.csv",sep = '\t',row.names = F)#replace with local paths
