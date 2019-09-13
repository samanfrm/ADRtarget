
library(rentrez)
library(org.Hs.eg.db)
library(dplyr)
library(tidyr)
require(biomaRt)

nkey="485af51aa11dbfc71be592e6da677c416b09"
ADR2mesh=read.csv(file = "~/Desktop/pubmed_novartis/ADR_mesh_list_all.txt",header = T,sep = '\t')
ADR2mesh=ADR2mesh %>% distinct()
ADR_dic=read.csv(file="~/Desktop/pubmed_novartis/ meddr1a_full_12092006.csv",header = T,sep=',')
ADR2mesh=left_join(ADR2mesh,ADR_dic,by=c("meddra_code"="HLGT_Code")) %>% 
  dplyr::select(c('meddra_code',"mesh_term",'HLGT_TXT')) %>% distinct()


results_ADR=data.frame("HLGT"='',"mesh_term"='',"mesh_N"='',"mesh_pmids"='')

for (i in c(1:nrow(ADR2mesh))) {
  ADR=as.character(ADR2mesh[i,"HLGT_TXT"])
  ADR_mesh=as.character(ADR2mesh[i,"mesh_term"])
  query_term=paste(ADR_mesh,' [MeSH Terms]')
  dbs = entrez_search(db="pubmed", term=query_term,retmax="1000000000",api_key=nkey)
  if(length(dbs$ids)){
    results_ADR=rbind(results_ADR,data.frame("HLGT"=ADR,"mesh_term"=ADR_mesh,
                                             "mesh_N"=as.factor(length(dbs$ids)),
                                             "mesh_pmids"=paste(dbs$ids,collapse = ";")))
    
    # results[i,"gene_N"]=length(dbs$ids)
    # results[i,"gene_pmids"]=paste(dbs$ids,collapse = ";")
  }
  else{
    results_ADR=rbind(results_ADR,data.frame("HLGT"=ADR,"mesh_term"=ADR_mesh,
                                             "mesh_N"=as.factor(0),
                                             "mesh_pmids"=''))
    # results[i,"gene_N"]=0
    # results[i,"gene_pmids"]=''
  }
  print(paste0("Writing the results for ADR ",i,"_",length(dbs$ids)))
}
print("Writing the final results...")
results_ADR=results_ADR[-1,]

#results=results %>% mutate(HLGT_N=NO_PMIDs,HLGT_pmids=PMIDs) %>% dplyr::select(-c(NO_PMIDs,PMIDs))

################ Summarise mesh terms for each unique ADR/HLGT

unique_HLGT=unique(results_ADR$HLGT)
results_ADR_uniq=data.frame("HLGT"='',"HLGT_N"='',"HLGT_pmids"='')
for (i in c(1:length(unique_HLGT))){
  temp=results_ADR %>% filter(HLGT==unique_HLGT[i])
  ccl=''
  for (j in c(1:nrow(temp))) {
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

gene2mesh=read.csv(file = "~/Desktop/pubmed_novartis/final_gene_mesh_all.txt",header = T,sep='\t')
gene2mesh=gene2mesh %>% dplyr::select(c("gene_symbol","mesh_term")) %>% distinct()

results_gene=data.frame("gene"='',"mesh_tr"='',"mesh_N"='',"mesh_pmids"='')

for (i in c(1:nrow(gene2mesh))) {
  gene=as.character(gene2mesh[i,"gene_symbol"])
  gene_mesh=as.character(gene2mesh[i,"mesh_term"])
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
full_ADR_target_summ$HLGT=tolower(full_ADR_target_summ$HLGT)

##tag positive with True
positive=read.csv(file = "~/Desktop/pubmed_novartis/predicted_ADR_target_final.txt",header = T,sep = '\t')
positive=positive %>% mutate(meddra_name=MedDRA.term.name,meddra_id=MedDRA.ID,gene_symbol=Entrez.Gene.Symbol.for.target,target_mesh=MeSH.term.for.target)
positive=positive %>% dplyr::select(c('meddra_name','gene_symbol')) 
positive$flag=TRUE
## join
full_ADR_target_summ=left_join(full_ADR_target_summ,positive,by=c("HLGT"="meddra_name","gene"="gene_symbol"))
full_ADR_target_summ$flag[is.na(full_ADR_target_summ$flag)]=FALSE


write.table(x = full_ADR_target_summ,file = "~/Desktop/Full_result_ADRmesh_Genemesh_pubmed_NEW_106.csv",sep = '\t',row.names = F)


