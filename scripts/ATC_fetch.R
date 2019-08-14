#install.packages('httr')
library(httr)
#install.packages("RCurl")
library(RCurl)
#install.packages("RJSONIO")
library(RJSONIO)
#install.packages("plyr")
library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)

getATC=function(name){
  return (tryCatch({
  name=gsub(" ","%20",name)
  pubchem_url_in = paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/',name,'/JSON')
  pubchempage = GET(pubchem_url_in)
  if(pubchempage[[2]]!=200){
    return(pubchempage[[2]])
  }
  page_text = content(pubchempage,as='text')
  page_test2 = fromJSON(page_text)
  cid=as.character(unlist(page_test2)[1])
  
  pubchem_url2 = paste('https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/',cid,'/JSON',sep='')
  pubchempage2 = GET(pubchem_url2) 
  page_text = content(pubchempage2,as='text')
  page_test2 = fromJSON(page_text)
  ATC_parse1  = grep('www.whocc.no',unlist(page_test2),value=T)
  ATC_parse2 = ATC_parse1[2]
  ATCout1 = strsplit(strsplit(ATC_parse2,'code=')[[1]][2],'&showdescription')[[1]][1]
  return(ATCout1)
  },error = function(e){
    return (as.character(e))
  }))
}

MDDB=read.delim("~/Desktop/ATC_mining/MDDB_clean.csv",sep=',')
pname=read.delim("~/Desktop/ATC_mining/MDDB_pname.csv",sep=',')
MDDB=left_join(MDDB,pname,"number")

MDDB_ATC= MDDB %>% filter(ATC_code!='' & !is.na(ATC_code))
MDDB_ATC=MDDB_ATC %>% select(number,compound,ATC=ATC_code,ACT_pname)

MDDB_noATC=MDDB %>% filter(ATC_code=='' | is.na(ATC_code))
MDDB_noATC=MDDB_noATC %>% select(number,compound)


MDDB_noATC=MDDB_noATC %>% separate(compound,into = paste0("c",c(1:10)) ,sep = ';',remove = T)
MDDB_noATC=MDDB_noATC %>% melt(paste0("c",c(1:10)),id.vars=1) %>% na.omit() 
MDDB_noATC2=MDDB_noATC %>% select(number,compound=value) 
MDDB_noATC2$compound=gsub("[^[:alnum:] ]","",MDDB_noATC2$compound)
MDDB_noATC2$compound=gsub("[0-9]+","",MDDB_noATC2$compound)
MDDB_noATC2=MDDB_noATC2 %>% mutate(ATC=unlist(lapply(compound, getATC)))
MDDB_noATC2=left_join(MDDB_noATC2,MDDB,by="number")
MDDB_noATC2=MDDB_noATC2 %>% select(number,compound=compound.x,ATC,ACT_pname) 

MDDB_matched=MDDB_noATC2 %>% filter(ATC!='404' & ATC!="400" & !is.na(ATC))
MDDB_no_matched=MDDB_noATC2 %>% filter(ATC=='404' | is.na(ATC) | ATC=="400")
MDDB_no_matched=MDDB_no_matched %>% mutate(ATC=unlist(lapply(ACT_pname, getATC)))
MDDB_all=rbind(rbind(MDDB_matched,MDDB_no_matched),MDDB_ATC)

MDDB_all=MDDB_all %>% mutate(ATC=replace(ATC,ATC=="400",NA))
MDDB_all=MDDB_all %>% mutate(ATC=replace(ATC,ATC=="404",NA))
MDDB_all=MDDB_all %>% mutate(ATC=replace(ATC,ATC=="503",NA))
MDDB_all=MDDB_all %>% mutate(ATC=replace(ATC,is.na(ATC),""))
MDDB_all=MDDB_all[!colnames(MDDB_all) %in% c("ACT_pname","compound")]
nrow(filter(MDDB_all,ATC==""))
MDDB_all$ATC=gsub(";;",";",MDDB_all$ATC)

MDDB_all=MDDB_all %>% separate(ATC,into = paste0("c",c(1:11)) ,sep = ';',remove = T)
MDDB_all=MDDB_all %>% melt(paste0("c",c(1:11)),id.vars=1) %>% na.omit() 
nrow(filter(MDDB_all,value==""))

MDDB_all=MDDB_all %>% extract(value,c("1st"),"(^.{1})",remove = F)
MDDB_all=MDDB_all %>% extract(value,c("2st"),"(^.{1,3})",remove = F)


ATC_names=read.delim("~/Desktop/ATC_mining/ATC_names.txt",sep=",",header = 0)

MDDB_all_with_names=left_join(MDDB_all,ATC_names,by=c("1st"="V1"))
MDDB_all_with_names=left_join(MDDB_all_with_names,ATC_names,by=c("2st"="V1"))
write.csv(x = MDDB_all_with_names,"~/Desktop/ATC_mining/MDDB_all_names.csv")