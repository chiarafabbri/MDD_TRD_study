### R code to extract diagnoses and prescription primary care records for depression, psychotic, bipolar and substance use disorders and antidepressant prescriptions
### using the UK Biobank resource (https://www.ukbiobank.ac.uk) or similar electronic health records of primary care
### the code may need to be adapted (e.g. file names, file formats etc)

### information on primary care data can be found in Resource 591 (http://biobank.ndph.ox.ac.uk/showcase/refer.cgi?id=591)
### Information on codes used in the primary care data can be found in Resource 592 (http://biobank.ndph.ox.ac.uk/showcase/refer.cgi?id=592)

### this script was used for the analyses of the study available at: https://doi.org/10.1101/2020.08.24.20178715

################################################
########## DIAGNOSIS DATA EXTRACTION ###########
################################################

## R version 3.6.0 was used for the study described at https://doi.org/10.1101/2020.08.24.20178715 

# load R packages (you may need to install these packages)
library(tidyverse)
library(lubridate)
library(readxl)
library(dplyr)
library(knitr)
library(xtable)
library(data.table)
library(ggplot2)
library(readr)

# read clinical data:
clinical <- fread("gp_diagnoses",header=T, data.table=F, na.strings="") # this is the table with diagnostic records that can be downloaded from UK biobank website (https://www.ukbiobank.ac.uk) by investigators having access to the data (field ID 42040; note that the name of the file as you downloaded is probably different)

# open Read2 and Read3 code files (these include all codes for depression, psychosis, bipolar disorder and substance abuse disorders)
readv2 <- fread("readv2_mood_psych_sub",header=T,data.table=F) 
# text table with the following columns: readv2_code, term_v2, type, readv3_code and term_v3 for the diagnoses of interest, where type is the type of disorder (e.g. depression, psychosis); these information can be extracted using documentation provided at: http://biobank.ndph.ox.ac.uk/showcase/showcase/auxdata/primarycare_codings.zip 
# alternatively, see https://doi.org/10.1101/2020.08.24.20178715 and/or contact chiara.fabbri@yahoo.it to ask for a copy of this table used in: doi: https://doi.org/10.1101/2020.08.24.20178715

readv3 <- fread("readv3_mood_psych_sub",header=T,data.table=F) 
# text table with the following columns: readv3_code, term_v3 and type for the diagnoses of interest; see comments reported for readv2 


# separately extract read 2 and then map these to read3 code before extracting read3 and binding rows

# extract read 2 codes and add mapped read3 codes and terms
clinical %>%
  filter(clinical$read_2 %in% readv2$readv2_code) %>%
  left_join(readv2, by=c("read_2"="readv2_code")) %>%
  select(-read_3,-value1,-value2,-value3) %>%
  rename(read_3 = readv3_code) -> v2_psych_codes

# extract read 3 code and add read2 codes
clinical %>%
  filter(clinical$read_3 %in% readv3$readv3_code) %>%
  left_join(readv3, by=c("read_3"="readv3_code")) %>%
  select(-value1,-value2,-value3) %>%
  bind_rows(v2_psych_codes) %>%
  filter(!is.na(event_dt)) -> psych_codes

psych_codes$event_dt <- as.Date(psych_codes$event_dt, "%d/%m/%Y")

# check and set to missing the following dates (see http://biobank.ndph.ox.ac.uk/showcase/showcase/docs/primary_care_data.pdf):
# where clinical event or prescription date precedes participant date of birth it has been altered to 01/01/1901
# Where the date matches participant date of birth it has been altered to 02/02/1902
# Where the date follows participant date of birth but is in the year of their birth it has been altered to 03/03/1903
# Where the date was in the future this has been changed to 07/07/2037 as these are likely to have been entered as a place-holder or other system default

dates_rm<-which(psych_codes$event_dt=="1901-01-01"|psych_codes$event_dt=="1902-02-02"|psych_codes$event_dt=="1903-03-03"|psych_codes$event_dt=="2037-07-07") 
psych_codes<-psych_codes[-dates_rm,]

# mark date but keep record for all rows where the date is 01/01 (this is likely some issue with the system defaulting to 01/01, but the year may be correct):
psych_codes$diagn_month_day <- format(as.Date(psych_codes$event_dt), "%m-%d")
table(psych_codes$diagn_month_day=='01-01')
psych_codes$diagn_year <- format(as.Date(psych_codes$event_dt),"%Y")
psych_codes$diagn_year <- as.numeric(psych_codes$diagn_year)


# write output
fwrite(psych_codes,"mood_BP_SU_codes.txt",col.names=T,row.names=F,quote=F,sep="\t")


################################################
######### PRESCRIPTION DATA EXTRACTION #########
################################################

## This part of the code will extract medication codes of interest
## Read 2, BNF, dm+d

# open full medication data (table with prescription records (field ID 42039) that can be downloaded from UK biobank website (https://www.ukbiobank.ac.uk) by investigators having access to the data, check file name since it is probably different from the one in this script)
meds <-
  fread("gp_prescriptions",
    header = T,
    data.table = F,
    colClasses = c(
      "read_2" = "character",
      "bnf_code" = "character",
      "dmd_code" = "character"), 
    na.strings="")

#####################
#### BNF codes ######
#####################

## extract all BNF codes (this will also include part of read2 codes, as some individuals will have both codes)
bnf_codes <- c(
  "0403","04.03")

bnf_data <- data.frame()
for (i in bnf_codes){
  print(paste("started BNF code",i, "at", Sys.time()))
  temp <- grep(paste0("^",i,"."),meds$bnf_code,ignore.case=T)
  temp2 <- meds[temp,]
  bnf_data <- rbind(bnf_data,temp2)
  print(dim(bnf_data))
  rm(temp2)
  }

# remove rows that have 040300 as this is an anti-epileptic not antidepressant (pregabalin)
pregabalin_04030000=which(bnf_data$bnf_code=="04030000")
bnf_data<-bnf_data[-pregabalin_04030000,]

# identify number of duplicates
bnf_data$unique <- !(duplicated(bnf_data) | duplicated(bnf_data, fromLast = TRUE))
table(bnf_data$unique)

# remove duplicates that have been added due to the way we extract the data
bnf_data %>%
select(-unique) %>%
distinct() -> bnf_nodup

#######################
#### Read2 codes ######
#######################

# load clinical read2 codes:
read_2_codes <- read.table("ADs_read_2_codes",header=T) 
# this corresponds to a text table with the following column: read2_code that includes read2 codes for antidepressant medications
# these codes can be extracted using documentation provided at: http://biobank.ndph.ox.ac.uk/showcase/showcase/auxdata/primarycare_codings.zip 
# alternatively, see https://doi.org/10.1101/2020.08.24.20178715 and/or contact chiara.fabbri@yahoo.it to ask for a copy of this table used in: doi: https://doi.org/10.1101/2020.08.24.20178715


# filter medication file to include read2 codes but no people with read2 and bnf codes as these have already been extracted
meds %>%
filter(!is.na(read_2) & is.na(bnf_code)) -> meds_no_bnf

# extract read 2 codes from medication data
read_2_data <- meds_no_bnf[which(meds_no_bnf$read_2 %in% read_2_codes$read2_code),]

# this identifies all duplicates, not just one of duplicate pair
read_2_data$unique <- !(duplicated(read_2_data) | duplicated(read_2_data, fromLast = TRUE))
table(read_2_data$unique)

# remove duplicates
read_2_data %>%
select(-unique) %>%
distinct() -> read_2_nodup


#######################
##### dm+d codes ######
#######################


# medication names : 
dmd_names <-c("agomelatine","allegron","alventa","amfebutamone","amitriptyline","amoxapine","amphero","anafranil","apclaven","asendis","bolvidon","bonilux","brintellix","bupropion","butriptyline","cipralex","cipramil","citalopram","clomipramine","cymbalta","dapoxetine","defanyl","depalta","depefex","desipramine","domical","dosulepin","dothapax","dothiepin","doxepin","duciltia","duloxetine","dutonin","dutor","edronax","efexor","escitalopram","faverin","felicium","feprapax","fetzima","fluanxol","fluoxetine","flupenthixol","flupentixol","fluvoxamine","gamanil","imipramine","iprindole","iproniazid","isocarboxazid","lentizol","levomilnacipran","lofepramine","loferpramine","lomont","ludiomil","lustral","majoven","manerix","maprotiline","marplan","mianserin","milnacipran","mirtazapine","mirtazepine","mirtazipine","moclobemide","molipaxin","motipress","motival","nardil","nefazodone","nortriptyline","norval","olena","optimax","oxactin","parnate","paroxetine","parstelin","perphenazine","phenelzine","politid","prepadine","priligy","prothiaden","protriptyline","prozac","prozep","ranfaxine","ranflutin","reboxetine","rodomel","seroxat","sertraline","sinepin","sinequan","sunveniz","surmontil","tardcaps","thaden","tifaxin","tofranil","tonpular","tranylcypromine","trazadone","trazodone","trimipramine","triptafen","triptafen-m","trixat","tryptizol","tryptophan","valdoxan","vaxalin","venaxx","vencarm","venlablue","venladex","venlafaxin","venlafaxine","venlalic","venlaneo","venlasov","vensir","venzip","vexarin","viepax","viloxazine","vivactil","vivalan","vortioxetine","winfex","yentreve","zispin","zyban")

# subset to only have dmd codes otherwise I'll extract a lot of read 2 individuals as well (this is purely to have smaller file for now)
meds %>%
filter(!is.na(dmd_code) & is.na(read_2)) -> meds_dmd_only

dmd_data <- data.frame()
for (i in dmd_names){
  print(paste("started dm+d code",i, "at", Sys.time()))
  temp <- grep(paste0("^",i,"."),meds_dmd_only$drug_name,ignore.case=T)
  temp2 <- meds_dmd_only[temp,]
  dmd_data <- rbind(dmd_data,temp2)
  print(dim(dmd_data))
  }

dmd_data$unique <- !(duplicated(dmd_data) | duplicated(dmd_data, fromLast = TRUE))
table(dmd_data$unique)

dmd_data %>%
select(-unique) %>%
distinct() -> dmd_nodup

## combine files together
# make new read2 variable that removes the .00 addition from all codes
# use ifelse statement to find all rows that end with 00 using grepl ($ indicates that 00 is at the end,
# this returns logical vector instead of matched items)
# then use substr to remove last characters from those with 00, for others use normal read2 code

bnf_nodup %>%
bind_rows(read_2_nodup,dmd_nodup) %>%
distinct() %>%
mutate(read_2_no00 = ifelse(grepl("00$",read_2),
  substr(read_2,1,nchar(read_2)-2),
  read_2)) -> bnf_read2_dmd

bnf_read2_dmd$unique <- !(duplicated(bnf_read2_dmd) | duplicated(bnf_read2_dmd, fromLast = TRUE))
table(bnf_read2_dmd$unique)

## map chemical name and functional class to codes:
## annotated files with all unique codes:
bnf_unique_drugnames_mapped <- read.table("bnf_drug_names") 
# this is a text file with the following columns: bnf_code, drug_name, chem_name, functional_class; bnf codes and drug name are usually available in the prescription record data, while chemical name and functional class have to be added; see https://doi.org/10.1101/2020.08.24.20178715 for a copy of this table and/or contact chiara.fabbri@yahoo.it to ask for a copy of this table used in: doi: https://doi.org/10.1101/2020.08.24.20178715  


dmd_unique_drugnames_mapped <- read.table("dmd_drug_names",
  col_types = cols(dmd_code = col_character())) 
# this is a text file with the following columns: dmd_code, drug_name, chem_name, functional_class; dmd codes and drug name are usually available in the prescription record data, while chemical name and functional class have to be added; see https://doi.org/10.1101/2020.08.24.20178715 and/or contact chiara.fabbri@yahoo.it for a copy of this table used in: doi: https://doi.org/10.1101/2020.08.24.20178715

read2_unique_drugnames_mapped <- read.table("read2_drug_names") 
# this is a text file with the following columns: read_2_no00, drug_name, chem_name, functional_class; read_2_no00 codes were derived as showed above, for the other variables see comments above  

names(read2_unique_drugnames_mapped)[1]<-"code"
read2_unique_drugnames_mapped$code_system <- "read_2"
names(bnf_unique_drugnames_mapped)[1]<-"code"
bnf_unique_drugnames_mapped$code_system <- "bnf"
names(dmd_unique_drugnames_mapped)[1]<-"code"
dmd_unique_drugnames_mapped$code_system <- "dmd"

## map antidep_codes_mapping to full medication file to get chemical name and functional class for all rows
## for read2 not all drugs have a drug name so map on code only, 
## therefore I make a file which does not have duplicated read2 codes

read2_unique_drugnames_mapped %>%
distinct(code, .keep_all=T) %>%
select(-drug_name,-n) -> read2_distinct_codes_mapped

# map read2 codes by code only
# map bnf and dmd code by both code and drug name to get unique combinations as some codes have more than 1 medication name
bnf_read2_dmd %>%
left_join(read2_distinct_codes_mapped, by=c("read_2_no00"="code")) %>%
left_join(bnf_unique_drugnames_mapped, by=c("bnf_code"="code","drug_name"="drug_name")) %>%
left_join(dmd_unique_drugnames_mapped, by=c("dmd_code"="code","drug_name"="drug_name")) %>%
mutate(
  chem_name = coalesce(chem_name.x, chem_name.y,chem_name),
  func_class = coalesce(functional_class.x, functional_class.y,functional_class),
  code_system = coalesce(code_system.x,code_system.y,code_system)) %>%
select(eid,data_provider,issue_date,read_2_no00,bnf_code,dmd_code,drug_name,quantity,chem_name,func_class,code_system) %>%
filter(!is.na(issue_date)) -> bnf_read2_dmd_mapped

# write output
fwrite(bnf_read2_dmd_mapped, "ADs_extracted.txt",col.names=T,row.names=F,quote=F,sep="\t")


###########################################
### MERGE DIAGNOSIS AND MEDICATION DATA ###
###########################################

# adding prescription data to diagnosis data
psych_codes %>%
  bind_rows(bnf_read2_dmd_mapped) -> dep_med

# create new phenotypes:
# xx_readcode: specific read code for depression/psychosis/bipolar/substance abuse disorders
# xx_ncode: count the total number of codes per individual per disorder
# xx_code_distinct: count number of distinct code per individual per disorder
# dep_diagnosis_noqc: depression diagnosis (0|1) without removing individuals with psychosis/bipolar disorder/substance abuse disorders
# dep_diagnosis_qc: depression diagnosis (0|1) excluding individuals with psychosis/bipolar disorder/substance abuse disorders


dep_med %>%
  group_by(eid) %>%
  mutate(depression_readcode=ifelse(type=="depression", read_3, NA),
         psychosis_readcode=ifelse(type=="psychosis", read_3, NA),
         bipolar_readcode=ifelse(type=="bipolar_disorder", read_3, NA),
         sub_abuse_readcode=ifelse(type=="substance_abuse", read_3, NA),
         depression_ncode=sum(!is.na(depression_readcode)),
         psychosis_ncode=sum(!is.na(psychosis_readcode)),
         bipolar_ncode=sum(!is.na(bipolar_readcode)),
         sub_abuse_ncode=sum(!is.na(sub_abuse_readcode)),
         chem_name_ncode = sum(!is.na(chem_name)),
         depression_code_distinct=n_distinct(depression_readcode,na.rm=TRUE),
         psychosis_code_distinct=n_distinct(psychosis_readcode,na.rm=TRUE),
         bipolar_code_distinct=n_distinct(bipolar_readcode,na.rm=TRUE),
         sub_abuse_code_distinct=n_distinct(sub_abuse_readcode,na.rm=TRUE),
         chem_name_code_distinct = n_distinct(chem_name, na.rm=TRUE),
         dep_vs_med= ifelse(!is.na(type), "depression",
                            ifelse(!is.na(code_system), "medication", NA)),
         dep_diagnosis_noqc = ifelse((depression_ncode >=2 | depression_code_distinct >=2),1,0),
         dep_diagnosis_qc = ifelse(((depression_ncode >=2 | depression_code_distinct >=2) & 
                                      (sub_abuse_code_distinct==0 & bipolar_code_distinct==0 & psychosis_code_distinct==0)),1,0)) %>% 
  ungroup() -> dep_med_2

# filter down to individuals who have depression read code
dep_med_2 %>%
  group_by(eid) %>%
  filter(any(!is.na(depression_readcode))) %>%
  ungroup() -> gp_depression

gp_depression$issue_date <- as.Date(gp_depression$issue_date, "%d/%m/%Y")

# dep_vs_med : column for distinguishing lines reporting a diagnosis (e.g. depression, psychosis) vs. antidepressant prescriptions

# write the output:
fwrite(gp_depression, "dep_ADs_data.txt", col.names=T,row.names=F,quote=F,sep="\t")


########################################
####### CREATE TRD PHENOTYPE ###########
########################################

dep <- read.delim("dep_ADs_data.txt", h=T,stringsAsFactors=F) ## this file includes diagnosis of depression and antidepressant medications prescriptions, see output created above

# selections of individuals with depression diagnosis (dep_diagnosis_qc) and rows containing antidepressant prescription data:
dep<-dep[dep$dep_diagnosis_qc == 1,]
dep<-dep[!dep$chem_name=='',] ## see also column dep_vs_med to select rows with antidepressant prescription data

# delete non-needed columns:
dep <- dep[,-c(2:10,20:36)] # note that if you modified output files in a different way you probably have different columns and this code may need to be adapted

# in one case prescription date is 1902-02-02, changed to missing:
dep$issue_date[dep$issue_date == '1902-02-02' ] <- NA
dep <- dep[!is.na(dep$issue_date),]

# remove duplicated rows:
dep <- distinct(dep)

# thryptophan was not considered as an antidepressant prescription valid for the definition of TRD
dep <- dep[!dep$func_class == 'amino_acid_supplement',]

# recode mirtazapine class to NaSSA:
dep$func_class <- ifelse(dep$chem_nome == 'mirtazapine', 'NaSSA', dep$func_class)

# create drug prescription episodes (sometimes the same antidepressant is prescribed in periods of different years)
# 6 months between prescriptions of the same antidepressant to distinguish different prescription episodes and add a unique number to distinguish different prescription episodes of the same antidepressant in the same individual

dep %>%
  group_by(eid, chem_name) %>%
  arrange(issue_date) %>%
  mutate(prev_drug_date = dplyr::lag(issue_date, n = 1, default = NA)) %>%
  mutate(diff_weeks_drug = as.numeric(difftime(issue_date, prev_drug_date, units = "weeks"))) %>%
  mutate(prescription_episode = ifelse (diff_weeks_drug > 26, seq_along(diff_weeks_drug), 1)) %>%
  mutate(prescription_episode = ifelse (is.na(prescription_episode), 999, prescription_episode)) %>%
  mutate(prescription_episode = ifelse ( prescription_episode == 1,
  NA, prescription_episode)) %>%
  fill(prescription_episode) %>%
  mutate(prescription_episode = ifelse (prescription_episode == 999, 1, prescription_episode)) %>%
  ungroup() -> dep2

# add variables to check proportion of adequate prescription intervals (<= 14 weeks between two consecutive prescriptions)
dep2 %>%
  group_by(eid,chem_name,prescription_episode) %>%
  arrange(issue_date) %>%
  mutate(duration_prescr_ep = as.numeric(difftime(lead(issue_date,1), issue_date, units = "weeks"))) %>%
  mutate(prescription_episode_drug = sum(as.numeric(duration_prescr_ep),na.rm=T)) %>%
  mutate(adequate_prescr_period = ifelse((diff_weeks_drug <= 14 & !is.na(diff_weeks_drug)), 1,
  ifelse(is.na(diff_weeks_drug), NA, 0))) %>% # if there are > 14 weeks between two consecutive prescriptions the adherence to the drug could be low
  mutate(n_adequate_prescr = sum(adequate_prescr_period,na.rm=T)) %>%
  mutate(n_prescr = length(!is.na(adequate_prescr_period))) %>%
  mutate(adequate_prescr_prop = ifelse((!is.na(adequate_prescr_period) & !is.na(n_adequate_prescr) & n_adequate_prescr >0), n_adequate_prescr / n_prescr, NA)) %>%
  mutate(last_drug_prescr = ifelse(issue_date == max(issue_date), 'yes', 'no')) %>% # this says when each prescription episode for a drug ends
  ungroup() -> dep2

dep2 %>%
  group_by(eid) %>%
  mutate(adequate_prescr_mean = mean(adequate_prescr_prop, na.rm=T)) %>%
  mutate(adequate_prescr_sd = sd(adequate_prescr_prop, na.rm=T)) %>%
  ungroup() -> dep2

# define TRD as at least two antidepressant switches, each drug prescribed for >= 6 weeks and no longer than 14 weeks between the two drugs
dep2 %>%
  group_by(eid) %>%
  arrange(issue_date) %>%
  mutate(drug_distinct = n_distinct(chem_name,na.rm=T)) %>%
  mutate(drug_switch = ifelse((chem_name != lead(chem_name,1) & last_drug_prescr == 'yes'), 'yes', 'no')) %>%
  mutate(prescription_episode_drug = ifelse((prescription_episode_drug == 0 & as.numeric(difftime(lead(issue_date,1), issue_date, units = "weeks")) <=14), as.numeric(difftime(lead(issue_date,1), issue_date, units = "weeks")), prescription_episode_drug)) %>%
  mutate(between_drugs_weeks = ifelse(drug_switch == 'yes', abs(as.numeric(difftime(issue_date, lead(issue_date,1), units = "weeks"))), 999)) %>%
  mutate(treatment_resistance_switch = ifelse((prescription_episode_drug >= 6 & between_drugs_weeks <= 14 & between_drugs_weeks > 0 &  drug_switch == 'yes' & drug_distinct > 2), 1, 0)) %>%
  mutate(treatment_resistance_drug = ifelse(sum(as.numeric(treatment_resistance_switch),na.rm=T) >= 2, 1, 0)) %>%
  ungroup() -> dep2


# check if some individuals switched starting from the same drug more than once (if so, this should be counted as one switch only):
dep3 <- dep2[dep2$treatment_resistance_switch == 1,]
dep3 <- dep3[!is.na(dep3$chem_name),]
ids <- subset(dep3$eid, !duplicated(dep3$eid))

library(foreach) (you may need to install this package)
count <- foreach(i=1:length(ids),.combine='rbind') %dopar% data.frame(eid=ids[i], n_rep=length(subset(dep3[dep3$eid == ids[i],]$chem_name, duplicated(dep3[dep3$eid == ids[i],]$chem_name))), tot=nrow(dep3[dep3$eid == ids[i],]))
count$diff <- ifelse(count$n_rep > 0, count$tot - count$n_rep, NA)
sub <- count[count$diff <= 2 & !is.na(count$diff),]
# if count$diff = 2, it means two valid (non-repeated drugs) switches were done
# if count$diff = 1, only one non-repeated drug switch was done, these should not be considered as TRD
del_drug <- sub[sub$diff == 1,] # to check number of individuals
colnames(del_drug) <- c('eid', 'n_rep_switch_drug', 'tot_drug_swtiches', 'diff_drug_sw')
dep2 <- merge(dep2, del_drug, by='eid', all.x=T)

dep2$treatment_resistance_drug <- ifelse((dep2$diff_drug_sw==1 & !is.na(dep2$diff_drug_sw)), 0, dep2$treatment_resistance_drug)

# write the output:
fwrite(dep2, "TRD_phenotype.txt", col.names=T,row.names=F,quote=F,sep="\t")

