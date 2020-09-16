### R code to extract all psychiatric disorders except mood, psychotic and substance use disorders (extracted using the script extract_diagn_ADs_TRD_pheno.R), seizure disorders and prescription primary care data for antipsychotics, mood stabilizers (only lithium, lamotrigine, valproate and pregabalin are considered), anxiolytics/hypnotics 
### using UK Biobank resource (https://www.ukbiobank.ac.uk) or similar electronic health records of primary care; consider to exclude those with a diagnosis of seizure disorders when extracting some mood stabilisers and benzodiazepines in individuals with mood disorders ###

### information on primary care data can be found in Resource 591 (http://biobank.ndph.ox.ac.uk/showcase/refer.cgi?id=591)
### Information on codes used in the primary care data can be found in Resource 592 (http://biobank.ndph.ox.ac.uk/showcase/refer.cgi?id=592)

### the code may need to be adapted (e.g. file names, file formats)

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


# open clinical data
clinical <- fread("gp_diagnoses",header=T, data.table=F, na.strings="") # table with diagnostic records that can be downloaded from UK Biobank website (https://www.ukbiobank.ac.uk) by investigators having access to the data (field ID 42040; check the name of the file you downloaded as it is probably different)

# open Read2 and Read3 code files (this includes all codes for anxiety and other disorders except mood, psychotic and substance related disorders)
readv2 <- fread("read_2_codes_other",header=T,data.table=F) 
# this is a text table with the following columns: readv2_code, term_v2, readv3_code, term_v3, description, diagnostic_group for the diagnoses of interest; these information can be extracted using documentation provided at: http://biobank.ndph.ox.ac.uk/showcase/showcase/auxdata/primarycare_codings.zip 
# alternatively, see https://doi.org/10.1101/2020.08.24.20178715 for a copy of the table and/or contact chiara.fabbri@yahoo.it to ask for a copy of this table used in: doi: https://doi.org/10.1101/2020.08.24.20178715

readv3 <- fread("read_3_codes_other",header=T,data.table=F) 
# this is a text table with the following columns: readv3_code, term_v3, description, diagnostic_group for the diagnoses of interest; these information can be extracted using documentation provided at: http://biobank.ndph.ox.ac.uk/showcase/showcase/auxdata/primarycare_codings.zip 
# alternatively, see https://doi.org/10.1101/2020.08.24.20178715 for a copy of the table and/or contact chiara.fabbri@yahoo.it to ask for a copy of this table used in: doi: https://doi.org/10.1101/2020.08.24.20178715  


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

# mark all rows where the date is 01/01 (this is likely some issue with the system defaulting to 01/01, but the year may be correct)
psych_codes$diagn_month_day <- format(as.Date(psych_codes$event_dt), "%m-%d")
table(psych_codes$diagn_month_day=='01-01')
psych_codes$diagn_year <- format(as.Date(psych_codes$event_dt),"%Y")

# remove duplicated lines:
psych_codes <- distinct(psych_codes)

# recode some categories (this was done just because of some inconsistencies in the table used for annotation of diagnoses, so you may not need this part):
psych_codes$diagnostic_group[psych_codes$diagnostic_group=='habit_impulse_disorder'] <- 'habit_impulse_disorders'
psych_codes$diagnostic_group[psych_codes$diagnostic_group=='manic_symptoms'] <- 'hypo_manic_symptoms'
psych_codes$description <- ifelse(psych_codes$term_v3 == 'Panic attack', 'panic_disorder', psych_codes$description) ## panic attack was wrongly classified as phobia
psych_codes$diagnostic_group <- ifelse(psych_codes$description == 'mixed_anxiety_depressive_disorder', 'depression_with_anxiety', psych_codes$diagnostic_group)


# count the number of diagnostic codes per individual per the main diagnostic categories (these may vary depending from how you annotated the diagnostic categories):
psych_codes %>%
  group_by(eid) %>%
  mutate(anx_ncode=sum(diagnostic_group=='anxiety_disorders')) %>% 
  mutate(anx_symp_ncode=sum(diagnostic_group=='anxiety_no_diagnosis'|diagnostic_group=='anxiety_symptoms')) %>% # this category included anxiety symptoms that didn't satisfy the criteria for a disorder
  mutate(dep_anx_ncode=sum(diagnostic_group=='depression_with_anxiety')) %>%
  mutate(eat_ncode=sum(diagnostic_group=='eating_disorders')) %>%
  mutate(OCD_ncode=sum(diagnostic_group=='obsessive_compulsive_and_related_disorders')) %>%
  mutate(pers_dis_ncode=sum(diagnostic_group=='personality_disorders')) %>%
  mutate(somatoform_ncode=sum(diagnostic_group=='somatoform_disorders'|diagnostic_group=='somatoform_anxiety_disorders')) %>%
  mutate(stress_dis_ncode=sum(diagnostic_group=='stress_related_disorders')) %>%
  mutate(sleep_dis_ncode=sum(diagnostic_group=='sleep_disorders')) %>%
  mutate(suic_self_harm_ncode=sum(diagnostic_group=='attempted_suicide_self_harm'|diagnostic_group=='self_harm_suicide')) %>%
  ungroup() -> psych_codes2

# write the output:
fwrite(psych_codes2, "other_psych_diagnoses.txt",col.names=T,row.names=F,quote=F,sep="\t")



################################################
######### PRESCRIPTION DATA EXTRACTION #########
################################################

## This part of the code extracts psychotropic medication codes from primary care data (antipsychotics, mood stabilizers, anxiolytics/hypnotics)
## used medication codes: Read 2, BNF, dm+d

# open full medication data
meds <- 
  fread("gp_prescriptions", # this is the table with prescription records (field ID 42039) that can be downloaded from UK Biobank website (https://www.ukbiobank.ac.uk) by investigators having access to the data, check file name as it is probably different to the one in this script
  header=T, data.table=F, na.strings="", 
  colClasses = c(
      "read_2" = "character",
      "bnf_code" = "character",
      "dmd_code" = "character"))


#####################
###### BNF ##########
#####################

## extract all BNF codes (this will also include part of read2 codes, as some individuals will have both codes)
## 0408 -> lamotrigine
## 0402 -> lithium and valproic acid, antipsychotics
## 0401 -> anxiolytics/hypnotics

bnf_codes <- c(
  "0408","04.08","0402", "04.02", "0401", "04.01") ## these capture also other drugs that then are deleted

bnf_data <- data.frame()
for (i in bnf_codes){
  print(paste("started BNF code",i, "at", Sys.time()))
  temp <- grep(paste0("^",i,"."),meds$bnf_code,ignore.case=T)
  temp2 <- meds[temp,]
  bnf_data <- rbind(bnf_data,temp2)
  print(dim(bnf_data))
  rm(temp2)
  }

# identify number of duplicates
bnf_data$unique <- !(duplicated(bnf_data) | duplicated(bnf_data, fromLast = TRUE))
table(bnf_data$unique) 

# remove duplicates that have been added due to the way we extract the data
bnf_data %>%
select(-unique) %>%
distinct() -> bnf_nodup

#######################
###### Read2 ##########
#######################

# load clinical read2 codes for anxiolytics/hypnotics, antipsychotics and mood stabilizers (antipsychotic formulations used for other indications were removed)
read_2_codes <- read.table("read2_med_codes",header=T) 
# this corresponds to a text file with the following column: read2_code, with the read2 codes for the medications of interest except the ones specified below; this can be can be extracted using documentation provided at: http://biobank.ndph.ox.ac.uk/showcase/showcase/auxdata/primarycare_codings.zip
# alternatively, see https://doi.org/10.1101/2020.08.24.20178715 for a copy of the table and/or contact chiara.fabbri@yahoo.it to ask for a copy of this table used in: doi: https://doi.org/10.1101/2020.08.24.20178715

read_2_codes2 <- read.table("read2_valpr_codes",header=T) # same as previous, but only for valproate

read_2_codes3 <- read.delim("read2_clon_preg_codes", h=T) # same as previous, but only for clonazepam and pregabalin (these medications read 2 codes are split just because they were extracted at different times of the analyses)
read_2_codes3<-data.frame(read2_code=read_2_codes3[,1]) # this table had a second column including drug name that was not needed at that point and it was removed

read_2_codes <- rbind(read_2_codes, read_2_codes2, read_2_codes3) # put the read 2 codes of these medications together
read_2_codes <- read_2_codes[!duplicated(read_2_codes[,1]),] # remove duplicated codes
read_2_codes <- data.frame(read2_code=read_2_codes)

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
###### dm+d ###########
#######################


# medication names: 

dmd_names_antipsychotics <-c("abilify", "amisulpride", "anquil", "aripiprazole", "arpoya", "benperidol", "cariprazine", "reagila", "chlorpromazine", "largactil", "chlorprothixene", "truxal", "clotiapine", "clozapine", "clozaril", "denzapine", "zaponex", "droperidol", "xomolix", "depixol", "fluanxol", "flupenthixol", "psytixol", "fluphenazine", "modecate", "moditen", "dozic", "haldol", "haloperidol", "serenace", "levomepromazine", "nozinan", "latuda", "lurasidone", "melperon", "melperone", "olanzapine", "zalasta", "zypadhera", "zyprexa", "invega", "paliperidone", "trevicta", "xeplion", "neulactil", "pericyazine", "fentazin", "perphenazine", "orap", "pimozide", "prazine", "promazine", "atrolak", "biquelle", "brancico", "mintreleq", "psyquet", "quetiapine", "seotiapim", "seroquel", "sondate", "tenprolide", "zaluron", "risperdal", "risperidone", "serdolect", "sertindole", "dolmatil", "sulpiride", "sulpitil", "melleril", "rideril", "thioridazine", "stelazine", "terrazine", "trifluoperazine", "zeldox", "ziprasidone", "ciatyl-z", "clopixol", "zuclopenthixol", "asenapine", "sycrest")  

dmd_names_mood_stabilizers <- c ("lithium", "priadel", "camcolit", "liskonum", "li-Liquid", "litarex", "efalith", "valproic", "depakote", "convulex", "lamotrigine", "lamictal")  

dmd_names_anxiolytics_hypno <- c ("alprazolam", "xanax", "buspirone", "buspar", "chlordiazepoxide", "librium", "tropium", "libraxin", "chlormezanone", "trancopal", "bromazepam", "lexotan", "diazepam", "tensium", "diazemuls", "valium", "solis", "stesolid", "atensine", "alupram", "rimapam", "dialar", "ketazolam", "anxon", "lorazepam", "ativan", "medazepam", "nobrium", "meprobamate", "equanil", "meprate", "oxazepam", "prazepam", "centrax", "clorazepate", "tranxene", "slenyto", "s.gard", "melatonin", "circadin", "chloral hydrate", "noctec", "welldorm", "cloral betaine", "clomethiazole", "heminevrin", "flunitrazepam", "rohypnol", "flurazepam", "dalmane", "loprazolam", "dormonoct", "lormetazepam", "nitrazepam", "remnos", "mogadon", "nitrados", "somnite", "unisomnia", "potassium bromide", "dibro-be", "temazepam", "normison", "triazolam", "halcion", "zaleplon", "sonata", "triclofos sodium", "zolpidem", "stilnoct", "zopiclone", "zimovane", "zileze", "midazolam", "epistatus", "buccolam", "hypnovel", "dormicum")   

dmd_names_others <- c('pregabalin', 'alzain', 'axalid', 'lecaent', 'lyrica', 'rewisca', 'clonazepam', 'rivotril', 'klonopin','solian', 'benquil', 'chlorazin', 'flupentixol', 'halkid', 'levinan', 'arkolamyl', 'periciazine', 'sulpor', 'lithonate','belvo','syonell', 'librax', 'ozalin')

dmd_names <- c(dmd_names_antipsychotics, dmd_names_mood_stabilizers, dmd_names_anxiolytics_hypno, dmd_names_others)

# subset to only have dmd codes otherwise I'll extract a lot of read 2 people as well (this is purely to have smaller file)
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
table(bnf_read2_dmd$unique) ## no duplicates

# write the output:
fwrite(bnf_read2_dmd, "extracted_med_data_raw.txt",col.names=T,row.names=F,quote=F,sep="\t") 


# remove medications of no interest having non-missing drug_name: 
meds_to_excl <- c('gabapentin', 'carbamazepine', 'perampanel', 'phenobarbital', 'oxcarbazepine', 'phenytoin', 'phenobarbitone', 'tiagabine', 'topiramate', 'primidone', 'zonisamide', 'vigabatrin', 'lacosamide', 'ethosuximide', 'levetiracetam', 'retigabine', 'rufinamide', 'eslicarbazepine', 'keppra', 'mysoline', 'neurontin', 'sabril', 'zonegran', 'vimpat', 'trileptal','fycompa', 'carbagen', 'tegretol','epanutin', 'Sodium Amytal')  

meds_excl <- data.frame()
for (i in meds_to_excl){
  print(paste("started dm+d code",i, "at", Sys.time()))
  temp <- grep(paste0("^",i,"."),bnf_read2_dmd$drug_name,ignore.case=T)
  temp2 <- bnf_read2_dmd[temp,]
  meds_excl <- rbind(meds_excl, temp2)
  }

meds <- bnf_read2_dmd[!bnf_read2_dmd$drug_name %in% meds_excl$drug_name,]

# annotated files to add chemical name and drug class/category:
ma <- read.delim("other_med_drug_names.txt",h=T) 
# this is a text file including the following columns: drug_name, read_2, bnf_code, dmd_code, read_2_no00, chem_name, category and class; medication codes and drug name are usually available in prescription records, other information needs to be added; see https://doi.org/10.1101/2020.08.24.20178715 for a copy of the table and/or contact chiara.fabbri@yahoo.it to ask for a copy of this table used in: doi: https://doi.org/10.1101/2020.08.24.20178715 

mb <- read.delim("other_med_drug_read_2_no00.txt",h=T) 
# this is a text file including the following columns: drug_name, read_2, bnf_code, dmd_code, read_2_no00, chem_name, category and class; see previous line for other info. this second file was created to include medications having no drug name but only read 2 codes available in the prescription records 

# select rows with non-missing drug name:
meds1 <- meds[!is.na(meds$drug_name),]
# or meds1 <- meds[meds$drug_name != '',] # this in case the missing are read as ''

# select rows with missing drug name that will be annotated using read2 codes:
meds2 <- meds[is.na(meds$drug_name),]
# or meds2 <- meds[meds$drug_name == '',] # this in case the missing are read as ''

# exclude some medication of no interest for this study (anticonvulsants):
meds1 <- meds1[!(meds1$drug_name == 'Brivaracetam  Tablets  25 mg' | meds1$drug_name == 'Brivaracetam  Tablets  50 mg' | meds1$drug_name == 'EFALITH cream 8% + 0.05%' | meds1$drug_name == 'Zarontin 250mg capsules (Pfizer Ltd)'),]
meds1_ex <- grep("^topamax", meds1$drug_name, ignore.case=T)
meds1 <- meds1[-meds1_ex,]

ma <- ma[ma$drug_name %in% meds1$drug_name,]
ma <- ma[,-c(2:5)] # check before removing columns
meds1 %>%
  left_join(ma, by=c("drug_name"="drug_name")) -> meds1
meds1 <- meds1[,-10] # check before removing columns

mb <- mb[mb$read_2_no00 %in% meds2$read_2_no00,]
mb <- mb[,-c(2:4)] # check before removing columns
meds2 <- meds2[,-c(7,10)] # check before removing columns
meds2 %>%
  left_join(mb, by=c("read_2_no00"="read_2_no00")) -> meds2
meds2 <- data.frame(eid=meds2$eid, data_provider=meds2$data_provider, issue_date=meds2$issue_date, read_2=meds2$read_2, bnf_code=meds2$bnf_code, dmd_code=meds2$dmd_code, drug_name=meds2$drug_name, quantity=meds2$quantity, read_2_no00=meds2$read_2_no00, chem_name=meds2$chem_name, category=meds2$category, class=meds2$class) # to rearrange column order

meds <- rbind(meds1, meds2) # bind the two data.frames

meds$class <- ifelse(meds$chem_name == 'flunitrazepam', 'benzodiazepine', meds$class) # add missing class for 2 flunitrazepam entries

# write the output:
fwrite(meds, "extracted_med_data_annotated.txt",col.names=T,row.names=F,quote=F,sep="\t")


#####################################################
########## Antidepressant combinations ##############
#####################################################

dep2 <- read.delim("TRD_phenotype.txt",h=T) # see script extract_diagn_ADs_TRD_pheno.R of this repository for the creation of this file that includes antidepressant prescriptions

dep2 <- dep2[!is.na(dep2$issue_date),] # remove lines with missing date

library(timetools)
# This first step only identifies individuals who received antidepressant combinations to reduce the number of data included in the following step, but it is not necessary and code could be adapted to run step 2 directly 

dep2$issue_date<-as.Date(dep2$issue_date, "%Y-%m-%d")

dep2 %>%
  group_by(eid,chem_name,prescription_episode) %>%
  arrange(issue_date) %>%
  mutate(AD_start = min(issue_date, na.rm=T)) %>%
  mutate(AD_end = max(issue_date, na.rm=T)) %>%
  ungroup() -> dep2

time_data <- data.frame(eid=dep2$eid, chem_name=dep2$chem_name, func_class=dep2$func_class, AD_start=as.POSIXct(dep2$AD_start, "%Y-%m-%d"), AD_end=as.POSIXct(dep2$AD_end, "%Y-%m-%d"))
time_data <- distinct(time_data)
time_data <- time_data[order(time_data$eid, time_data$chem_name),]

ids<-unique(time_data$eid)
l<-list()
for(i in 1:length(unique(time_data$eid))) l[[i]] <- time_data[time_data$eid ==ids[i],]

intervals<-list()
for(i in 1:length(unique(time_data$eid))) intervals[[i]] <- TimeIntervalDataFrame(l[[i]]$AD_start, l[[i]]$AD_end)

overlap<-vector()
for(i in 1:length(unique(time_data$eid))) overlap[i] <- overlapping(intervals[[i]])

data<-data.frame(eid=unique(time_data$eid), overlap=overlap)

AD_comb <- time_data[time_data$eid %in% data[data[,2]=='TRUE',]$eid,] 

# This second step identifies which antidepressants were prescribed in combination and for how long; you could adapt the script to run directly step 2 without the first step described above
detach("package:timetools", unload=TRUE)
AD_comb$AD_start <- as.character(as.Date(AD_comb$AD_start, "%Y-%m-%d"))
AD_comb$AD_end <- as.character(as.Date(AD_comb$AD_end, "%Y-%m-%d"))
AD_comb$AD_interval <- interval(AD_comb$AD_start, AD_comb$AD_end,tzone="UTC")

AD_comb$AD_interval_2 <- AD_comb$AD_interval

ids<-unique(AD_comb$eid)
l<-list()
for(i in 1:length(unique(AD_comb$eid))) l[[i]] <- crossing(i1=AD_comb[AD_comb$eid == ids[i],]$AD_interval, i2=AD_comb[AD_comb$eid == ids[i],]$AD_interval_2) # identifies combinations of prescription intervals
for(i in 1:length(unique(AD_comb$eid))) l[[i]]$intersect <- day(as.period(intersect(l[[i]]$i1, l[[i]]$i2), "days")) # number of overlapping days
for(i in 1:length(unique(AD_comb$eid))) l[[i]] <- l[[i]] [!l[[i]][,1] == l[[i]][,2],] # excludes lines where the prescription interval is identical (same drug)
for(i in 1:length(unique(AD_comb$eid))) l[[i]] <- l[[i]] [!is.na(l[[i]]$intersect),]
for(i in 1:length(unique(AD_comb$eid))) colnames(l[[i]])[1] <- 'AD_interval'
for(i in 1:length(unique(AD_comb$eid))) l[[i]] <- merge(AD_comb[AD_comb$eid == ids[i],], l[[i]], by='AD_interval')

AD_comb_2<-do.call(rbind.data.frame, l)

AD_comb_2$TRD_AD_comb <- ifelse(AD_comb_2$intersect > 30, 1, 0) ## this selects individuals taking both antidepressants for longer than 30 days, in order to exclude co-prescription for cross-tapering between drugs (during switch from a drug to another)

AD_comb_2 <- AD_comb_2[,-7]
colnames(AD_comb_2)[7] <- "AD_2_interval"
AD_comb_2$AD_interval <- as.character(AD_comb_2$AD_interval) # conversion of the interval to character to avoid writing errors
AD_comb_2$AD_2_interval <- as.character(AD_comb_2$AD_2_interval)

fwrite(AD_comb_2, "AD_comb.txt",col.names=T,row.names=F,quote=F,sep="\t")


#####################################################
#### Antidepressant - antipsychotic combinations ####
#####################################################

meds <- read.delim("extracted_med_data_annotated.txt", h=T, stringsAsFactors=F) # see output created previously
ap <- meds[meds$category == 'typical_AP' | meds$category == 'atypical_AP',] # select lines including prescription records of antipsychotics
ap <- ap[ap$eid %in% dep2$eid,] 

ap$issue_date <- dmy(ap$issue_date)
ap$issue_date <- as.Date(ap$issue_date, "%Y-%m-%d")

ap <- ap[,-(5:6)] # before removing columns, check if that's what you'd like to do and consider that if you modified your files in a different way compared to these scripts you'd need to adapt the scripts to your data
colnames(ap)[9] <- 'func_class' # column 9 here corresponds to the classification of antipsychotics in atypical and typical
ap<-ap[,-10] # before removing columns, check if that's what you'd like to do and consider that if you modified your files in a different way compared to these scripts you'd need to adapt the scripts to your data

# create prescription episodes for antipsychotics (see also extract_diagn_ADs_TRD_pheno.R script):
ap %>%
  group_by(eid, chem_name) %>%
  arrange(issue_date) %>%
  mutate(prev_drug_date = dplyr::lag(issue_date, n = 1, default = NA)) %>%
  mutate(diff_weeks_drug = as.numeric(difftime(issue_date, prev_drug_date, units = "weeks"))) %>%
  mutate(prescription_episode = ifelse (diff_weeks_drug > 24, seq_along(diff_weeks_drug), 1)) %>%
  mutate(prescription_episode = ifelse (is.na(prescription_episode), 999, prescription_episode)) %>%
  mutate(prescription_episode = ifelse ( prescription_episode == 1,
  NA, prescription_episode)) %>%
  fill(prescription_episode) %>%
  mutate(prescription_episode = ifelse (prescription_episode == 999, 1, prescription_episode)) %>%
  ungroup() -> ap

# add start and end date for each prescription episode: 
ap %>%
  group_by(eid,chem_name,prescription_episode) %>%
  arrange(issue_date) %>%
  mutate(AP_start = min(issue_date, na.rm=T)) %>%
  mutate(AP_end = max(issue_date, na.rm=T)) %>%
  ungroup() -> ap

# remove suppositories (formulation used for other disorders):
sup <- grep('supp', ap$drug_name, ignore.case = T)
ap <- ap[-sup,]

dep2 <- read.delim("TRD_phenotype.txt",h=T,stringsAsFactors=F) # see extract_diagn_ADs_TRD_pheno.R of this repository for the creation of this file, it includes antidepressant prescriptions
dep2$prev_drug_date <- as.Date(dep2$prev_drug_date, "%Y-%m-%d")
dep2$issue_date <-as.Date(dep2$issue_date, "%Y-%m-%d")
dep2<-dep2[!is.na(dep2$issue_date),]

dep2 %>%
  group_by(eid,chem_name,prescription_episode) %>%
  arrange(issue_date) %>%
  mutate(AD_start = min(issue_date, na.rm=T)) %>%
  mutate(AD_end = max(issue_date, na.rm=T)) %>%
  ungroup() -> dep2

# extract only individuals taking antipsychotics:
dep2 <- dep2[dep2$eid %in% ap$eid,]

dep2 %>% bind_rows(ap) -> dep3

time_data <- data.frame(eid=dep3$eid, chem_name=dep3$chem_name, func_class=dep3$func_class, AD_start=as.character(as.Date(dep3$AD_start, "%Y-%m-%d")), AD_end=as.character(as.Date(dep3$AD_end, "%Y-%m-%d")),AP_start=as.character(as.Date(dep3$AP_start, "%Y-%m-%d")),AP_end=as.character(as.Date(dep3$AP_end, "%Y-%m-%d")))
time_data <- distinct(time_data)
time_data <- time_data[order(time_data$eid, time_data$chem_name),]

time_data$AD_interval <- interval(time_data$AD_start, time_data$AD_end,tzone="UTC")
time_data$AP_interval <- interval(time_data$AP_start, time_data$AP_end,tzone="UTC")


ids<-unique(time_data$eid)
l<-list()
for(i in 1:length(unique(time_data$eid))) l[[i]] <- crossing(i1=time_data[time_data$eid == ids[i],]$AD_interval, i2=time_data[time_data$eid == ids[i],]$AP_interval)
for(i in 1:length(unique(time_data$eid))) l[[i]]$intersect <- day(as.period(intersect(l[[i]]$i1, l[[i]]$i2), "days")) # number of overlapping days
for(i in 1:length(unique(time_data$eid))) l[[i]] <- data.frame(l[[i]])
for(i in 1:length(unique(time_data$eid))) l[[i]] <- l[[i]] [!is.na(l[[i]]$i1) & !is.na(l[[i]]$i2),] # excludes lines where there are missing time intervals
for(i in 1:length(unique(time_data$eid))) l[[i]] <- l[[i]] [!is.na(l[[i]]$intersect),]
for(i in 1:length(unique(time_data$eid))) colnames(l[[i]])[1] <- 'AD_interval'
for(i in 1:length(unique(time_data$eid))) l[[i]] <- merge(time_data[time_data$eid == ids[i],], l[[i]], by='AD_interval')

comb<-do.call(rbind.data.frame, l)

comb$TRD_AP_comb <- ifelse(comb$intersect > 30, 1, 0) # this selects individuals taking both medications for > 30 days

# extract from time_data the corresponding antipsychotic prescriptions:
app <- time_data[time_data$eid %in% comb[comb$TRD_AP_comb==1,]$eid & !is.na(time_data$AP_start),]
app<-app[,-c(4,5,8)] # before removing columns check if that's what you'd like to do based on your data

comb<-comb[,-c(7:10)] # check before removing columns

# combine them with antidepressant prescriptions:
comb %>%
  left_join(app, by=c("eid"="eid")) -> comb2

comb2 <- comb2[comb2$eid %in% comb2[comb2$TRD_AP_comb == 1,]$eid,]
comb2 <- data.frame(comb2); comb2 <- comb2[order(comb2$eid),]

# check and write combinations:
comb2 %>%
 group_by(eid) %>%
 mutate(overlap = ifelse(int_overlaps(AD_interval,AP_interval) == 'TRUE', 'yes', 'no')) %>% 
 mutate(AD_AP = ifelse(overlap == 'yes', paste(chem_name.x, chem_name.y, sep='_'), NA)) %>%
 mutate(AD_AP_class =  ifelse(overlap == 'yes', paste(func_class.x, func_class.y, sep='_'), NA)) %>%
 ungroup() -> comb2

comb2$AD_interval <- as.character(comb2$AD_interval) # otherwise time intervals do not get written properly in the output table
comb2$AP_interval <- as.character(comb2$AP_interval)

# write output:
fwrite(comb2, "AD_AP_comb.txt",col.names=T,row.names=F,quote=F,sep="\t")


#####################################################
### Antidepressant - mood stabilizers combinations ##
#####################################################

#### considering only antidepressant augmentation with lithium, lamotrigine, valproate or pregabalin #####

dep2 <- read.delim("TRD_phenotype.txt",h=T,stringsAsFactors=F) # see extract_diagn_ADs_TRD_pheno.R of this project for the creation of this file, it includes antidepressant prescriptions

meds <- read.delim("extracted_med_data_annotated.txt", h=T, stringsAsFactors=F) # see output created previously
ms <- meds[meds$category == 'mood_stabilizer',] # extract lines including prescription records of mood stabilizers
rm(meds)
ms <- ms[ms$eid %in% dep2$eid,]

ep <- read.delim("epilepsy_diagnoses.txt", h=T,stringsAsFactors=F) # individuals with seizure disorders according to primary care data, data were extracted using the same approach used for psychiatric diagnoses

ms1 <- ms[ms$chem_name == 'lithium',]
ms2 <- ms[ms$chem_name != 'lithium',]
ms2 <- ms2[!ms2$eid %in% ep$eid,] # exclusion of individuals with epilepsy taking lamotrigine, valproate or pregabalin since epilepsy was likely to be the disorder they were prescribed these drugs for
ms <- rbind(ms1, ms2)

ms$issue_date <- dmy(ms$issue_date)
ms$issue_date <- as.Date(ms$issue_date, "%Y-%m-%d")

ms <- ms[,-c(2,4:6,12)] # before removing columns, check if that's what you'd like to do and consider that if you modified your files in a different way compared to these scripts you'd need to adapt the scripts to your data

ms <- ms[!is.na(ms$issue_date),]

# create prescription episodes for mood stabilizers (see also extract_diagn_ADs_TRD_pheno.R script):
ms %>%
  group_by(eid, chem_name) %>%
  arrange(issue_date) %>%
  mutate(prev_drug_date = dplyr::lag(issue_date, n = 1, default = NA)) %>%
  mutate(diff_weeks_drug = as.numeric(difftime(issue_date, prev_drug_date, units = "weeks"))) %>%
  mutate(prescription_episode = ifelse (diff_weeks_drug > 24, seq_along(diff_weeks_drug), 1)) %>%
  mutate(prescription_episode = ifelse (is.na(prescription_episode), 999, prescription_episode)) %>%
  mutate(prescription_episode = ifelse ( prescription_episode == 1,
  NA, prescription_episode)) %>%
  fill(prescription_episode) %>%
  mutate(prescription_episode = ifelse (prescription_episode == 999, 1, prescription_episode)) %>%
  ungroup() -> ms

# add start and end date for each prescription episode: 
ms %>%
  group_by(eid,chem_name,prescription_episode) %>%
  arrange(issue_date) %>%
  mutate(MS_start = min(issue_date, na.rm=T)) %>%
  mutate(MS_end = max(issue_date, na.rm=T)) %>%
  ungroup() -> ms

dep2$prev_drug_date <- as.Date(dep2$prev_drug_date, "%Y-%m-%d")
dep2$issue_date <-as.Date(dep2$issue_date, "%Y-%m-%d")
dep2<-dep2[!is.na(dep2$issue_date),]

dep2 %>%
  group_by(eid,chem_name,prescription_episode) %>%
  arrange(issue_date) %>%
  mutate(AD_start = min(issue_date, na.rm=T)) %>%
  mutate(AD_end = max(issue_date, na.rm=T)) %>%
  ungroup() -> dep2

# extract only patients taking mood stabilizers:
dep2 <- dep2[dep2$eid %in% ms$eid,]

dep2 %>% bind_rows(ms) -> dep3

time_data <- data.frame(eid=dep3$eid, chem_name=dep3$chem_name, func_class=dep3$func_class, AD_start=as.character(as.Date(dep3$AD_start, "%Y-%m-%d")), AD_end=as.character(as.Date(dep3$AD_end, "%Y-%m-%d")),MS_start=as.character(as.Date(dep3$MS_start, "%Y-%m-%d")),MS_end=as.character(as.Date(dep3$MS_end, "%Y-%m-%d")))
time_data <- distinct(time_data)
time_data <- time_data[order(time_data$eid, time_data$chem_name),]

time_data$AD_interval <- interval(time_data$AD_start, time_data$AD_end,tzone="UTC")
time_data$MS_interval <- interval(time_data$MS_start, time_data$MS_end,tzone="UTC")


ids<-unique(time_data$eid)
l<-list()
for(i in 1:length(unique(time_data$eid))) l[[i]] <- crossing(i1=time_data[time_data$eid == ids[i],]$AD_interval, i2=time_data[time_data$eid == ids[i],]$MS_interval)
for(i in 1:length(unique(time_data$eid))) l[[i]]$intersect <- day(as.period(intersect(l[[i]]$i1, l[[i]]$i2), "days")) # number of overlapping days
for(i in 1:length(unique(time_data$eid))) l[[i]] <- data.frame(l[[i]])
for(i in 1:length(unique(time_data$eid))) l[[i]] <- l[[i]] [!is.na(l[[i]]$i1) & !is.na(l[[i]]$i2),] # excludes lines where there are missing time intervals
for(i in 1:length(unique(time_data$eid))) l[[i]] <- l[[i]] [!is.na(l[[i]]$intersect),]
for(i in 1:length(unique(time_data$eid))) colnames(l[[i]])[1] <- 'AD_interval'
for(i in 1:length(unique(time_data$eid))) l[[i]] <- merge(time_data[time_data$eid == ids[i],], l[[i]], by='AD_interval')

comb<-do.call(rbind.data.frame, l)

comb$TRD_MS_comb <- ifelse(comb$intersect > 30, 1, 0)

# I extract from time_data the corresponding mood stabilizer prescriptions:
msp <- time_data[time_data$eid %in% comb[comb$TRD_MS_comb==1,]$eid & !is.na(time_data$MS_start),]
msp<-msp[,-c(4,5,8)] # check before removing colums

comb<-comb[,-c(7:10)] # before removing columns, check if that's what you'd like to do and consider that if you modified your files in a different way compared to these scripts you'd need to adapt the scripts to your data

# combine with antidepressant prescriptions:
comb %>%
  left_join(msp, by=c("eid"="eid")) -> comb2

comb2 <- comb2[comb2$eid %in% comb2[comb2$TRD_MS_comb == 1,]$eid,]
comb2 <- data.frame(comb2); comb2 <- comb2[order(comb2$eid),]

# check and write combinations:
comb2 %>%
 group_by(eid) %>%
 mutate(overlap = ifelse(int_overlaps(AD_interval,MS_interval) == 'TRUE', 'yes', 'no')) %>%
 mutate(AD_MS = ifelse(overlap == 'yes', paste(chem_name.x, chem_name.y, sep='_'), NA)) %>%
 ungroup() -> comb2

comb2$AD_interval <- as.character(comb2$AD_interval) # otherwise time intervals do not get written properly in the output table
comb2$MS_interval <- as.character(comb2$MS_interval)

# write the output:
fwrite(comb2, "AD_MS_comb.txt",col.names=T,row.names=F,quote=F,sep="\t")

