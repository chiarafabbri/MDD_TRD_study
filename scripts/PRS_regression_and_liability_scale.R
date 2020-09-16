##################################################################################################################
################## R code for regression models of phenotype and PRS of traits of interest #######################
##################################################################################################################

### this script was used for the analyses of the study: https://doi.org/10.1101/2020.08.24.20178715

## R version 3.6.0 was used for the cited study

setwd("~/PRS_directory") ### set directory with PRS files obtained by PRSice

file.names <- dir(pattern ="*.all.score") ### read PRS files created with PRSice, where the first two columns are FID and IID, and the following columns correspond to PRS at 11 p thresholds; PRS at different p thresholds were standardised to calculate effect size for easier interpretation

scores<-list()
for (k in 1:length(file.names)){scores[[k]] <- read.table(file.names[k],h=T)}

pheno<-read.table("phenotype_file", h=T) ### read phenotype file (two columns called FID and third column pheno)
cov<-read.table("covariate_file", h=T) ### read covariate file (FID, first six population principal components, batch and assessment_centre)

for (k in 1:length(file.names)) {scores[[k]]<-merge(scores[[k]],pheno, by='FID') } ### merging of PRS and phenotype files by individual ID (FID)

for (k in 1:length(file.names)) {scores[[k]]<-merge(scores[[k]],cov, by='FID') } ### merging of PRS-phenotype and covariate files by individual ID (FID)

for (k in 1:length(file.names)) {scores[[k]][,3:13] <- apply(scores[[k]][3:13], 2, scale)} ### scaling of the scores 

m1 <-list()

# Regression at 11 p threshold for all PRS (assuming that PRS at different threshold are in columns 3-13 of each element in the list called scores), covariates are the first 6 population principal components (PC), batch (coded as factor) and assessment_centre (factor)

for (i in 1:length(scores))
m1[[i]]<-apply(scores[[i]][3:13],2,FUN <- function(x) glm(pheno ~ x + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + batch + as.factor(assessment.centre), family='binomial',data=scores[[i]]))

results<-list()

for (i in 1:length(scores))
results[[i]] <- data.frame (do.call(rbind,lapply(m1[[i]], function(z)summary(z)$coefficients[2,])))

for (i in 1:length(results)) results[[i]]$OR <- exp(results[[i]]$Estimate) ### adds OR

for (i in 1:length(results)) results[[i]]$ci_l <- exp(results[[i]]$Estimate-(1.96 * results[[i]]$Std..Error)) ### adds lower limit of 95% CI

for (i in 1:length(results)) results[[i]]$up_l <- exp(results[[i]]$Estimate+(1.96 * results[[i]]$Std..Error)) ### adds upper limit of 95% CI

# calculate Nagelkerke R2

library(sizeMat) # you may need to install this package

# Nagelkerke R2 of the whole model
r2 <- list()
for (i in 1:length(results))
r2[[i]]<-data.frame(do.call(rbind,lapply(m1[[i]],function(z)nagelkerkeR2(z))))

# null model (including only covariates)
n <- list()
for (i in 1:length(results)) n[[i]] <- glm(pheno ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + batch + as.factor(assessment.centre), family='binomial',data=scores[[i]]))

# difference between whole model and null model
for (i in 1:length(results)) r2[[i]]$prs.r2 <- r2[[i]][,1] - nagelkerkeR2(n[[i]])

for (i in 1:length(results)) results[[i]]$full.r2 <- r2[[i]][,1]
for (i in 1:length(results)) results[[i]]$prs.r2 <- r2[[i]]$prs.r2

# write the results
setwd("~/results_directory/") # set directory you want to write the results to

for (i in 1:length(results)) 
write.table(results[[i]], paste(file.names[i], "results", sep='_'), col.names=T,row.names=T, sep=' ', quote=F) ### writes the results as separate tables for each PRS in the set directory


############################################################
####### convert Nagelkerke R2 to the liability scale #######
############################################################

### R function to covert Nagelkerke R2 to the liability scale

h2l_R2N <- function(k, r2n, p) {
  # k baseline disease risk
  # r2n Nagelkerke's attributable to genomic profile risk score
  # proportion of sample that are cases
  # calculates proportion of variance explained on the liability scale
  #from ABC at http://www.complextraitgenomics.com/software/
  #Lee SH, Goddard ME, Wray NR, Visscher PM. (2012) A better coefficient of determination for genetic profile analysis. Genet Epidemiol. 2012 Apr;36(3):214-24.
  x <- qnorm(1 - k)
  z <- dnorm(x)
  i <- z / k
  cc <- k * (1 - k) * k * (1 - k) / (z^2 * p * (1 - p))
  theta <- i * ((p - k)/(1 - k)) * (i * ((p - k) / ( 1 - k)) - x)
  e <- 1 - p^(2 * p) * (1 - p)^(2 * (1 - p))
  h2l_R2N <- cc * e * r2n / (1 + cc * e * theta * r2n)
}
