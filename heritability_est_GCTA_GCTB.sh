#!/bin/bash

##### the following code was used to estimate SNP-based heritability using GCTA GREML and GCTB bayes S
##### this was run as bash script, you may need to adapt it depending from your OS and download/install/load the required software 

##############################################################
#######################  GCTA GREML ##########################
##############################################################

## GCTA software was developed by Jian Yang with supports from Peter Visscher, Mike Goddard and Hong Lee (https://cnsgenomics.com/software/gcta/#Overview)
## version 1.93.1beta of GCTA was used
## the software supports the use of multiple threads using the option --threads or --thread-num that can be added to the following code 


# create the genetic relationship matrix using 10 subsets of the genetic data in plink binary format:
for i in $(seq 1 10); do 
~/gcta64 \
--make-grm-part 10 ${i} \
--bfile ~/plink_data \
--keep ~/cases.vs.controls \ ## keep subset of subjects of interest passing quality control
--extract ~/post_qc_snp_list.txt \ ## extract SNPs that passed quality control
--out data_gcta
done

# concatenate the different subsets:
cat data_gcta.part_10*.grm.id > data_GCTA.grm.id
cat data_gcta.part_10*.grm.bin > data_GCTA.grm.bin
cat data_gcta.part_10*.grm.N.bin > data_GCTA.grm.N.bin

# remove the subsets:
rm data_gcta.part_10*

# adjust GRM for incomplete tagging of causal SNPs:
~/gcta64 \
--grm data_GCTA \
--grm-adj 0 \
--make-grm \
--out data_GTCA.adjusted \

# remove the previous data:
rm data_GCTA.grm.*

# Remove related subjects using grm-cutoff 0.05:
~/gcta64 \
--grm data_GTCA.adjusted \
--grm-cutoff 0.05 \
--make-grm \
--out data_GTCA.adjusted.unrel \

# Run GREML-SC:
~/gcta64 \
--grm data_GTCA.adjusted.unrel \
--reml \
--reml-no-constrain \
--pheno ~/pheno \
--covar ~/categorical_cov \ ## file with categorical covariates
--qcovar ~/quantitative_cov \ ## file with quantitative covariates
--out data_GCTA.adjust.unrel.GREML_SC ## output file


##############################################################
##################### GCTB analysis ##########################
##############################################################

## GCTB software was developed by Jian Zeng with supports from Jian Yang, Futao Zhang and Zhili Zheng (https://cnsgenomics.com/software/gctb/#Overview)

## MPI version 2.0 of GCTB was used (https://cnsgenomics.com/software/gctb/#Download)
## This requires MPI and two C++ libraries i.e. Eigen3 and Boost

mpirun -np 5 ~/gctb \ # np is the number of CPU cores used for the analysis
--bfile ~/plink_data \ # plink format input data
--pheno ~/phenotype_file \ # phenotype file including only subjects after quality control
--extract ~/post_qc_snp_list.txt \ # extract SNPs that passed quality control
--covar ~/cov_file \ # name of covariate file, categorical covariates were converted to dummy variables
--bayes S \ # Bayesian alphabet for the analysis
--pi 0.05 \ # initial value of pi (standard value) 
--hsq 0.5 \ # initial value of hsq (standard value)
--S 0 \ # initial value of S (standard value)
--chain-length 21000 \ # number of simulations
--burn-in 1000 \ # number of burn-in simulations
--out-freq 500 \ # frequency of output printing in the .log file
--seed 23 \ # set seed for reproducibility of the results
--out ~/dep_heritability # name of output
