### This repository includes the R code used to create and analyse the phenotypes major depressive disorder and treatment-resistant depression in the study: 
### Fabbri et al. 2020. Genetic and clinical characteristics of treatment-resistant depression using primary care records in two UK cohorts. medRxiv 
### available at doi: https://doi.org/10.1101/2020.08.24.20178715


### Description of contents:

## scripts folder: includes the following R scripts:

	# extract_diagn_ADs_TRD_pheno.R: R code to: 1) extract and annotate diagnoses of depression, bipolar disorder, psychotic disorders and substance use disorders; 2) extract and annotate antidepressant medications; 3) create the phenotype major depressive disorder; 4) create the phenotype treatment-resistant depression

	# extract_other_diagn_med.R: R code to: 1) extract and annotate other psychiatric diagnoses; 2) extract and annotate other psychotropic drugs prescriptions; 3) identify combinations of antidepressants, antidepressants-antipsychotics and antidepressants-mood stabilisers prescribed for >= 30 days

	# PRS_regress_liability_scale.R: R code to analyse polygenic risk scores (PRS) and convert observed Nagelkerke R2 to the liability scale

	# SNP_h2_liability_scale.R: R code for conversion of observed heritability to the liability scale

## heritability_est_GCTA_GCTB.sh:
	# code used to estimate SNP-based heritability
