##################################################################################
################# heritability conversion to liability scale #####################
##################################################################################

### this script was used for the analyses of the study: https://doi.org/10.1101/2020.08.24.20178715

#### R code to convert observed heritability to the liability scale
#### see Yap et al. 2018, https://doi.org/10.1038/s41467-018-04807-3

h2liab_twothresholds <- function(Ku, Kl, h2o, p) {
  #Ku = population prevalence of cases
  #Kl = population prevalence of controls
  #h2o = observed scale heritability
  #p = sample prevalence
  #case_threshold = the position of the case threshold on the standard normal distribution
  #control_threshold = the position of the control threshold on the standard normal distribution
  #Zu = height of the standard normal distribution at the case prevalence
  #Zl = height of the standard normal distribution at the control prevalence
  #C = adjustment factor
  case_threshold = qnorm(1-Ku)

  control_threshold = qnorm(Kl)

  Zu = dnorm(case_threshold)

  Zl = dnorm(control_threshold)

  C = p*(1-p)*(((Zu/Ku)+(Zl/Kl))^2)

  h2l = h2o/C
}
