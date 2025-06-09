#R script master thesis Deborah van Putten
#This script fits a Bayesian multilevel model to estimate genetic effects 
#on two life-history traits (pupal mass and development rate) in Helicoverpa armigera 
#exposed to different biopesticide doses.
#The goal is to estimate genetic variation in these traits and the genetic correlations
#within and between traits across different doses to find potential trade-offs and
#gene-by-environment interactions

#load libraries####
library(readxl)
library(rethinking)
library(dplyr) 
library(tidyverse)


#load data####
data <- read_excel("moth_exp_data_master_thesis.xlsx", sheet = 1)

##data preparation####
#filter data set for no dosage information
filtered_data <- data %>%
  filter(!is.na(Dosage))

#data set for development and pupal mass analysis
#filter out no pupation, no pupal weight, and malformation
filtered_data_after_mortality <- filtered_data %>%
  filter(!is.na(Pupation_date) & 
           !is.na(Pupal_weight) & 
           is.na(Malformation) )

#filter to only keep males
filtered_data_males <- filtered_data_after_mortality %>%
  filter(Pupae_sex == "m")


#model####
## data preparation: encoding dose and trait variables####

#assign a numeric value to each dose level
dose_mapping <- c("0" = 1, "0.0625" = 2, "0.125" = 3, "0.25" = 4, "0.5" = 5, "1" = 6)

filtered_data_males$dose_numeric <- as.integer(dose_mapping[as.character(filtered_data_males$Dosage)])
filtered_data_males$Pupal_weight <- as.numeric(filtered_data_males$Pupal_weight)


#assinging a number 1 to 12 depending on the dose and trait
#pupal mass (trait 1) trait-dose 1 to 6 depending on dose
#development rate (trait 2) trait-dose 7 to 12 depending on dose
filtered_data_males$trait_dose <- NA
filtered_data_males$trait_dose[filtered_data_males$Pupal_weight > 0] <- filtered_data_males$dose_numeric
filtered_data_males$trait_dose[filtered_data_males$Age_at_pupation > 0] <- filtered_data_males$dose_numeric + 6


##log transform and  traits####
# pupal mass and development rate are log transformed and standardized to allow for comparability
mass_data <- filtered_data_males %>%
  mutate(trait = 1,  # Trait 1 is pupal mass
         standardized_value_trait = standardize(log(Pupal_weight)),  # Standardized log(pupal mass)
         trait_dose = dose_numeric) %>%
  select(Sire, Dam, dose_numeric, standardized_value_trait, trait_dose, trait)

#data for development rate, log transform and standardizing
#development rate is calculated as the inverse of the amount of days of egg hatching ot pupation
dev_data <- filtered_data_males %>%
  mutate(trait = 2,  # trait 2 is development rate
         standardized_value_trait = standardize(log(1 / Age_at_pupation)),  
         trait_dose = dose_numeric + 6) %>%  # Shift dose for development rate (7-12)
  select(Sire, Dam, dose_numeric, standardized_value_trait, trait_dose, trait)

#combine both datasets into one 
combined_data <- bind_rows(mass_data, dev_data)


##data frame to put into the model####
larvae_filtered_data_male_mass_dev <- list(
  sire = as.integer(factor(combined_data$Sire)),  # Sire index
  dam = as.integer(factor(combined_data$Dam)),    # Dam index
  dose = as.integer(factor(combined_data$dose_numeric)),  # Numeric dose
  stand_trait = combined_data$standardized_value_trait,
  trait_dose = combined_data$trait_dose,  # trait_dose index (1 to 12)
  N_sires = length(unique(combined_data$Sire)),  # Number of unique sires
  N_dams = length(unique(combined_data$Dam)),  # Number of unique dams
  N_trait_dose = length(unique(combined_data$trait_dose))
)


##ulam model ####
#models the standardized trait values as normally distributed
m_male_stand <- ulam(
  alist(
    #likelihood standardized traits 
    stand_trait ~ dnorm(mu_trait, sigma),
    
    #linear model 
    mu_trait <- baseline + dose_effect[trait_dose] + sire_effect[sire, trait_dose] + dam_effect[dam],
    
    #fixed intercept
    baseline ~ dnorm(0, 0.100), 
    
    #fixed dose effect varying by trait
    dose_effect[trait_dose] ~ dnorm(0,0.1),
    
    #dam effect not varying by treatment
    dam_effect[dam] ~ dnorm(0,sigma_dam), 
    
    #prior SD dam effect
    sigma_dam ~ dexp(1),  
    
    #random slopes per sire across doses 
    transpars> matrix[N_sires, N_trait_dose]:sire_effect <-
      compose_noncentered(sigma_sire, L_rho_sire, z_sire),  
    matrix[N_trait_dose, N_sires]:z_sire ~ normal(0,1),  
    
    #priors for random effects
    vector[N_trait_dose]:sigma_sire ~ dexp(1),  # SD of slopes per dose
    cholesky_factor_corr[N_trait_dose]:L_rho_sire ~ lkj_corr_cholesky(1),  # Correlation between dose effects per sire
    
    #correlation matrix sire effect of dose
    gq> matrix[N_trait_dose, N_trait_dose]:rho_sire <<- Chol_to_Corr(L_rho_sire),
    
    #residual SD
    sigma ~ dexp(1)
  ),
  data = larvae_filtered_data_male_mass_dev,
  chains = 4, cores = 4 
)


##after running model diagnostics####
precis(m_male_stand, depth=3)
trankplot(m_male_stand)
traceplot_ulam(m_male_stand)


##extract posteriors####
#contains the posterior samples from the model including
#dose-effects, sire and dam effects, and the correlation matrix
post_mass_dev_male <- extract.samples(m_male_stand)