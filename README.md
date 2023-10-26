# rmst-simulations
This repository contains the code used to simulate results for the paper **The use of restricted mean survival time to estimate treatment
effect - a simulation study**. 

## Background 
This project evaluates the power and type I error of the Restricted Mean Survival Time (RMST) estimator when treatment effect is explained by multiple covariates. Operating characteristics of each estimator were evaluated when RMST was calculated using a non-parametric, fully specified parametric, and misspecified parametric approach. 

## Description of R files
* **functions_missurvival.R** - code for section 4: functions required to generate trials and associated Z-statistics for RMST estimators of treatment effect when efficacy is explained across multiple covariates.
  
* **simulations_missurvival.R** - code for section 4: code required to assess power and type I error for the RMST estimator. Functions required to run this code are defined in `functions_missurvival.R`.
