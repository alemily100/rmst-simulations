# rmst-simulations
This repository contains the code used to simulate results for the paper **The use of restricted mean survival time to estimate joint treatment effects in the presence of interaction terms under model misspecification - a simulation study**. 

## Background 
This project evaluates the power and type I error of the Restricted Mean Survival Time (RMST) estimator when treatment effect is explained by multiple covariates. Operating characteristics of each estimator were evaluated when RMST was calculated using a non-parametric, fully specified parametric, and misspecified parametric approach. 

## Description of R files
* **functions_missurvival.R** - functions required to generate trials and associated Z-statistics for RMST estimators of treatment effect when efficacy is explained across multiple covariates.
  
* **simulations_missurvival.R** - code required to assess power and type I error for the RMST estimator. Functions required to run this code are defined in `functions_missurvival.R`.

* **functions_crossing.R** - functions required to generate trials and associated Z-statistics for RMST estimators of treatment effect when survival curves are known to cross.
  
* **simulations_crossing.R** - code required to assess power and type I error for the RMST estimator with crossing survival curves. Functions required to run this code are defined in `functions_crossing.R`.
