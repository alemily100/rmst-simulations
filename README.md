# rmst-simulations
This repository contains the code used to simulate results for the paper **The use of restricted mean survival time to summarise treatment
effect - a simulation study**. 

## Background 
This project looked to evaluate the power and Type I error of the Restricted Mean Survival Time (RMST) endpoint when treatment effect is explained by multiple covariates. Operating characteristics of each estimand were evaluated when RMST was evaluated using a non-parametric, fully-specified parametric, and mis-specified parametric approach. 

## Description of R files
**functions_interaction.R** - code for section 4: functions required to generate trials and associated z-values for RMST estimands of treatment effect when efficacy is explained across multiple covariates. 
**simulations_interaction.R** - code for section 4: code required to access power and Type I error for the RMST estimands. Functions required to run this code are contained in `functions_interaction.R`.
