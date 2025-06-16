# Hong Kong Omicron Epidemic Modeling 

This repository contains scripts and data to model the Omicron outbreak in Hong Kong, focusing on the description of the foundational model structure. It employs a deterministic SEIIR (Susceptible-Exposed-Infectious-Infectious-Recovered) model with observational processes to estimate the incidence of severe health outcomes, including hospitalizations, severe cases, deaths, and total infections. The model integrates vaccination status, age-specific dynamics, and social mixing patterns to provide a comprehensive framework for outbreak analysis.


## Repository Structure
### 1. Scripts
#### 1.1 Model_Inputs.R
This script initializes the model by:

Setting up demographic inputs (age groups, vaccination doses, population data).
Defining compartments for the SEIR model and their interactions.
Loading necessary data, such as population census and vaccination rates.
Preparing parameters

#### 1.2 Model_TransObs_ftns.R
This script contains the core functions for the epidemic model framework with the two major components the SEIR model (seiir.model) and observational processes (lines 604-682)

seiir.model: Implements the deterministic SEIR system with vaccination and infection dynamics.
TransRep: Combines the SEIR model with observational processes (severe disease delay pmfs) to simulate severe health outcomes (hospitalizations, severe cases, deaths) and outputs daily incidence rates.
Includes utility functions for age group conversions and slope calculations for severe outcome probabilities.

#### 1.3 Model_run.R
This script runs the model and visualizes the results:

Model Execution: Calls the initialization (Model_Inputs.R) and SEIR/observation functions (Model_TransObs_ftns.R) to simulate daily infections and severe health outcomes.
Visualization: Generates plots of daily incidences of infections, hospitalizations, severe cases, and deaths, stratified by age and vaccination status.

### 1.4 Sims/Sim_ftn.R
Creates EpiWaveSimulation function, which runs the epidemic model (TransRep) with a various inputs, balise coverages, vaccinations rates, or time_varying_iota. 

Outputs: The specific severe disease and infectino incidence.


### 1.5 Sims/Sims_VacBasChg.R

Sets up counterfactual simuations, requires  model inputs, vaccination data and other model inputs. Loops over ages, doses, coverage change using EpiWaveSimulation function
Ouput: Counterfactul results


### 1.6 Fig/Fig_RisksVEDelt.R

Plots a range of calibrated model parameters, baseline risk of severe disease, vaccine effectiveness, and the peak period amplificaiton factor\

### 1.7 Fig/Fig_BaseChange.R

Plots coutnerfactual simulations results from the Sims_VacBasChg.R
Output:geom_tile plot. 


### 2. Data
temp_data/: Contains input datasets required to run the model, such as:
Population census data (Clean_age_sex_10AgeBandsDemoHk_end21.csv).
Vaccination proportions (N_vacp.rds) and vaccination rates
Probabilities mass funcitons severe outcomes delay, for peak and non-peak periods (Peak_pmf_mat.rds, NonPeak_pmf_mat.rds).
Social mixing matrices (contact_school.rds, contact_work.rds, etc.).

### 3. Outputs
Figure of simulated daily incidences of infections, hospitalizations, severe cases, and deaths.
Counterfactual simuations for baseline vaccination coverages
Figures of severe dissease risk related parameters and counterfactual vaccination scenario results. 

## Counterfactual simulations
Counterfactual scerarios were simulated as part of the study to observe the impact of the vaccinations in Hong Kong
Starting vaccination coverage proportions +/- 80% in 10% increments were run to observe the impact of the differece
Timing, shifting the observed coverage rates +/- 8 weeks 

## Awknowlegdements
Some of the functions in this modelling framework were adapted from work published by Corbella et al. (2018) DOI: 10.1186/s12889-018-5671-7. We thank the authors for their contributions to open science.
