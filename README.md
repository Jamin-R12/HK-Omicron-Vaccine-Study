# Omicron Hong Kong Modeling Using SEIR Framework to study vaccination modelling

This repository contains basic scripts and data to simulate the Omicron outbreak using a deterministic SEIIR model with observational processes for severe health outcomes incidence. 
The model is designed to estimate hospitalizations, severe cases, and deaths incidence and the total infecitons, incorporating vaccination and social mixing dynamics.

## Repository Structure
### 1. Scripts
1.1 Model_Inputs.R
This script initializes the model by:

Setting up demographic inputs (age groups, vaccination doses, population data).
Defining compartments for the SEIR model and their interactions.
Loading necessary data, such as population census and vaccination rates.
Preparing parameters

1.2 Model_TransObs_ftns.R
This script contains the core functions for the SEIR model and observational processes:

seiir.model: Implements the deterministic SEIR system with vaccination and infection dynamics.
TransRep: Combines the SEIR model with observational processes to simulate severe health outcomes (hospitalizations, severe cases, deaths) and outputs daily incidence rates.
Includes utility functions for age group conversions and slope calculations for severe outcome probabilities.

1.3 Model_run.R
This script runs the model and visualizes the results:

Model Execution: Calls the initialization (Model_Inputs.R) and SEIR/observation functions (Model_TransObs_ftns.R) to simulate daily infections and severe health outcomes.
Visualization: Generates plots of daily incidences of infections, hospitalizations, severe cases, and deaths, stratified by age and vaccination status.

### 2. Data
temp_data/: Contains input datasets required to run the model, such as:
Population census data (Clean_age_sex_10AgeBandsDemoHk_end21.csv).
Vaccination proportions (N_vacp.rds) and vaccination rates
Probabilities mass funcitons severe outcomes delay, for peak and non-peak periods (Peak_pmf_mat.rds, NonPeak_pmf_mat.rds).
Social mixing matrices (contact_school.rds, contact_work.rds, etc.).

3. Outputs
Simulated daily incidences of infections, hospitalizations, severe cases, and deaths.

## Counterfactual simulations
Counterfactual scerarios were simulated as part of the study to observe the impact of the vaccinations in Hong Kong
Starting vaccination coverage proportions +/- 80% in 10% increments were run to observe the impact of the differece
Timing, shifting the observed coverage rates +/- 8 weeks 

## Awknowlegdements
Some of the functions in this modelling framework were adapted from work published by Corbella et al. (2018) DOI: 10.1186/s12889-018-5671-7. 
We thank the authors for their contributions to open science.
