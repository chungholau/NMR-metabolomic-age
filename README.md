# NMR-metabolomic-age
This repository contains supporting materials (data/codes) for the manuscript "NMR metabolomic modelling of age and lifespan: a multi-cohort analysis".

## Introduction
Metabolomic age models have been proposed for the study of biological aging, however they have not been widely validated. We aimed to assess the performance of newly developed and existing nuclear magnetic resonance spectroscopy (NMR) metabolomic age models for prediction of chronological age, mortality, and age-related disease. New statistical models of biological age were created in R using NMR metabolomics data. These metabolomics models may be applied to new NMR datasets similarly generated by Brainshake Ltd./Nightingale Health.

## File description
**metaboAgeDB_models.Rdata**: This R file contains the metabolomic age model objects which are required to generate new model predictions. This file includes four new models that were generated and presented in the manuscript: "Elastic Net", "MARS", "Study mortality score", and "phenotypic ageing".

**pickUpUKBSampleTOCalibrate.R**: This R script was used to match and calibrate metabolic variables in UK Biobank to Whitehall II based on matching on demongraphic characteristics (age, sex and BMI). This can be used as an example script for calibrating metabolic variables in any new dataset to our study reference dataset (Whitehall II)

**script_generate_scores_UKB.R**: This R script was used to generate model ageing scores in the UK Biobank data. This can be used as an example script for generating new model predictions in other nightingale dataset.

**nightingale_name_conversion.csv**: This csv file gives information on how metabolic variable names were converted between cohorts/ dataset versions.

## Instruction to use/Contact

1. In order to apply these metabolomics models generated in the R environment, users first need to check their variable names are identical to those listed under the "UKB" column in the *nightingale_name_conversion.csv*.
2. Cohort calibration / matching will be needed prior to generated the model predictions, expecially on Nightigale data generated post-2020. This can be performed using *pickUpUKBSampleTOCalibrate.R* as a template.
3. Once *metaboAgeDB_models.Rdata* is saved in the working directory, Metabolomic ageing scores can be generated by using *script_generate_scores_UKB.R* as a template

NB. Some adaptations in the script/ additional data wrangling will likly be required to make it work for your datasets, and the scripts and information provided should be used as a guide.

Please direct any enquires to chungho.lau@imperial.ac.uk
Date: 06/11/2023
