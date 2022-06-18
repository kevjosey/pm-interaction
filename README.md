# Effects of PM2.5 and Corticosteroid Use on Cardiovascular and Thromboembolic Events Among Older Adults: Evidence of Drug-Environment Interaction

This repository contains the code used to fit history adjusted marginal structural Cox models with a time to event outcome. We use this model to estimate the risks to cardiovascular and thromboembolic events (CTE) caused by exposre to fine particulate matter (PM2.5) and corticosteroids (steroids) for high-risk older adults enrolled in Medicare services. Below is a detailed description of the various scripts, which we have divided into two folders. The Functions folder contains 

#### Functions
[`Functions/interaction.R`](https://github.com/xiaodan-zhou/covid_wildfire/tree/master/src/bayes) include the BH-ZINB-DL model. 

[`Functions/ipw_fun.R`](https://github.com/xiaodan-zhou/covid_wildfire/tree/master/src/bayes) include the BH-ZINB-DL model. 

#### Analysis
[`Analysis/descriptives.R`](https://github.com/xiaodan-zhou/covid_wildfire/tree/master/src/bayes) include the BH-ZINB-DL model. 

[`Analysis/ipw_models.R`](https://github.com/xiaodan-zhou/covid_wildfire/tree/master/src/bayes) include the BH-ZINB-DL model. 

[`Analysis/survival_models.R`](https://github.com/xiaodan-zhou/covid_wildfire/tree/master/src/bayes) include the BH-ZINB-DL model. 

[`Analysis/aalen_output.R`](https://github.com/xiaodan-zhou/covid_wildfire/tree/master/src/bayes) include the BH-ZINB-DL model. 

[`Analysis/cox_output.R`](https://github.com/xiaodan-zhou/covid_wildfire/tree/master/src/bayes) include the BH-ZINB-DL model. 

#### Data
For data privacy reasons, the Medicare data used in this study cannot be made publicly available, but interested parties can request access by applying through the US Centers for Medicare and Medicaid Services. The PM2.5 exposure data are publicly available at the following link: https://doi.org/10.7927/0rvr-4538. Area-level covariates used herein are also publicly available from the US Census Bureau website.
