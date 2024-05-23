# Evaluating PM2.5 Interactions

This repository contains code used to fit history adjusted marginal structural Cox models with a time-to-event outcome. We use this model to estimate the risks to cardiovascular and thromboembolic events (CTE) caused by exposure to fine particulate matter (PM2.5) and either corticosteroids ([`Steroids`](https://github.com/kevjosey/pm-interaction/blob/main/Steroids/)) or Anticoagulants ([`Anticoagulants`](https://github.com/kevjosey/pm-interaction/blob/main/Anticoagulants/)) among high-risk older adults enrolled in Medicare. 

Detailed descriptions of the various scripts are available within their respective folders. The R files in the root directory contain the necessary functions for constructing the inverse probability weights and for evaluating additive interactions with either a Cox proportional hazards model or an Aalen additive hazards model.

### Functions

[`interaction.R`](https://github.com/kevjosey/pm-interaction/blob/main/Functions/interaction.R) Functions for evaluating additive and multiplicative interactions. 

[`ipw_fun.R`](https://github.com/kevjosey/pm-interaction/blob/main/Functions/ipw_fun.R) Functions for generating inverse probability weights (with Random Forest, XGBoost, or GEE). 

