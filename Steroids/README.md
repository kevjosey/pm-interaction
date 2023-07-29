# Effects of PM2.5 and Corticosteroid Use on Cardiovascular and Thromboembolic Events Among Older Adults: Evidence of Drug-Environment Interaction

[`Steroids/descriptives.R`](https://github.com/kevjosey/pm-steroid/blob/main/Steroids/descriptives.R) Contains code for finding descriptive statistics that populate Tables 1, 2, and 3. 

`Steroids/data_cleaning.R`](https://github.com/kevjosey/pm-steroid/blob/main/Steroids/ipw_models.R) Code for transforming data (wide-format) into a counting process with season/drug-quarter intervals (long-format).

[`Steroids/ipw_models.R`](https://github.com/kevjosey/pm-steroid/blob/main/Steroids/ipw_models.R) Code used for fitting the gradient boosted forests which model the mean of the (generalized) propensity scores and are subsequently used to construct the inverse probability weights.

[`Steroids/cox_models.R`](https://github.com/kevjosey/pm-steroid/blob/main/Steroids/cox_models.R) Code for fitting the proportional hazard models relating PM2.5 and steroids to each of the six health outcomes that we examined in the manuscript.

[`Steroids/bootstrap.R`](https://github.com/kevjosey/pm-steroid/blob/main/Steroids/bootstrap.R) Fits m-out-of-n bootstrap SEs.

[`Steroids/output.R`](https://github.com/kevjosey/pm-steroid/blob/main/Steroids/output.R) Creates tables and figures from the fitted proportional hazard model.

### Data

For data privacy reasons, the Medicare data used in this study cannot be made publicly available, but interested parties can request access by applying through the US Centers for Medicare and Medicaid Services. The PM2.5 exposure data are publicly available at the following link: https://doi.org/10.7927/0rvr-4538. Area-level covariates used herein are also publicly available from the US Census Bureau website.
