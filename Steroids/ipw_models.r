# Load libraries
library(haven)
library(data.table)
library(lubridate)
library(stringr)

source("M:/External Users/KevinJos/code/ipw_fun.R")

# Specify the drug class to be analyzed
medclass <- "oralsteroid"

# Specify the outcome variables to be analyzed
outvar_all <- c('mi_acs','iscstroke_tia','newvte','newhf','fib','death')

filenames <- list.files("M:/External Users/KevinJos/data/steroids/", full.names = TRUE)
fnames <- list.files("M:/External Users/KevinJos/data/steroids/", full.names = FALSE)

# Start loop over all outcome variables
for (i in 1:length(outvar_all)){
  
  # Set the outcome variable for the current iteration
  outvar <- outvar_all[i]

  # Load the data for the current drug class and outcome variable
  cpFit <- get(load(file =  filenames[i]))
  setDT(cpFit)

  # Select variables that include 'dx' and 'hx' in their names
  select_dx <- names(cpFit)[grep("dx", names(cpFit))] 
  select_hx <- names(cpFit)[grep("hx", names(cpFit))]

  # IPW Analysis ------------------------------------------------------------

  # Prepare the data for PM IPW weights
  # Remove duplicates and set formula for denominator
  pmFit <- cpFit[!duplicated(zip_ssn)]
  
  # Set formula for the denominator in the propensity score model for pm
  fmla.denom.pm <- formula(paste0('~ pm.lag1+pm.lag2+pm.lag3+pm.lag4+
                                    season+yr_ssn+poverty+popdensity+
                                    medianhousevalue+pct_owner_occ+education+
                                    medhouseholdincome+pct_hisp+pct_blk+pct_white'))
  
  # Compute PM IPW weights using the ipwtm_xgb function
  pm_w <- ipwtm_xgb(exposure = "pm", denominator = fmla.denom.pm, 
                    id = "zip",timevar = "ssn_time", data = pmFit, trunc = NULL,
                    continuous = TRUE, num.trees = 500, max.depth = 6, num.threads = 8)
  
  save_pm <- data.table(pmFit[,c("zip","ssn_time")], ipw_pm = pm_w$ipw.weights)
  
  # Free up memory
  rm(pmFit, fmla.denom.pm, pm_w); gc()

  # Compute Steroid IPW weights using similar steps as above
  medFit <- cpFit[cpFit[,.I[onMeds == max(onMeds)], by = bene_id_drug]$V1]
  medFit <- medFit[!duplicated(bene_id_drug)][order(bene_id_drug)]
  medFit <- medFit[,onMeds.lag:=shift(onMeds, n = 1, type = "lag"), by = bene_id]
  medFit <- medFit[onMeds.lag == 0 | is.na(onMeds.lag)]
  
  rm(cpFit); gc()
  
  fmla.denom.med <- formula(paste0('~ pm+age_tm+sex+race+dualeligible+season+yr_ssn+poverty+
                                    popdensity+medianhousevalue+pct_owner_occ+education+
                                    medhouseholdincome+pct_hisp+pct_blk+pct_white+',
                                    paste0(select_hx, collapse = '+'), '+',
                                    paste0(select_dx, collapse = '+')))

  med_w <- ipwtm_xgb(exposure = "onMeds", denominator = fmla.denom.med, numerator = formula("~ pm"),
                     id = "bene_id",timevar = "drug_time", data = medFit, trunc = NULL, 
                     continuous = FALSE, censor = FALSE, num.trees = 200, max.depth = 6, num.threads = 8)

  save_med <- data.table(medFit[,c("bene_id","drug_time")], ipw_med = med_w$ipw.weights)
  
  rm(medFit, fmla.denom.med, med_w); gc()

  # Compute censor IPW weights using similar steps as above
  cpFit <- get(load(file =  filenames[i]))
  setDT(cpFit)
  
  if (outvar != "death")
    cpFit[, censored := as.numeric(cpFit$censored | cpFit$died)]
  
  censorFit <- cpFit[cpFit[,.I[censored == max(censored)], by = bene_id_drug]$V1][!duplicated(bene_id_drug)]

  rm(cpFit); gc()
  
  fmla.denom.censor <- formula(paste0('~ pm+onMeds+age_tm+sex+race+dualeligible+season+yr_ssn+poverty+
                                      popdensity+medianhousevalue+pct_owner_occ+education+
                                      medhouseholdincome+pct_hisp+pct_blk+pct_white+',
                                     paste0(select_hx, collapse = '+'), '+',
                                     paste0(select_dx, collapse = '+')))

  censor_w <- ipwtm_xgb(exposure = "censored", denominator = fmla.denom.censor, numerator = formula("~ onMeds*pm"),
                        id = "bene_id", timevar = "drug_time", data = censorFit, trunc = NULL,
                        continuous = FALSE, censor = TRUE, num.trees = 200, max.depth = 6, num.threads = 8)
  
  save_censor <- data.table(censorFit[,c("bene_id","drug_time")], ipw_censor = censor_w$ipw.weights)
  
  rm(censorFit, fmla.denom.censor, censor_w); gc()
  
  # Merge IPW weights back to the original dataset and handle NA's
  cpFit <- get(load(file = filenames[i]))
  setDT(cpFit)
               
  cpFit <- merge(cpFit, save_pm, by = c("zip","ssn_time"), all.x = TRUE)
  cpFit <- merge(cpFit, save_med, by = c("bene_id","drug_time"), all.x = TRUE)
  cpFit <- merge(cpFit, save_censor, by = c("bene_id","drug_time"), all.x = TRUE)

  cpFit[is.na(ipw_pm), ipw_pm := 1]
  cpFit[is.na(ipw_med), ipw_med := 1]
  cpFit[is.na(ipw_censor), ipw_censor := 1]

  # Prepare the dataset for saving
  dat <- subset(cpFit, select = c(bene_id, date, zip, season, yr_ssn, ssn_time, drug_time,
                                  qc_type, age, onMeds, pm, event, failed, died, censored,
                                  time0, time1, ipw_pm, ipw_med, ipw_censor))

  dat <- dat[order(bene_id,time0)]
  dat[, time0 := time0/365.25 + age]
  dat[, time1 := time1/365.25 + age]
  dat[, pm_med := ifelse(dat$onMeds == 1, dat$pm, 8)]
  dat[, pm_nomed := ifelse(dat$onMeds == 0, dat$pm, 8)]

  # Save the data
  save(dat, file=paste0("M:/External Users/KevinJos/output/ipw/steroids/", fnames[i], '.RData'))

  # Clear memory for the next iteration
  rm(dat, cpFit, save_censor, save_med, save_pm); gc()
  
}
