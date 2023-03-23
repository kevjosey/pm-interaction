library(haven)
library(data.table)
library(lubridate)
library(stringr)
library(ranger)
library(geepack)
library(survival)
library(plyr)
library(table1)

source("M:/External Users/KevinJos/code/ipw_fun.R")

medclass <- c("oralantic", "oralantip")
outvar_all <- c("bleeding_allgi", "bleeding_intracranial", "epistaxis", "death") 

for (i in 2:length(outvar_all)){
  
  for (j in 2) { # only looking at anticoagulants for now.
    
    # read the dataset
    med <- medclass[j]
    outvar <- outvar_all[i]
    cpFit <- get(load(file =  paste0("M:/External Users/KevinJos/data/", med, "_", outvar, ".RData")))
    
    select_dx <- names(cpFit)[grep("dx", names(cpFit))] 
    select_hx <- names(cpFit)[grep("hx", names(cpFit))]
    
    ################################### PM IPW WEIGHTS ######################################
    setDT(cpFit)
    pmFit <- cpFit[!duplicated(zip_ssn)] # remove duplicates
    setDF(pmFit) # make dataframe
    gc()
    
    # formula for denominator of IPW weights -- adjusting for neighborhood level confounders
    fmla.denom.pm <- formula("~ pm.lag1 + pm.lag2 + pm.lag3 + pm.lag4 + season + yr_ssn +
                             poverty + popdensity + medianhousevalue + pct_owner_occ +
                             education + medhouseholdincome + pct_hisp + pct_blk + pct_white")
    
    # use lag variables because we are looking at effect up until now, so we are including the history
    fmla.num.pm <- formula("~ pm.lag1 + pm.lag2 + pm.lag3 + pm.lag4")
    
    # run random forest model to get the weights
    pm_w <- ipwtm_ranger(exposure = "pm", denominator = fmla.denom.pm, numerator = fmla.num.pm,
                         id = "zip", timevar = "ssn_time", data = pmFit, trunc = NULL, continuous = TRUE,
                         num.trees = 200, max.depth = 8, num.threads = 16)
    
    # add the weights to the dataset
    cpFit <- merge(cpFit, data.frame(pmFit[, c("zip", "ssn_time")], ipw_pm = pm_w$ipw.weights), by = c("zip", "ssn_time"), all.x = TRUE)
    
    # if the weights are NA, make them equal to 1
    # there shouldn't be any NAs
    cpFit$ipw_pm[is.na(cpFit$ipw_pm)] <- 1
    
    rm(fmla.denom.pm, fmla.num.pm, pmFit, pm_w)
    gc()
    
    
    ################################ ANTIC/ANTIP IPW WEIGHTS ################################################
    
    setDT(cpFit)
    
    # for people who are on medication
    medFit <- cpFit[cpFit[, .I[onMeds == max(onMeds)], by = bene_id_drug]$V1]
    medFit <- medFit[!duplicated(bene_id_drug)]
    medFit <- medFit[order(bene_id_drug)]
    # create lag variable for medication
    medFit <- medFit[, onMeds.lag := shift(onMeds, n = 1, type = "lag"), by = bene_id]
    # subsetting by those who are not on medication
    medFit <- subset(medFit, onMeds.lag == 0 | is.na(onMeds.lag)) 
    newuser <- paste0("newuser_",medclass[j],"_365d")
    setDF(medFit)
    gc()
    
    # formula for denominator -- adjusting for neighborhood and individual level confounders
    fmla.denom.med <- formula(paste0("~pm + age_tm + sex + race + dualeligible + season + 
                                     yr_ssn + poverty + popdensity + medianhousevalue + pct_owner_occ +
                                     education + medhouseholdincome + pct_hisp + pct_blk + pct_white +",
                                     paste0(select_hx, collapse = "+"), "+", 
                                     paste0(select_dx, collapse = "+")))
    
    #formula for numerator
    fmla.num.med <- formula("~pm")
    
    #run RF model to get the IPW weights for medication
    med_w <- ipwtm_ranger(exposure = "onMeds", denominator = fmla.denom.med, numerator = fmla.num.med, 
                          id = "bene_id", timevar = "drug_time", data = medFit, trunc = NULL, continuous = FALSE, 
                          num.trees = 200, max.depth = 8, num.threads = 8)
    
    # adding the weights to the dataset
    cpFit <- merge(cpFit, data.frame(medFit[,c("bene_id", "drug_time")], ipw_med = med_w$ipw.weights), by = c("bene_id", "drug_time"), all.x = TRUE)
    
    cpFit$ipw_med[is.na(cpFit$ipw_med)] <- 1
    rm(fmla.denom.med, fmla.num.med, medFit, med_w)
    gc()
    
    ################################ CENSORING IPW WEIGHTS #######################################
    
    setDT(cpFit)
    
    if (outvar != "death")
      cpFit$censored <- as.numeric(cpFit$censored | cpFit$died)
    
    censorFit <- cpFit[cpFit[,.I[censored == max(censored)], by = bene_id_drug]$V1]
    censorFit <- censorFit[!duplicated(bene_id_drug)]
    setDF(censorFit); gc()
    
    fmla.denom.censor <- formula(paste0('~ pm+onMeds+age_tm+sex+race+dualeligible+season+yr_ssn+poverty+
                                      popdensity+medianhousevalue+pct_owner_occ+education+
                                      medhouseholdincome+pct_hisp+pct_blk+pct_white+',
                                        paste0(select_hx, collapse = '+'), '+',
                                        paste0(select_dx, collapse = "+")))
    
    censor_w <- ipwtm_ranger(exposure = "censored", denominator = fmla.denom.censor, numerator = formula("~ onMeds+pm"),
                             id = "bene_id", timevar = "drug_time", data = censorFit, trunc = NULL,
                             continuous = FALSE, censor = TRUE, num.trees = 200, max.depth = 8, num.threads = 8)
    
    cpFit <- merge(cpFit, data.frame(censorFit[,c("bene_id","drug_time")], 
                                     ipw_censor = censor_w$ipw.weights),
                   by = c("bene_id","drug_time"), all.x = TRUE)
    
    cpFit$ipw_censor[is.na(cpFit$ipw_censor)] <- 1
    
    rm(fmla.denom.censor, fmla.num.censor, censorFit, censor_w); gc()
    
    # selecting only the variables that we need to fit the survival models
    dat <- setDT(subset(cpFit, select = c(bene_id, zip, ssn_time, drug_time, age,
                                          onMeds, pm, event, censored, died, failed,
                                          time0, time1, ipw_pm, ipw_med, ipw_censor)))
    
    rm(cpFit)
    gc()
    
    dat <- dat[order(bene_id, time0)]
    dat$time0 <- with(dat, time0/365.25 + age) # using age as timescale
    dat$time1 <- with(dat, time1/365.25 + age)
    dat$pm_med <- ifelse(dat$onMeds == 1, dat$pm, 8) 
    dat$pm_nomed <- ifelse(dat$onMeds == 0, dat$pm, 8)
    
    save(dat, file=paste0("M:/External Users/KevinJos/output/ipw/", med, '_', outvar, '.RData'))
    
    rm(dat)
    gc()
    
  }
  
}