library(haven)
library(data.table)
library(lubridate)
library(stringr)
library(parallel)
library(geepack)
library(xgboost)
library(ranger)
library(survival)

source("M:/External Users/KevinJos/code/ipw_fun.R")
source("M:/External Users/KevinJos/code/interaction.R")

filenames <- list.files("M:/External Users/KevinJos/data/steroids/", full.names = TRUE)
fnames <- list.files("M:/External Users/KevinJos/data/steroids/", full.names = FALSE)

# get the full name of the drug in the dataset (medclass)
medclass <- "oralsteroid"
outvar_all <- c('mi_acs','iscstroke_tia','newvte','newhf','fib','death')
boot.iter <- 100

bootstrap_data <- function(data, index, u.zip) {
  
  n.zip <- length(u.zip)
  boot <- data.frame()
  
  aa <- u.zip[index]
  aa <- aa[which(aa %in% data$zip)]
  bb <- table(aa)
  
  for (j in 1:max(bb)) {
    
    cc <- data[data$zip %in% names(bb[bb == j]),]
    
    for (k in 1:j) {
      
      cc$boot.bene_id <- paste(cc$bene_id, k, sep = "-")
      cc$boot.bene_id_drug <- paste(cc$bene_id_drug, k, sep = "-")
      boot <- rbind(boot, cc)
      
    }
    
  }
  
  return(boot)
  
}

# For loop for each outcome

for (i in 1:length(outvar_all)){
  
  outvar <- outvar_all[i]
  cpFit <- get(load(file =  filenames[i]))
  
  select_dx <- names(cpFit)[grep("dx", names(cpFit))] 
  select_hx <- names(cpFit)[grep("hx", names(cpFit))]
  
  u.zip <- unique(cpFit$zip)
  m <- length(u.zip)/log(length(u.zip))

  boot_list <- mclapply(1:boot.iter, function(b, ...) {
    
    print(b)
    
    index <- sample(1:length(u.zip), m, replace = TRUE)
    cpBoot <- bootstrap_data(data = cpFit, index = index, u.zip = u.zip)
    
    ## PM IPW WEIGHTS
    setDT(cpBoot)
    pmFit <- cpBoot[!duplicated(zip_ssn)]
    
    fmla.denom.pm<-formula(paste0('~ pm.lag1+pm.lag2+pm.lag3+pm.lag4+
                                    season+yr_ssn+poverty+popdensity+
                                    medianhousevalue+pct_owner_occ+education+
                                    medhouseholdincome+pct_hisp+pct_blk+pct_white'))
    
    pm_w <- ipwtm_xgb(exposure = "pm", denominator = fmla.denom.pm,
                     id = "zip",timevar = "ssn_time", data = pmFit, trunc = NULL, 
                     continuous = TRUE, num.trees = 50, max.depth = 8, num.threads = 16)
    
    cpBoot <- merge(cpBoot, data.frame(pmFit[,c("zip","ssn_time")], ipw_pm = pm_w$ipw.weights),
                   by = c("zip","ssn_time"), all.x = TRUE)
    
    cpBoot$ipw_pm[is.na(cpBoot$ipw_pm)] <- 1
    
    rm(fmla.denom.pm, pmFit, pm_w) 
    
    ## STEROID IPW WEIGHTS
    
    medFit <- cpBoot[cpBoot[,.I[onMeds == max(onMeds)], by = boot.bene_id_drug]$V1]
    medFit <- medFit[!duplicated(boot.bene_id_drug)]
    medFit <- medFit[order(boot.bene_id_drug)]
    medFit <- medFit[,onMeds.lag:=shift(onMeds, n = 1, type = "lag"), by = boot.bene_id]
    medFit <- subset(medFit, onMeds.lag == 0 | is.na(onMeds.lag))
    
    fmla.denom.med <- formula(paste0('~ pm+age_tm+sex+race+dualeligible+season+yr_ssn+poverty+
                                    popdensity+medianhousevalue+pct_owner_occ+education+
                                    medhouseholdincome+pct_hisp+pct_blk+pct_white+',
                                     paste0(select_hx, collapse = '+'), '+',
                                     paste0(select_dx, collapse = '+')))
    
    med_w <- ipwtm_xgb(exposure = "onMeds", denominator = fmla.denom.med, numerator = formula("~ pm"),
                      id = "boot.bene_id",timevar = "drug_time", data = medFit, trunc = NULL, 
                      continuous = FALSE, censor = FALSE, num.trees = 50, max.depth = 8, num.threads = 16)
    
    cpBoot <- merge(cpBoot, data.frame(medFit[,c("boot.bene_id","drug_time")], ipw_med = med_w$ipw.weights),
                    by = c("boot.bene_id","drug_time"), all.x = TRUE)
    
    cpBoot$ipw_med[is.na(cpBoot$ipw_med)] <- 1
    
    rm(med_w, fmla.denom.med, medFit)
    
    ## CENSOR IPW WEIGHTS
    
    if (outvar != "death")
      cpBoot$censored <- as.numeric(cpBoot$censored | cpBoot$died)
    
    censorFit <- cpBoot[cpBoot[,.I[censored == max(censored)], by = boot.bene_id_drug]$V1]
    censorFit <- censorFit[!duplicated(boot.bene_id_drug)]
    
    fmla.denom.censor <- formula(paste0('~ pm+onMeds+age_tm+sex+race+dualeligible+season+yr_ssn+poverty+
                                      popdensity+medianhousevalue+pct_owner_occ+education+
                                      medhouseholdincome+pct_hisp+pct_blk+pct_white+',
                                        paste0(select_hx, collapse = '+'), '+',
                                        paste0(select_dx, collapse = '+')))
    
    censor_w <- ipwtm_xgb(exposure = "censored", denominator = fmla.denom.censor, numerator = formula("~ onMeds*pm"),
                          id = "boot.bene_id", timevar = "drug_time", data = censorFit, trunc = NULL,
                          continuous = FALSE, censor = TRUE, num.trees = 50, max.depth = 8, num.threads = 16)
    
    cpBoot <- merge(cpBoot, data.frame(censorFit[,c("boot.bene_id","drug_time")], ipw_censor = censor_w$ipw.weights),
                    by = c("boot.bene_id","drug_time"), all.x = TRUE)
    
    cpBoot$ipw_censor[is.na(cpBoot$ipw_censor)] <- 1
    
    rm(censor_w, fmla.denom.censor, censorFit) 
    
    # put it all together
    cpBoot <- cpBoot[order(bene_id,time0)]
    cpBoot[, time0 := time0/365.25 + age]
    cpBoot[, time1 := time1/365.25 + age]
    cpBoot[, pm_med := ifelse(cpBoot$onMeds == 1, cpBoot$pm, 8)]
    cpBoot[, pm_nomed := ifelse(cpBoot$onMeds == 0, cpBoot$pm, 8)]
    
    ## msm code
    
    medFit <- cpBoot[!duplicated(bene_id_drug)]
    medFit <- medFit[order(medFit$bene_id, medFit$time0)]
    medFit$ipw_temp <- unsplit(lapply(split(with(medFit, ipw_med*ipw_censor), medFit$bene_id),
                                      cumprod), medFit$bene_id)
    
    cpBoot <- merge(cpBoot, data.frame(medFit[,c("bene_id","drug_time","ipw_temp")]),
                     by = c("bene_id","drug_time"), all.x = TRUE)
    
    cpBoot$ipw <- cpBoot$ipw_temp*cpBoot$ipw_pm
    rm(medFit)
    
    # truncation
    cpBoot$ipw.trunc <- cpBoot$ipw
    cpBoot$ipw.trunc[cpBoot$ipw > quantile(cpBoot$ipw, 0.975)] <- quantile(cpBoot$ipw, 0.975)
    cpBoot$ipw.trunc[cpBoot$ipw < quantile(cpBoot$ipw, 0.025)] <- quantile(cpBoot$ipw, 0.025)
    
    ## spline cox model
    cpBoot <- setDF(cpBoot[order(cpBoot$bene_id, cpBoot$time0),])
    model_ns <- coxph(Surv(time0, time1, failed) ~ onMeds + pspline(pm_nomed, 4) + pspline(pm_med, 4) +
                        strata(qc_type) + cluster(bene_id), id = bene_id, weights = ipw.trunc, 
                      data = cpBoot, robust = TRUE, model = TRUE)
    
    # outcome coefficients 
    coef_out <- data.frame(outcome = outvar[i], 
                           coefs = summary(model_ns)$coefficients,
                           vars = rownames(summary(model_ns)$coefficients),
                           iter = b)
    
    # reri components
    hr_tmp1 <- pm_contrast(model_ns, pm0 = 5, pm1 = 10)
    hr_tmp2 <- pm_contrast(model_ns, pm0 = 8, pm1 = 12)
    
    hr_out <- data.frame(rbind(hr_tmp1, hr_tmp2))
    hr_out$descript <- rep(c('Corticosteroid Use','Increasing PM while off Corticosteroids','Increasing PM while on Corticosteroids'), 2)
    hr_out$outcome <- outvar[i]
    hr_out$iter <- b
    
    reri_tmp1 <- data.frame(t(sapply(seq(8, 13, by = 0.1), function(z, ...)
      add_interact_cox(model_ns, pm0 = 8, pm1 = z))))
    reri_tmp2 <- data.frame(t(sapply(seq(5, 13, by = 0.1), function(z, ...)
      add_interact_cox(model_ns, pm0 = 5, pm1 = z))))
    
    reri_out <- data.frame(rbind(reri_tmp1, reri_tmp2))
    reri_out <- reri_out[order(reri_out$pm0, reri_out$pm1),]
    
    reri_out$outcome <- outvar[i]
    reri_out$iter <- b
    
    colnames(reri_out) <- c("est","lower","upper","pm0","pm1","outcome","iter")
    
    rm(model_ns)
    return(list(coef = coef_out, reri = reri_out, hr = hr_out))
    
  }, mc.cores = 16)
  
  coefs <- do.call(rbind, lapply(boot_list, function(item, ...) item$coefs))
  reri <- do.call(rbind, lapply(boot_list, function(item, ...) item$est))
  hr <- do.call(rbind, lapply(boot_list, function(item, ...) item$hr))
  
  coefs_se <- coefs %>% group_by(outcome, vars) %>% summarise(boot_se = (m/length(u.zip))*sd(coefs))
  hr_se <- coefs %>% group_by(outcome, descript, pm0) %>% summarise(boot_se = (m/length(u.zip))*sd(log.hr))
  reri_se <- coefs %>% group_by(outcome, pm0, pm1) %>% summarise(boot_se = (m/length(u.zip))*sd(est))
  
  out_list <- list(coefs = coefs, coefs_se = coefs_se,
                   reri = reri, reri_se = reri_se, 
                   hr = hr, hr_se = hr_se)
  
  ## save models
  save(out_list, file=paste0("M:/External Users/KevinJos/output/age_time/cox_spline/boot/", fnames[i]))
  
  rm(boot_list, out_list); gc()
  
}
