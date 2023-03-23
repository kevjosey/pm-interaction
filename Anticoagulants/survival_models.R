library(data.table)
library(survival)
library(splines)
library(haven)
library(timereg)

filenames <- list.files("M:/External Users/KevinJos/output/ipw/", full.names = TRUE)
fnames <- list.files("M:/External Users/KevinJos/output/ipw/", full.names = FALSE)

for (i in 1:length(fnames)){
  
  print(fnames[i])
  dat <- get(load(filenames[i]))
  setDT(dat)
  
  dat$bene_id_drug <- paste(dat$bene_id, dat$drug_time)
  
  #remove duplicate combinations of bene_id and drug_time
  # medFit <- dat[!duplicated(bene_id_drug)]
  #order by bene id and start date of interval
  # medFit <- medFit[order(medFit$bene_id, medFit$time0)]
  #take the cumulative product of the med*censor weights
  # medFit$ipw_temp <- unsplit(lapply(split(with(medFit, ipw_med*ipw_censor), medFit$bene_id), cumprod), medFit$bene_id)
  
  #add these cumulative weights to the dataset
  # dat_mod <- merge(dat, data.frame(medFit[, c("bene_id", "drug_time", "ipw_temp")]), by = c("bene_id", "drug_time"), all.x = TRUE)
  dat_mod <- dat
  
  #multiply the IPWs to get a final IPW
  dat_mod$ipw <- dat_mod$ipw_pm*dat_mod$ipw_med*dat_mod$ipw_censor
  
  # truncation
  dat_mod$ipw.trunc <- dat_mod$ipw
  dat_mod$ipw.trunc[dat_mod$ipw > quantile(dat_mod$ipw, 0.975)] <- quantile(dat_mod$ipw, 0.975)
  dat_mod$ipw.trunc[dat_mod$ipw < quantile(dat_mod$ipw, 0.025)] <- quantile(dat_mod$ipw, 0.025)
  
  #order by bene id and start date of interval
  dat_mod <- setDF(dat_mod[order(dat_mod$bene_id, dat_mod$time0), ])
  
  #fit the model
  aalen_model <- aalen(formula = Surv(time0, time1, failed) ~ onMeds*pm,
                       data = dat_mod, start.time = 65.1, max.time = 95,
                       clusters = dat_mod$bene_id, id = dat_mod$bene_id,
                       robust = 1, covariance = 1, n.sim = 1000)
  
  cox_model <- coxph(Surv(time0, time1, failed) ~ onMeds*pm + cluster(bene_id),
                     id = bene_id, weights = ipw.trunc, data = dat_mod, robust = TRUE, model = TRUE)
  
  #save the output
  save(aalen_model, file = paste0("M:/External Users/KevinJos/output/age_time/aalen/", fnames[i]))
  save(cox_model, file=paste0("M:/External Users/KevinJos/output/age_time/cox/",fnames[i]))
  
  rm(cox_model, aalen_model, dat, dat_mod, medFit)
  gc()
  
}
