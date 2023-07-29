library(data.table)
library(survival)
library(splines)
library(haven)

### Analysis

filenames <- list.files("M:/External Users/KevinJos/output/ipw/steroids/", full.names = TRUE)[-1]
fnames <- list.files("M:/External Users/KevinJos/output/ipw/steroids/", full.names = FALSE)[-1]

## cox spline models
for (i in 1) {
  
  print(fnames[i])
  load(filenames[i])
  
  setDT(dat)
  dat$bene_id_drug <- paste(dat$bene_id, dat$drug_time)
  
  ## msm code
  medFit <- dat[!duplicated(bene_id_drug)]
  medFit <- medFit[order(medFit$bene_id, medFit$time0)]
  medFit$ipw_temp <- unsplit(lapply(split(with(medFit, ipw_med*ipw_censor), medFit$bene_id),
                                    cumprod), medFit$bene_id)

  dat_mod <- merge(dat, data.frame(medFit[,c("bene_id","drug_time","ipw_temp")]),
                    by = c("bene_id","drug_time"), all.x = TRUE)

  dat_mod$ipw <- dat_mod$ipw_temp*dat_mod$ipw_pm
  
  ## hamsm code
  # dat_mod <- dat
  # dat_mod$ipw <- dat_mod$ipw_censor*dat_mod$ipw_med*dat_mod$ipw_pm
  
  # truncation
  dat_mod$ipw.trunc <- dat_mod$ipw
  dat_mod$ipw.trunc[dat_mod$ipw > quantile(dat_mod$ipw, 0.95)] <- quantile(dat_mod$ipw, 0.95)
  dat_mod$ipw.trunc[dat_mod$ipw < quantile(dat_mod$ipw, 0.05)] <- quantile(dat_mod$ipw, 0.05)
  
  ## index time scale
  # dat_mod$time0 <- with(dat_mod,(time0 - age)*365.25)
  # dat_mod$time1 <- with(dat_mod,(time1 - age)*365.25)
  # 
  # dat_tmp0 <- setDT(subset(dat_mod, onMeds == 0))
  # dat_tmp0$id <- paste0(dat_tmp0$bene_id, "_0")
  # 
  # dat_tmp1 <- setDT(subset(dat_mod, onMeds == 1))
  # dat_meds <- dat_tmp1[dat_tmp1[,.I[time0 == min(time0)], by = bene_id]$V1][!duplicated(bene_id)]
  # dat_tmp1 <- merge(dat_tmp1, data.frame(begin_meds = dat_meds$time0, bene_id = dat_meds$bene_id), by = "bene_id", all.x = TRUE)
  # dat_tmp1$id <- paste0(dat_tmp1$bene_id, "_1")
  # dat_tmp1$time0 <- with(dat_tmp1, time0 - begin_meds)
  # dat_tmp1$time0 <- with(dat_tmp1, time0 - begin_meds)
  # dat_tmp1$begin_meds <- NULL
  # 
  # dat_mod <- rbind(dat_tmp0, dat_tmp1)
  # dat_mod <- setDF(setDT(dat_mod)[order(id, time0)])
  # 
  # rm(dat_tmp0, dat_tmp1, dat_meds); gc()
  
  ## spline cox model
  dat_mod <- setDF(dat_mod[order(dat_mod$bene_id, dat_mod$time0),])
  model_ns <- coxph(Surv(time0, time1, failed) ~ onMeds + pspline(pm_nomed, 4) + pspline(pm_med, 4) +
                       strata(qc_type) + cluster(bene_id), id = bene_id, weights = ipw.trunc,
                    data = dat_mod, robust = TRUE, model = TRUE)
  
  ## save models
  save(model_ns, file=paste0("M:/External Users/KevinJos/output/age_time/cox_spline/",fnames[i]))
  
  rm(model_ns, dat, dat_mod, medFit)
  gc()
  
}
