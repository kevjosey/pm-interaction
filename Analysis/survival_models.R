library(data.table)
library(survival)
library(splines)
library(timereg)

### Analysis

filenames <- list.files("M:/External Users/KevinJos/output/fit_data/", full.names = TRUE)
fnames <- list.files("M:/External Users/KevinJos/output/fit_data/", full.names = FALSE)

# aalen model
for (i in 1:length(fnames)) {
  
  print(fnames[i])
  load(filenames[i])
  
  ## msm code
  dat$bene_id_drug <- paste(dat$bene_id, dat$drug_time)

  if(is.null(dat$ipw_death)) {

    medFit <- dat[!duplicated(bene_id_drug)]
    medFit <- medFit[order(medFit$bene_id, medFit$time0)]
    medFit$ipw_temp <- unsplit(lapply(split(with(medFit, ipw_med), medFit$bene_id), cumprod), medFit$bene_id)

    dat <- merge(dat, data.frame(medFit[,c("bene_id","drug_time","ipw_temp")]),
                 by = c("bene_id","drug_time"), all.x = TRUE)

  } else {

    medFit <- dat[!duplicated(bene_id_drug)]
    medFit <- medFit[order(medFit$bene_id, medFit$time0)]
    medFit$ipw_temp <- unsplit(lapply(split(with(medFit, ipw_med*ipw_death), medFit$bene_id), cumprod), medFit$bene_id)

    dat <- merge(dat, data.frame(medFit[,c("bene_id","drug_time","ipw_temp")]),
                 by = c("bene_id","drug_time"), all.x = TRUE)

  }

  dat$ipw <- dat$ipw_temp*dat$ipw_pm
  
  ## truncate
  dat$weights.trunc <- dat$ipw
  dat$weights.trunc[dat$ipw < quantile(dat$ipw, 0.025)] <- quantile(dat$ipw, 0.025)
  dat$weights.trunc[dat$ipw > quantile(dat$ipw, 0.975)] <- quantile(dat$ipw, 0.975)
  
  ## medication time scale
  # dat$time0 <- with(dat,(time0 - age)*365.25)
  # dat$time1 <- with(dat,(time1 - age)*365.25)
  # 
  # dat_tmp0 <- setDT(subset(dat, onMeds == 0))
  # dat_tmp0$id <- paste0(dat_tmp0$bene_id, "_0")
  # 
  # dat_tmp1 <- setDT(subset(dat, onMeds == 1))
  # dat_meds <- dat_tmp1[dat_tmp1[,.I[time0 == min(time0)], by = bene_id]$V1][!duplicated(bene_id)]
  # dat_tmp1 <- merge(dat_tmp1, data.frame(begin_meds = dat_meds$time0, bene_id = dat_meds$bene_id), by = "bene_id", all.x = TRUE)
  # dat_tmp1$id <- paste0(dat_tmp1$bene_id, "_1")
  # dat_tmp1$time0 <- with(dat_tmp1, time0 - begin_meds)
  # dat_tmp1$time1 <- with(dat_tmp1, time1 - begin_meds)
  # dat_tmp1$begin_meds <- NULL
  # 
  # dat <- rbind(dat_tmp0, dat_tmp1)
  # dat <- setDF(setDT(dat)[order(id, time0)])
  # 
  # rm(dat_tmp0, dat_tmp1, dat_meds);gc()
  
  ## additive hazard model
  dat <- setDF(dat[order(dat$bene_id, dat$time0),])
  aalen_model <- aalen(formula = Surv(time0, time1, failed) ~ onMeds + pm + onMeds:pm,
                       data = dat, start.time = 66, max.time = 95, clusters = dat$bene_id, id = dat$bene_id, 
                       weights = dat$weights.trunc, robust = 1, covariance = 1, n.sim = 1000)
  
  ## save models
  save(aalen_model, file=paste0("M:/External Users/KevinJos/output/age_time/aalen_msm/",fnames[i]))
  
  rm(aalen_model, dat)
  gc()
  
}

# cox spline models
for (i in 1:length(fnames)) {
  
  print(fnames[i])
  load(filenames[i])
  
  ## msm code
  dat$bene_id_drug <- paste(dat$bene_id, dat$drug_time)

  if(is.null(dat$ipw_death)) {

    medFit <- dat[!duplicated(bene_id_drug)]
    medFit <- medFit[order(medFit$bene_id, medFit$time0)]
    medFit$ipw_temp <- unsplit(lapply(split(with(medFit, ipw_med), medFit$bene_id), cumprod), medFit$bene_id)

    dat <- merge(dat, data.frame(medFit[,c("bene_id","drug_time","ipw_temp")]),
                 by = c("bene_id","drug_time"), all.x = TRUE)

  } else {

    medFit <- dat[!duplicated(bene_id_drug)]
    medFit <- medFit[order(medFit$bene_id, medFit$time0)]
    medFit$ipw_temp <- unsplit(lapply(split(with(medFit, ipw_med*ipw_death), medFit$bene_id), cumprod), medFit$bene_id)

    dat <- merge(dat, data.frame(medFit[,c("bene_id","drug_time","ipw_temp")]),
                 by = c("bene_id","drug_time"), all.x = TRUE)

  }

  dat$ipw <- dat$ipw_temp*dat$ipw_pm
  
  ## truncate
  dat$weights.trunc <- dat$ipw
  dat$weights.trunc[dat$ipw < quantile(dat$ipw, 0.025)] <- quantile(dat$ipw, 0.025)
  dat$weights.trunc[dat$ipw > quantile(dat$ipw, 0.975)] <- quantile(dat$ipw, 0.975)
  
  ## medication time scale
  # dat$time0 <- with(dat,(time0 - age)*365.25)
  # dat$time1 <- with(dat,(time1 - age)*365.25)
  # 
  # dat_tmp0 <- setDT(subset(dat, onMeds == 0))
  # dat_tmp0$id <- paste0(dat_tmp0$bene_id, "_0")
  # 
  # dat_tmp1 <- setDT(subset(dat, onMeds == 1))
  # dat_meds <- dat_tmp1[dat_tmp1[,.I[time0 == min(time0)], by = bene_id]$V1][!duplicated(bene_id)]
  # dat_tmp1 <- merge(dat_tmp1, data.frame(begin_meds = dat_meds$time0, bene_id = dat_meds$bene_id), by = "bene_id", all.x = TRUE)
  # dat_tmp1$id <- paste0(dat_tmp1$bene_id, "_1")
  # dat_tmp1$time0 <- with(dat_tmp1, time0 - begin_meds)
  # dat_tmp1$time1 <- with(dat_tmp1, time1 - begin_meds)
  # dat_tmp1$begin_meds <- NULL
  # 
  # dat <- rbind(dat_tmp0, dat_tmp1)
  # dat <- setDF(setDT(dat)[order(id, time0)])
  # 
  # rm(dat_tmp0, dat_tmp1, dat_meds);gc()
  
  ## spline cox model
  dat <- setDF(dat[order(dat$bene_id, dat$time0),])
  model_ns <- coxph(Surv(time0, time1, failed) ~ onMeds + pspline(pm_nomed, 4) + pspline(pm_med, 4) + cluster(bene_id),
                    data = dat, id = bene_id, weights = weights.trunc, robust = TRUE, model = TRUE)
  
  print(warnings())
  
  ## save models
  save(model_ns, file=paste0("M:/External Users/KevinJos/output/age_time/cox_hamsm/",fnames[i]))
  
  rm(model_ns, dat)
  gc()
  
}
