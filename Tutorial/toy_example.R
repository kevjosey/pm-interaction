### Interaction Between CVD and PM2.5 ###

library(data.table)
library(readr)
library(survival)
library(splines)
library(ggplot2)
library(ranger)

source("~/Github/pm-steroid/Functions/ipw_fun.R")
source("~/Github/pm-steroid/Tutorial/interaction.R")
dat <- read_csv("~/Github/pm-steroid/Tutorial/simulated_data.csv")
dat <- data.frame(dat)

# Create start/stop times
dat$time0 <- dat$time
dat$time1 <- dat$time + 1
dat$time <- NULL
dat$patient_time <- paste(dat$patient, dat$time1, sep = "-")

### IPW Weights

fmla.denom.pm<-formula('~ sex + age + x1 + x2 + x3 + v1 + v2 + v3') # consider lagged PM

dat$ipw_pm <- ipwtm_ranger(exposure = "PM", denominator = fmla.denom.pm, 
                     id = "patient", timevar = "time1", data = dat, trunc = NULL, 
                     continuous = TRUE, num.trees = 1000, max.depth = 8, num.threads = 4)$ipw.weights
dat$ipw_pm[is.na(dat$ipw_pm)] <- 1

## Steroid IPW Weights

fmla.denom.cvd <- formula('~ PM + sex + age + x1 + x2 + x3 + v1 + v2 + v3') # is CVD montonic (i.e. once CVD = 1 it never changes to 0)?

dat$ipw_cvd <- ipwtm_ranger(exposure = "CVD", denominator = fmla.denom.cvd, 
                            numerator = formula("~ PM"), id = "patient", timevar = "time1", 
                            data = dat, trunc = NULL, continuous = FALSE, num.trees = 1000,
                            max.depth = 8, num.threads = 4)$ipw.weights

dat$ipw_cvd[is.na(dat$cvd)] <- 1

# no censoring by death weights?

### MSM Analysis (with Cox model)

dat$ipw <- unsplit(lapply(split(with(dat, ipw_cvd*ipw_pm), dat$patient), cumprod), dat$patient)

dat$pm_cvd <- ifelse(dat$CVD == 1, dat$PM, 23) 
dat$pm_nocvd <- ifelse(dat$CVD == 0, dat$PM, 23)

## truncate (feel free to try out different levels)
dat$weights.trunc <- dat$ipw
dat$weights.trunc[dat$ipw < quantile(dat$ipw, 0.025)] <- quantile(dat$ipw, 0.025)
dat$weights.trunc[dat$ipw > quantile(dat$ipw, 0.975)] <- quantile(dat$ipw, 0.975)
  
## spline cox model (fails to converge - only 31 events? out of 12000 patients?)
dat <- dat[order(dat$patient, dat$time0),]
model <- coxph(Surv(time0, time1, event) ~ CVD + pspline(pm_nocvd, 4) + pspline(pm_cvd, 4) + cluster(patient),
               data = dat, id = patient, weights = weights.trunc, robust = TRUE, model = TRUE)

# reri components
hr <- data.frame(pm_contrast(model, pm0 = 20, pm1 = 25)) # 5 unit increase in PM
hr$contrast <- paste0(hr$pm1, " vs. ", hr$pm0)

# reri model
reri <- data.frame(t(sapply(seq(20, 26, by = 0.1), function(z, ...) 
  add_interact_cox(model, pm0 = 20, pm1 = z))))
reri$contrast <- paste0(reri$pm1, " vs. ", reri$pm0)
reri[,1:3] <- reri[,1:3]*100

colnames(reri) <- c("est","lower","upper","pm0","pm1","contrast")

g <- ggplot(reri, aes(x = pm1, y = est)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = "dotted")+
  geom_line(size = 1) +
  labs(title = "RERI", x = "PM2.5", y = "Risk Increase due to Interaction (%)")+
  geom_hline(yintercept = 0, color = "blue", linetype = "dashed")+ 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

g
