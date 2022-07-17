library(splines)
library(survival)
library(data.table)
library(ggplot2)
library(ggpubr)
library(survminer)

source("M:/External Users/KevinJos/code/interaction.R")

outvar<-c('mi_acs','iscstroke_tia','newhf','newvte','fib','death')
reri <- coef <- hr <- data.frame()

for (i in 1:length(outvar)){
  
  print(i)
  
  load(paste0("M:/External Users/KevinJos/output/age_time/cox_msm/oralsteroid_",outvar[i],'.RData'))
  load(paste0("M:/External Users/KevinJos/output/fit_data/oralsteroid_",outvar[i],'.RData'))

  # outcome coefficients 
  coef_val <- cbind(outcome = outvar[i], summary(model_ns)$coefficients)
  coef <- rbind(coef, coef_val)
  
  # reri components
  hr_tmp1 <- pm_contrast(model_ns, pm0 = 5, pm1 = 10)
  hr_tmp2 <- pm_contrast(model_ns, pm0 = 8, pm1 = 12)
  
  hr_out <- data.frame(rbind(hr_tmp1, hr_tmp2))
  hr_out$descript <- rep(c("-", "Off Corticteroids", "On Corticosteroids"), 2)
  hr_out$outcome <- outvar[i]
  hr_out$contrast <- paste0(hr_out$pm1, " vs. ", hr_out$pm0)
  hr <- rbind(hr, hr_out)
  
  # reri model
  reri_tmp1 <- data.frame(t(sapply(seq(8, 15, by = 0.1), function(z, ...)
    add_interact_cox(model_ns, pm0 = 8, pm1 = z))))
  reri_tmp2 <- data.frame(t(sapply(seq(5, 15, by = 0.1), function(z, ...)
    add_interact_cox(model_ns, pm0 = 5, pm1 = z))))
  
  reri_out <- data.frame(rbind(reri_tmp1, reri_tmp2))
  
  reri_out$outcome <- outvar[i]
  reri_out$contrast <- paste0(reri_out$pm1, " vs. ", reri_out$pm0)
  reri_out <- reri_out[order(reri_out$pm0, reri_out$pm1),]
  
  colnames(reri_out) <- c("est","lower","upper","pm0","pm1","outcome","contrast")
  reri <- rbind(reri, reri_out)
  
}

reri[,1:3] <- reri[,1:3]*100 # scale to percentile

write.csv(reri, "M:/External Users/KevinJos/output/age_time/main/reri_cox.csv")
write.csv(coef, "M:/External Users/KevinJos/output/age_time/main/coef_cox.csv")
write.csv(hr, "M:/External Users/KevinJos/output/age_time/main/hr_cox.csv")

## Plots

reri <- read.csv("M:/External Users/KevinJos/output/age_time/main/reri_cox.csv")
reri <- read.csv("M:/External Users/KevinJos/output/age_time/main/hr_cox.csv")

nice_names_1<-c('Myocardial Infarction or\nAcute Coronary Syndrome',
                'Ischemic Stroke or\nTransient Ischemic Attack',
                'Heart Failure',
                'Venous Thromboembolism',
                'Atrial Fibrillation',
                'All-Cause Mortality')

nice_names_2<-c('Myocardial Infarction or ACS',
                'Ischemic Stroke or TIA',
                'Heart Failure',
                'Venous Thromboembolism',
                'Atrial Fibrillation',
                'All-Cause Mortality')

cross <- cbind(nice_names = nice_names_1, 
               outcome = c("mi_acs", "iscstroke_tia", "newhf", "newvte", "fib", "death"))

hr_sub <- subset(hr, p(m0 == 5 & pm1 == 10 | pm0 == 8 & pm1 == 12) & descript != "-")
hr_sub$descript <- factor(hr_sub$descript, levels = c("Off Corticosteroids", "On Corticosteroids"))
hr_sub$contrast <- factor(hr_sub$contrast, levels = c("10 vs. 5", "12 vs. 8"))
hr_plot <- merge(hr_sub, cross, by = "outcome")

f <- ggplot(hr_plot, aes(x = nice_names, y = est, colour = contrast, shape = descript)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.25))+
  labs(title = "Hazard Ratio Estimates from Increasing PM2.5", x = "Outcome", 
       y = "Hazard Ratio", colour = "PM2.5 Contrast", shape = "Medication Status") +
  geom_hline(yintercept = 0, color = "blue", linetype = "dashed")+ 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

reri$pm0 <- factor(reri$pm0)

plotlist_1 <- list()

for (i in 1:length(outvar)) {
  
  g <- ggplot(subset(reri, outcome == outvar[i]), aes(x = pm1, y = est, colour = pm0)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = "dotted")+
    geom_line(size = 1) +
    labs(title = nice_names_2[i], x = "PM2.5", y = "Relative Excess Risk Increase due to Interaction (%)", colour = "Reference PM2.5")+
    geom_hline(yintercept = 0, color = "blue", linetype = "dashed")+ 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  
  plotlist_1[[i]] <- g
  
}

plotlist_2 <- list()

for (i in 1:length(outvar)) {
  
  print(i)

  load(paste0("M:/External Users/KevinJos/output/age_time/cox_msm/oralsteroid_",outvar[i],'.RData'))
  
  newdata0 <- data.frame(onMeds = c(0,0), pm_nomed = c(8,12), pm_med = c(8,8), weights.trunc = c(1,1))
  newdata1 <- data.frame(onMeds = c(1,1), pm_nomed = c(8,8), pm_med = c(8,12), weights.trunc = c(1,1))

  s.obj.0 <- survfit(model_ns, newdata = newdata0, stype = 1)
  s.obj.1 <- survfit(model_ns, newdata = newdata1, stype = 1)
  lp0 <- 1 - s.obj.0$surv
  lp1 <- 1 - s.obj.1$surv
  colnames(lp0) <- colnames(lp1) <- c('8', '12')
  rownames(lp0) <- s.obj.0$time
  rownames(lp1) <- s.obj.1$time

  tmp <- rbind(reshape2::melt(lp0, value_name = "survival"),
               reshape2::melt(lp1, value_name = "survival"))
  lp <- data.frame(tmp, status = rep(c("Off Steroids", "On Steroids"), each = nrow(tmp)/2))
  colnames(lp) <- c("age", "pm", "survival", "status")
  lp$pm <- factor(lp$pm, levels = c(8,12))
  lp$status <- factor(lp$status, levels = c("Off Steroids", "On Steroids"))
  lp <- lp[lp$age <= 95,]

  survplot <- ggplot(lp, aes(x = age, y = survival, colour = pm, linetype = status)) +
    geom_line(size = 1.2) +
    labs(title = nice_names_2[i], x = "Age (Years)", y = "Event Probability", 
         colour = "PM2.5", linetype = "Medication Status") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  
  plotlist_2[[i]] <- survplot

}

pdf("M:/External Users/KevinJos/output/age_time/main/hr_plot.pdf", width = 12, height = 8, onefile = FALSE)
f
dev.off()

pdf("M:/External Users/KevinJos/output/age_time/main/reri_plot.pdf", width = 12, height = 8, onefile = FALSE)
ggarrange(plotlist_1[[1]], plotlist_1[[2]], plotlist_1[[3]],
          plotlist_1[[4]], plotlist_1[[5]], plotlist_1[[6]],
          ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")
dev.off()

pdf("M:/External Users/KevinJos/output/age_time/main/cox_survival_plot.pdf", width = 12, height = 8, onefile = FALSE)
ggarrange(plotlist_2[[1]], plotlist_2[[2]], plotlist_2[[3]],
          plotlist_2[[4]], plotlist_2[[5]], plotlist_2[[6]],
          ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")
dev.off()