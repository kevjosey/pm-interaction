library(splines)
library(survival)
library(data.table)
library(ggplot2)
library(ggpubr)
library(survminer)
library(sandwich)

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
  hr_out$descript <- rep(c('Corticosteroid Use','Increasing PM while off Corticosteroids','Increasing PM while on Corticosteroids'), 2)
  hr_out$outcome <- outvar[i]
  hr_out$contrast <- paste0(hr_out$pm1, " vs. ", hr_out$pm0)
  hr <- rbind(hr, hr_out)
  
  reri_tmp1 <- data.frame(t(sapply(seq(8, 13, by = 0.1), function(z, ...)
    add_interact_cox(model_ns, pm0 = 8, pm1 = z))))
  reri_tmp2 <- data.frame(t(sapply(seq(5, 13, by = 0.1), function(z, ...)
    add_interact_cox(model_ns, pm0 = 5, pm1 = z))))
  
  reri_out <- data.frame(rbind(reri_tmp1, reri_tmp2))
  
  reri_out$outcome <- outvar[i]
  reri_out$contrast <- paste0(reri_out$pm1, " vs. ", reri_out$pm0)
  reri_out <- reri_out[order(reri_out$pm0, reri_out$pm1),]
  
  colnames(reri_out) <- c("est","lower","upper","pm0","pm1","outcome","contrast")
  reri <- rbind(reri, reri_out)
  
}

reri[,1:3] <- reri[,1:3]*100 # scale to percentile

write.csv(reri, "M:/External Users/KevinJos/output/age_time/sens_msm/reri_cox.csv")
write.csv(coef, "M:/External Users/KevinJos/output/age_time/sens_msm/coef_cox.csv")
write.csv(hr, "M:/External Users/KevinJos/output/age_time/sens_msm/hr_cox.csv")

## Plots

reri <- read.csv("M:/External Users/KevinJos/output/age_time/sens_msm/reri_cox.csv")
hr <- read.csv("M:/External Users/KevinJos/output/age_time/sens_msm/hr_cox.csv")

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

hr_sub <- subset(hr, (pm0 == 5 & pm1 == 10 | pm0 == 8 & pm1 == 12))
hr_sub <- subset(hr_sub, !(descript == "Corticosteroid Use" & pm0 == 5))
hr_sub$descript <- factor(hr_sub$descript, levels = c('Corticosteroid Use', 
                                                      'Increasing PM while off Corticosteroids', 
                                                      'Increasing PM while on Corticosteroids'))
hr_sub$contrast[hr_sub$descript == 'Corticosteroid Use'] <- "N/A"
hr_sub$contrast <- factor(hr_sub$contrast, levels = c("10 vs. 5", "12 vs. 8"))
hr_plot <- merge(hr_sub, cross, by = "outcome")

f <- ggplot(hr_plot, aes(x = nice_names, y = hr, colour = contrast)) +
  facet_grid(descript ~ ., scales = "free") +
  geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.25))+
  labs(title = "Hazard Ratio Estimates", y = "Hazard Ratio", colour = ~ PM[2.5]*" Contrast ("*mu*g*"/"*m^3*")") +
  geom_hline(yintercept = 1, color = "blue", linetype = "dashed") + 
  theme_bw() +
  scale_color_manual(breaks = c("10 vs. 5", "12 vs. 8"),
                     values = c("#F8766D", "#00BFC4")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title.x = element_blank()) +
  grids(linetype = "dashed")

reri$pm0 <- factor(reri$pm0)
cross <- cbind(nice_names = nice_names_2, 
               outcome = c("mi_acs", "iscstroke_tia", "newhf", "newvte", "fib", "death"))
reri_plot <- merge(reri, cross, by = "outcome")
reri_plot$nice_names <- factor(reri_plot$nice_names, levels = nice_names_2)

g <- ggplot(reri_plot, aes(x = pm1, y = est, colour = pm0)) +
  facet_wrap(nice_names ~ ., scales = "free") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = "dotted")+
  geom_line(size = 1) +
  labs(title = , x = ~ PM[2.5]*" ("*mu*g*"/"*m^3*")", y = "Relative Excess Risk Increase due to Interaction (%)", 
       colour = ~"Reference "*PM[2.5]*" Concentration ("*mu*g*"/"*m^3*")")+
  geom_hline(yintercept = 0, color = "blue", linetype = "dashed")+ 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom")

pdf("M:/External Users/KevinJos/output/age_time/sens_msm/hr_plot.pdf", width = 12, height = 8, onefile = FALSE)
f
dev.off()

pdf("M:/External Users/KevinJos/output/age_time/sens_msm/reri_plot.pdf", width = 12, height = 8, onefile = FALSE)
g
dev.off()
