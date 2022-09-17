library(survival)
library(data.table)
library(ggplot2)
library(ggpubr)
library(splines)

source("M:/External Users/KevinJos/code/interaction.R")

outvar<-c('mi_acs','iscstroke_tia','newhf','newvte','fib','death')
aeri1 <- aeri2 <- aeri3 <- aeri4 <- data.frame()

for (i in 1:length(outvar)){
  
  print(i)
  
  ### Main Analysis
  
  load(paste0("M:/External Users/KevinJos/output/age_time/aalen_msm/oralsteroid_",outvar[i],'.RData'))
  load(paste0("M:/External Users/KevinJos/output/fit_data/oralsteroid_",outvar[i],'.RData'))
  
  dat <- setDT(dat)[order(bene_id, time0)]
  mu_pm <- mean(dat$mu_pm)

  idx1 <- which.min(abs(aalen_model$cum[,1] - 70))  
  idx2 <- which.min(abs(aalen_model$cum[,1] - 80))
  idx3 <- which.min(abs(aalen_model$cum[,1] - 90))
  
  idx <- c(idx1, idx2, idx3)
  
  ## over time
  aeri_tmp1 <- data.frame(t(add_interact_aalen(aalen_model, pm0 = 8, pm1 = 12, mu_pm = mu_pm, idx = idx)))
  aeri_tmp2 <- data.frame(t(add_interact_aalen(aalen_model, pm0 = 5, pm1 = 10, mu_pm = mu_pm, idx = idx)))
  
  colnames(aeri_tmp1) <- colnames(aeri_tmp2) <-
    c("time", "pm0", "pm1", "est", "lower", "upper")
  aeri_out1 <- rbind(aeri_tmp1, aeri_tmp2)
  aeri_out1$contrast <- paste0(aeri_out1$pm1, " vs. ", aeri_out1$pm0)
  aeri_out1$age <- rep(c(70, 80, 90), 2)
  aeri_out1$outcome <- outvar[i]
  aeri_out1 <- aeri_out1[order(aeri_out1$pm0, aeri_out1$pm1),]
  
  aeri1 <- rbind(aeri1, aeri_out1)
  
  ## over pm

  aeri_tmp4 <- data.frame(t(sapply(seq(5, 15, by = 0.1), function(z, ...) 
    add_interact_aalen(aalen_model, pm0 = 8, pm1 = z, mu_pm = mu_pm, idx = idx1))))
  aeri_tmp5 <- data.frame(t(sapply(seq(5, 15, by = 0.1), function(z, ...)
    add_interact_aalen(aalen_model, pm0 = 8, pm1 = z, mu_pm = mu_pm, idx = idx2))))
  aeri_tmp6 <- data.frame(t(sapply(seq(5, 15, by = 0.1), function(z, ...)
    add_interact_aalen(aalen_model, pm0 = 8, pm1 = z, mu_pm = mu_pm, idx = idx3))))
  
  aeri_out2 <- data.frame(rbind(aeri_tmp4, aeri_tmp5, aeri_tmp6))
  colnames(aeri_out2) <- c("time", "pm0", "pm1", "est", "lower", "upper")
  aeri_out2$outcome <- outvar[i]
  aeri_out2$age <- rep(c(70, 80, 90), each = nrow(aeri_tmp4))
  aeri_out2 <- aeri_out2[order(aeri_out2$pm0, aeri_out2$pm1),]
  
  aeri2 <- rbind(aeri2, aeri_out2)
  
  ### Sensitivity Analyses
  
  ## all time
  
  max.idx <- which.min(abs(aalen_model$cum[,1] - 95))
  idx_seq <- round(seq(1, max.idx, length.out = max(length(1:max.idx), 1000)))
  
  aeri_tmp7 <- data.frame(t(add_interact_aalen(aalen_model, pm0 = 8, pm1 = 12, mu_pm = mu_pm, idx = idx_seq)))
  aeri_tmp8 <- data.frame(t(add_interact_aalen(aalen_model, pm0 = 5, pm1 = 10, mu_pm = mu_pm, idx = idx_seq)))
  
  colnames(aeri_tmp7) <- colnames(aeri_tmp8) <-
    c("time", "pm0", "pm1", "est", "lower", "upper")
  aeri_out3 <- data.frame(rbind(aeri_tmp7, aeri_tmp8))
  aeri_out3$outcome <- outvar[i]
  aeri_out3$contrast <- paste0(aeri_out3$pm1, " vs. ", aeri_out3$pm0)
  aeri_out3 <- aeri_out3[order(aeri_out3$pm0, aeri_out3$pm1),]
  
  aeri_out3 <- aeri_out3[aeri_out3$lower < aeri_out3$upper,]
  aeri3 <- rbind(aeri3, aeri_out3)
  
  ## montonicity
  
  aeri_tmp0 <- data.frame(t(add_interact_aalen(aalen_model, pm0 = 8, pm1 = 12, mu_pm = mu_pm, idx = idx, monotone = FALSE)))
  
  colnames(aeri_tmp0) <- c("time", "pm0", "pm1", "est", "lower", "upper")
  aeri_out4 <- rbind(aeri_tmp1, aeri_tmp0)
  aeri_out4$age <- rep(c(70, 80, 90), 2)
  aeri_out4$monotone <- rep(c("Monotonicity", "No Monotonicity"), each = 3)
  aeri_out4$outcome <- outvar[i]
  aeri_out4 <- aeri_out4[order(aeri_out4$pm0, aeri_out4$pm1),]
  
  aeri4 <- rbind(aeri4, aeri_out4)
  
}

write.csv(aeri1, "M:/External Users/KevinJos/output/age_time/main/aeri_age.csv")
write.csv(aeri2, "M:/External Users/KevinJos/output/age_time/main/aeri_pm.csv")
write.csv(aeri3, "M:/External Users/KevinJos/output/age_time/main/aeri_curve.csv")
write.csv(aeri4, "M:/External Users/KevinJos/output/age_time/main/aeri_monotone.csv")

## Plots

aeri1 <- read.csv("M:/External Users/KevinJos/output/age_time/main/aeri_age.csv")
aeri2 <- read.csv("M:/External Users/KevinJos/output/age_time/main/aeri_pm.csv")
aeri3 <- read.csv("M:/External Users/KevinJos/output/age_time/main/aeri_curve.csv")
aeri4 <- read.csv("M:/External Users/KevinJos/output/age_time/main/aeri_monotone.csv")

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

aeri1$contrast <- factor(aeri1$contrast, levels = c("10 vs. 5", "12 vs. 8"))
aeri1$age <- factor(aeri1$age, levels = c(70, 80, 90))
plotlist_1 <- list()

for (i in 1:length(outvar)) {
  
  f <- ggplot(subset(aeri1, outcome == outvar[i]), aes(x = age, y = est, colour = contrast)) +
    geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.25))+
    labs(title = nice_names_2[i], x = "Age (Years)", y = "Risk Increase", colour = "PM2.5 Contrast")+
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_hline(yintercept = 0, color = "blue", linetype = "dashed")+ 
    theme_bw()
  
  plotlist_1[[i]] <- f
  
}

aeri2$age <- factor(aeri2$age, levels = c(70, 80, 90))
plotlist_2 <- list()

for (i in 1:length(outvar)) {
  
  g <- ggplot(subset(aeri2, outcome == outvar[i]), aes(x = pm1, y = est, colour = age)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = "dotted")+
    geom_line(size = 1)+
    labs(title = nice_names_2[i], x = "PM2.5", y = "Excess Risk due to Interaction", colour = "Age (Years)")+
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_hline(yintercept = 0, color = "blue", linetype = "dashed")+ 
    theme_bw()
  
  plotlist_2[[i]] <- g
  
}

aeri3$contrast <- factor(aeri3$contrast, levels = c("10 vs. 5", "12 vs. 8"))
aeri3$time <- aeri3$time
plotlist_3 <- list()

for (i in 1:length(outvar)) {
  
  h <- ggplot(subset(aeri3, outcome == outvar[i]), aes(x = time, y = est, colour = contrast)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, linetype = "dotted")+
    geom_line(size = 1)+
    labs(title = nice_names_2[i], x = "Age (Years)", y = "Excess Risk due to Interaction", colour = "PM2.5 Contrast")+
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_hline(yintercept = 0, color = "blue", linetype = "dashed")+ 
    theme_bw()
  
  plotlist_3[[i]] <- h
  
}

aeri4$monotone <- factor(aeri4$monotone, levels = c("Monotonicity", "No Monotonicity"))
aeri4$age <- factor(aeri4$age, levels = c(70, 80, 90))
plotlist_4 <- list()

for (i in 1:length(outvar)) {
  
  a <- ggplot(subset(aeri4, outcome == outvar[i]), aes(x = age, y = est, colour = monotone)) +
    geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.3))+
    labs(title = nice_names_2[i], x = "Age (Years)", y = "Excess Risk due to Interaction", colour = "")+
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_hline(yintercept = 0, color = "blue", linetype = "dashed")+ 
    theme_bw()
  
  plotlist_4[[i]] <- a
  
}

plotlist_5 <- list()

for (i in 1:length(outvar)) {

  load(paste0("M:/External Users/KevinJos/output/age_time/aalen_msm/oralsteroid_",outvar[i],'.RData'))
  load(paste0("M:/External Users/KevinJos/output/fit_data/oralsteroid_",outvar[i],'.RData'))
  
  dat <- setDT(dat)[order(bene_id, time0)]
  
  lp0 <- 1 - exp(-aalen_model$cum[,-1]%*%rbind(c(1,1), c(0,0), c(8,12), c(0,0)))
  lp1 <- 1 - exp(-aalen_model$cum[,-1]%*%rbind(c(1,1), c(1,1), c(8,12), c(8,12)))
  
  colnames(lp0) <- colnames(lp1) <- c(8, 12)
  rownames(lp0) <- rownames(lp1) <- aalen_model$cum[,1]
  tmp <- rbind(reshape2::melt(lp0, value_name = "survival"),
               reshape2::melt(lp1, value_name = "survival"))
  lp <- data.frame(tmp,status = rep(c("Off Steroids", "On Steroids"), each = nrow(tmp)/2))
  colnames(lp) <- c("age", "pm", "survival", "status")
  lp$pm <- factor(lp$pm, levels = c(8,12))
  lp$status <- factor(lp$status, levels = c("Off Steroids", "On Steroids"))
  
  survplot <- ggplot(lp, aes(x = age, y = survival, colour = pm, linetype = status)) +
    geom_line(size = 1.2) +
    labs(title = nice_names_2[i], x = "Age (Years)", y = "Event Probability", 
         colour = "PM2.5", linetype = "Medication Status")+
    theme(plot.title = element_text(hjust = 0.5)) +
    theme_bw()
  
  plotlist_5[[i]] <- survplot
  
}

pdf("M:/External Users/KevinJos/output/age_time/main/aeri_time_plot.pdf", width = 12, height = 8, onefile = FALSE)
ggarrange(plotlist_1[[1]], plotlist_1[[2]], plotlist_1[[3]],
          plotlist_1[[4]], plotlist_1[[5]], plotlist_1[[6]],
          ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")
dev.off()

pdf("M:/External Users/KevinJos/output/age_time/main/aeri_pm_plot.pdf", width = 12, height = 8, onefile = FALSE)
ggarrange(plotlist_2[[1]], plotlist_2[[2]], plotlist_2[[3]],
          plotlist_2[[4]], plotlist_2[[5]], plotlist_2[[6]],
          ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")
dev.off()

pdf("M:/External Users/KevinJos/output/age_time/main/aeri_bigtime_plot.pdf", width = 12, height = 8, onefile = FALSE)
ggarrange(plotlist_3[[1]], plotlist_3[[2]], plotlist_3[[3]],
          plotlist_3[[4]], plotlist_3[[5]], plotlist_3[[6]],
          ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")
dev.off()

pdf("M:/External Users/KevinJos/output/age_time/main/aeri_monotone_plot.pdf", width = 12, height = 8, onefile = FALSE)
ggarrange(plotlist_4[[1]], plotlist_4[[2]], plotlist_4[[3]],
          plotlist_4[[4]], plotlist_4[[5]], plotlist_4[[6]],
          ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")
dev.off()

pdf("M:/External Users/KevinJos/output/age_time/main/aalen_survival_plot.pdf", width = 12, height = 8, onefile = FALSE)
ggarrange(plotlist_5[[1]], plotlist_5[[2]], plotlist_5[[3]],
          plotlist_5[[4]], plotlist_5[[5]], plotlist_5[[6]],
          ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")
dev.off()
