library(splines)
library(survival)
library(data.table)
library(ggplot2)
library(ggpubr)
library(survminer)
library(sandwich)

source("M:/External Users/KevinJos/code/interaction.R")

outvar <- c("bleeding_allgi", "bleeding_intracranial", "epistaxis", "death") 
reri <- coef <- hr <- data.frame()

for (i in 1:length(outvar)){
  
  print(i)
  
  load(paste0("M:/External Users/KevinJos/output/age_time/cox/oralantic_",outvar[i],'.RData'))
  load(paste0("M:/External Users/KevinJos/output/ipw/oralantic_",outvar[i],'.RData'))
  
  # outcome coefficients 
  coef_val <- cbind(outcome = outvar[i], summary(cox_model)$coefficients)
  coef <- rbind(coef, coef_val)
  
  # reri components
  hr_tmp1 <- pm_contrast(cox_model, pm0 = 5, pm1 = 10)
  hr_tmp2 <- pm_contrast(cox_model, pm0 = 8, pm1 = 12)
  
  hr_out <- data.frame(rbind(hr_tmp1, hr_tmp2))
  hr_out$descript <- rep(c('Anticoagulant Use','Increasing PM while off Anticoagulants','Increasing PM while on Anticoagulants'), 2)
  hr_out$outcome <- outvar[i]
  hr_out$contrast <- paste0(hr_out$pm1, " vs. ", hr_out$pm0)
  hr <- rbind(hr, hr_out)
  
  reri_tmp1 <- data.frame(t(sapply(seq(8, 13, by = 0.1), function(z, ...)
    add_interact_cox(cox_model, pm0 = 8, pm1 = z))))
  reri_tmp2 <- data.frame(t(sapply(seq(5, 13, by = 0.1), function(z, ...)
    add_interact_cox(cox_model, pm0 = 5, pm1 = z))))
  
  reri_out <- data.frame(rbind(reri_tmp1, reri_tmp2))
  
  reri_out$outcome <- outvar[i]
  reri_out$contrast <- paste0(reri_out$pm1, " vs. ", reri_out$pm0)
  reri_out <- reri_out[order(reri_out$pm0, reri_out$pm1),]
  
  colnames(reri_out) <- c("est","lower","upper","pm0","pm1","outcome","contrast")
  reri <- rbind(reri, reri_out)
  
}

reri[,1:3] <- reri[,1:3]*100 # scale to percentile

write.csv(reri, "M:/External Users/KevinJos/output/age_time/anticoagulants/reri_cox.csv")
write.csv(coef, "M:/External Users/KevinJos/output/age_time/anticoagulants/coef_cox.csv")
write.csv(hr, "M:/External Users/KevinJos/output/age_time/anticoagulants/hr_cox.csv")

## Plots

reri <- read.csv("M:/External Users/KevinJos/output/age_time/anticoagulants/reri_cox.csv")
hr <- read.csv("M:/External Users/KevinJos/output/age_time/anticoagulants/hr_cox.csv")

nice_names_1<-c('Gastrointestinal\n Bleeding',
                'Intracanial\n Bleeding',
                'Epistaxis',
                'Death')

nice_names_2<-c('Gastrointestinal Bleeding',
                'Intracanial Bleeding',
                'Epistaxis',
                'Death')

cross <- cbind(nice_names = nice_names_1, 
               outcome = c("bleeding_allgi", "bleeding_intracranial", "epistaxis", "death"))

hr_sub <- subset(hr, (pm0 == 5 & pm1 == 10 | pm0 == 8 & pm1 == 12))
hr_sub <- subset(hr_sub, !(descript == "Anticoagulant Use" & pm0 == 5))
hr_sub$descript <- factor(hr_sub$descript, levels = c('Anticoagulant Use', 
                                                      'Increasing PM while off Anticoagulants', 
                                                      'Increasing PM while on Anticoagulants'))
hr_sub$contrast[hr_sub$descript == 'Anticoagulant Use'] <- "N/A"
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
               outcome = c("bleeding_allgi", "bleeding_intracranial", "epistaxis", "death"))
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

# plotlist_2 <- list()
# 
# for (i in 1:length(outvar)) {
#   
#   print(i)
# 
#   load(paste0("M:/External Users/KevinJos/output/age_time/steroids/cox/oralsteroid_",outvar[i],'.RData'))
#   
#   newdata0 <- data.frame(onMeds = c(0,0), pm_nomed = c(8,12), pm_med = c(8,8), weights.trunc = c(1,1))
#   newdata1 <- data.frame(onMeds = c(1,1), pm_nomed = c(8,8), pm_med = c(8,12), weights.trunc = c(1,1))
# 
#   s.obj.0 <- survfit(cox_model, newdata = newdata0, stype = 1)
#   s.obj.1 <- survfit(cox_model, newdata = newdata1, stype = 1)
#   lp0 <- 1 - s.obj.0$surv
#   lp1 <- 1 - s.obj.1$surv
#   colnames(lp0) <- colnames(lp1) <- c('8', '12')
#   rownames(lp0) <- s.obj.0$time
#   rownames(lp1) <- s.obj.1$time
# 
#   tmp <- rbind(reshape2::melt(lp0, value_name = "survival"),
#                reshape2::melt(lp1, value_name = "survival"))
#   lp <- data.frame(tmp, status = rep(c("Off Steroids", "On Steroids"), each = nrow(tmp)/2))
#   colnames(lp) <- c("age", "pm", "survival", "status")
#   lp$pm <- factor(lp$pm, levels = c(8,12))
#   lp$status <- factor(lp$status, levels = c("Off Steroids", "On Steroids"))
#   lp <- lp[lp$age <= 95,]
# 
#   survplot <- ggplot(lp, aes(x = age, y = survival, colour = pm, linetype = status)) +
#     geom_line(size = 1.2) +
#     labs(title = nice_names_2[i], x = "Age (Years)", y = "Event Probability", 
#          colour = ~ PM[2.5]*" ("*mu*g*"/"*m^3*")", linetype = "Medication Status") +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5, face = "bold"))
#   
#   plotlist_2[[i]] <- survplot
# 
# }

pdf("M:/External Users/KevinJos/output/age_time/anticoagulants/hr_plot.pdf", width = 12, height = 8, onefile = FALSE)
f
dev.off()

pdf("M:/External Users/KevinJos/output/age_time/anticoagulants/reri_plot.pdf", width = 12, height = 8, onefile = FALSE)
g
dev.off()

# pdf("M:/External Users/KevinJos/output/age_time/steroids/cox_survival_plot.pdf", width = 12, height = 8, onefile = FALSE)
# ggarrange(plotlist_2[[1]], plotlist_2[[2]], plotlist_2[[3]],
#           plotlist_2[[4]], plotlist_2[[5]], plotlist_2[[6]],
#           ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")
# dev.off()