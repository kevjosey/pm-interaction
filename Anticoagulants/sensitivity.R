library(ggplot2)
library(survival)
library(gridExtra)

filenames <- list.files("M:/External Users/RindalaFayy/Temp/cox model/with spline/", full.names = TRUE)
fnames <- list.files("M:/External Users/RindalaFayy/Temp/cox model/with spline/", full.names = FALSE)

outcomes <- rep(c("Gastrointestinal Bleeding", "Intracranial Bleeding", "Coagulopathy", "All Cause Mortality", "Epistaxis"),2)
meds <- c(rep("Anticoagulant", 5), rep("Antiplatelet", 5))


p <- vector('list', 10)
q <- vector('list', 10)


for (i in 1:length(filenames)){
  
  mod1 <- get(load(filenames[i]))
  
  
  #onMeds = TRUE
  # newdat <- data.frame(onMeds = rep(1,25), pm_med = seq(4,16, by = 0.5), pm_nomed = rep(8,25))
  # HR <- exp(predict(mod1, newdata = newdat))
  # SE <- predict(mod1, newdata = newdat, se.fit = TRUE) 
  # 
  # lower <- exp(log(HR) - 1.96*SE$se.fit)
  # upper <- exp(log(HR) + 1.96*SE$se.fit)
  # 
  # temp_dat <- data.frame(HR = HR, Lower = lower, Upper = upper)
  # 
  # p[[i]] <- ggplot(data = temp_dat, aes(x = seq(4,16, by = 0.5), y = HR, group = 1))+
  #   geom_line(lwd = 1, color = "deeppink4")+
  #   geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.3, fill = "cornflowerblue")+
  #   xlab("PM 2.5 Exposure")+
  #   ylab("Hazard Ratio")+
  #   ggtitle(paste0("Drug = ", meds[i] ,"\nOutcome = ", outcomes[i]))+
  #   theme_light()
  # 
  
  
  
  #onMeds = FALSE
  newdat <- data.frame(onMeds = rep(0,25), pm_med = rep(8,25), pm_nomed = seq(4,16, by = 0.5))
  #HR <- exp(predict(mod1, newdata = newdat))
  log_HR <- predict(mod1, newdata = newdat)
  SE <- predict(mod1, newdata = newdat, se.fit = TRUE) 
  
  #lower <- exp(log(HR) - 1.96*SE$se.fit)
  #upper <- exp(log(HR) + 1.96*SE$se.fit)
  lower <- log_HR - 1.96*SE$se.fit
  upper <- log_HR + 1.96*SE$se.fit
  
  temp_dat <- data.frame(log_HR = log_HR, Lower = lower, Upper = upper)
  
  q[[i]] <- ggplot(data = temp_dat, aes(x = seq(4,16, by = 0.5), y = log_HR, group = 1))+
    geom_line(lwd = 1, color = "deeppink4")+
    geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.3, fill = "cornflowerblue")+
    xlab("PM 2.5 Exposure")+
    ylab("Log Hazard Ratio")+
    ggtitle(paste0("NO ", meds[i], " Exposure \nOutcome = ", outcomes[i]))+
    theme_light()
  
  
}



# pdf("spline_meds_True.pdf", height = 10, width = 20)
# grid.arrange(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], p[[6]], p[[7]], p[[8]], p[[9]], p[[10]], nrow = 2)
# dev.off()


pdf("spline_meds_False.pdf", height = 10, width = 20)
grid.arrange(q[[1]], q[[2]], q[[3]], q[[4]], q[[5]], q[[6]], q[[7]], q[[8]], q[[9]], q[[10]], nrow = 2)
dev.off()