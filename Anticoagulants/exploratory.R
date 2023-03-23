
library(haven)
library(data.table)
library(lubridate)
library(stringr)
library(ranger)
library(geepack)
library(survival)
library(plyr)
library(table1)
library(usmap)
library(maps)
library(ggplot2)
library(maps)
library(tidyverse)
library(flextable)
library(magrittr)
library(kableExtra)

cohort <- read_sas("M:/External Users/Shared/To RachelNet/cohortfinalsept2022_a2c2.sas7bdat")
backup <- cohort

#get data to draw map of US states
us_states <- map_data("state")

#read in data that gives Lat, Long and state data for each ZIP code
latlong <- read.csv("M:/External Users/RindalaFayy/Transfer/In/us_zips_states.csv")

#making all the state names lower case
latlong$state_name <- tolower(latlong$state_name)

#make all the zip codes 5 digitis long by adding zeros to the left
latlong$zip <- str_pad(latlong$zip, width = 5, side = "left", pad = "0")

#selecting only the relevant columns
latlong <- latlong %>% dplyr::select(zip, lat, lng, city, state_name, population, density)

cohort <- cohort %>% dplyr::left_join(latlong, by = "zip")

#only selecting relevant variables
cohort <- cohort %>% dplyr::select(BENE_ID, zip, dob, deathdate, sex, race, dtstart, dtend, enrol_mths, indexdate, first_date_oralantip, first_date_oralantic, last_date_oralantip, last_date_oralantic, death, age,lat, lng, city, state_name, population, density, deathdate, first_bleeding_allGI_date, First_bleeding_intracranial_date, First_coagulopathy_date, First_Epistaxis_date)



# 2010 Anticoagulants
prop_2010_antic <- cohort %>% dplyr::group_by(state_name) %>%
  dplyr::summarize(prop_antic_2010 = sum(year(ymd(first_date_oralantic)) == 2010, na.rm = TRUE)/
                     
                     
                     
                     sum(year(ymd(indexdate)) <= 2010 & 
                           
                           year(pmin(ymd(dtend), ymd(deathdate), ymd(first_bleeding_allGI_date), ymd(First_bleeding_intracranial_date), ymd(First_coagulopathy_date), ymd(First_Epistaxis_date), ymd(last_date_oralantic), na.rm = TRUE)) >= 2010 
                         
                         & (year(ymd(first_date_oralantic)) >= 2010 | is.na(first_date_oralantic))))





us_states <- us_states %>% dplyr::left_join(prop_2010_antic, by = c("region" = "state_name"))

us_states <- us_states %>%
  mutate(manual_fill_antic = cut(prop_antic_2010, breaks = c(0.00, 0.03, 0.06, 0.09, 0.12, 0.16), 
                                 labels = c("0.0%-2.9%", "3.0%-5.9%", "6.0%-8.9%", "9.0%-11.9%" , "12.0% - 15%"), 
                                 right = TRUE))


pal <- c("cornsilk", "antiquewhite", "bisque2", "bisque3", "burlywood3", "burlywood4")


prop_2010_map_antic <- us_states %>% 
  ggplot()+
  geom_polygon(aes(x = long, y = lat, group = group, fill = manual_fill_antic), color = "black")+
  scale_fill_manual(name = "Proportion", values = pal)+
  coord_fixed(1.5)+
  theme(panel.grid.major = element_blank(),
        panel.background =  element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  ggtitle(str_wrap("Proportion of population at risk who received Anticoagulants for the first time in 2010", 60))


# 2014 Anticoagulants
prop_2014_antic <- cohort %>% dplyr::group_by(state_name) %>%
  dplyr::summarize(prop_antic_2014 = sum(year(ymd(first_date_oralantic)) == 2014, na.rm = TRUE)/sum(year(ymd(indexdate)) <= 2014 & year(pmin(ymd(dtend), ymd(deathdate), ymd(first_bleeding_allGI_date), ymd(First_bleeding_intracranial_date), ymd(First_coagulopathy_date), ymd(First_Epistaxis_date), ymd(last_date_oralantic), na.rm = TRUE)) >= 2014 & (year(ymd(first_date_oralantic)) >= 2014 | is.na(first_date_oralantic))))


us_states <- us_states %>% dplyr::left_join(prop_2014_antic, by = c("region" = "state_name"))

us_states <- us_states %>%
  mutate(manual_fill_antic2 = cut(prop_antic_2014, breaks = c(0.00, 0.03, 0.06, 0.09, 0.12, 0.16), 
                                  labels = c("0.0%-2.9%", "3.0%-5.9%", "6.0%-8.9%", "9.0%-11.9%" , "12.0% - 15%"), 
                                  right = TRUE))


pal <- c("antiquewhite", "bisque2", "bisque3", "burlywood3", "burlywood4")


prop_2014_map_antic <- us_states %>% 
  ggplot()+
  geom_polygon(aes(x = long, y = lat, group = group, fill = manual_fill_antic2), color = "black")+
  scale_fill_manual(name = "Proportion", values = pal)+
  coord_fixed(1.5)+
  theme(panel.grid.major = element_blank(),
        panel.background =  element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  ggtitle(str_wrap("Proportion of population at risk who received Anticoagulants for the first time in 2014", 60))



# 2010 Antiplatelets
prop_2010_antip <- cohort %>% dplyr::group_by(state_name) %>%
  dplyr::summarize(prop_antip_2010 = sum(year(ymd(first_date_oralantip)) == 2010, na.rm = TRUE)/sum(year(ymd(indexdate)) <= 2010 & year(pmin(ymd(dtend), ymd(deathdate), ymd(first_bleeding_allGI_date), ymd(First_bleeding_intracranial_date), ymd(First_coagulopathy_date), ymd(First_Epistaxis_date), ymd(last_date_oralantip), na.rm = TRUE)) >= 2010 & (year(ymd(first_date_oralantip)) >= 2010 | is.na(first_date_oralantip))))


us_states <- us_states %>% dplyr::left_join(prop_2010_antip, by = c("region" = "state_name"))

us_states <- us_states %>%
  mutate(manual_fill_antip = cut(prop_antip_2010, breaks = c(0.00, 0.03, 0.06, 0.09, 0.12, 0.16), 
                                 labels = c("0.0%-2.9%", "3.0%-5.9%", "6.0%-8.9%", "9.0%-11.9%" , "12.0% - 15%"), 
                                 right = TRUE))


pal <- c("cornsilk", "antiquewhite", "bisque2", "bisque3", "burlywood3", "burlywood4")


prop_2010_map_antip <- us_states %>% 
  ggplot()+
  geom_polygon(aes(x = long, y = lat, group = group, fill = manual_fill_antip), color = "black")+
  scale_fill_manual(name = "Proportion", values = pal)+
  coord_fixed(1.5)+
  theme(panel.grid.major = element_blank(),
        panel.background =  element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  ggtitle(str_wrap("Proportion of population at risk who received Antiplatelets for the first time in 2010", 60))


# 2014 Antiplatelets
prop_2014_antip <- cohort %>% dplyr::group_by(state_name) %>%
  dplyr::summarize(prop_antip_2014 = sum(year(ymd(first_date_oralantip)) == 2014, na.rm = TRUE)/sum(year(ymd(indexdate)) <= 2014 & year(pmin(ymd(dtend), ymd(deathdate), ymd(first_bleeding_allGI_date), ymd(First_bleeding_intracranial_date), ymd(First_coagulopathy_date), ymd(First_Epistaxis_date), ymd(last_date_oralantip), na.rm = TRUE)) >= 2014 & (year(ymd(first_date_oralantip)) >= 2014 | is.na(first_date_oralantip))))


us_states <- us_states %>% dplyr::left_join(prop_2014_antip, by = c("region" = "state_name"))

us_states <- us_states %>%
  mutate(manual_fill_antip2 = cut(prop_antip_2014, breaks = c(0.00, 0.03, 0.06, 0.09, 0.12, 0.16), 
                                  labels = c("0.0%-2.9%", "3.0%-5.9%", "6.0%-8.9%", "9.0%-11.9%" , "12.0% - 15%"), 
                                  right = TRUE))


pal <- c("cornsilk", "antiquewhite", "bisque2", "bisque3", "burlywood3", "burlywood4")


prop_2014_map_antip <- us_states %>% 
  ggplot()+
  geom_polygon(aes(x = long, y = lat, group = group, fill = manual_fill_antip2), color = "black")+
  scale_fill_manual(name = "Proportion", values = pal)+
  coord_fixed(1.5)+
  theme(panel.grid.major = element_blank(),
        panel.background =  element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  ggtitle(str_wrap("Proportion of population at risk who received Antiplatelets for the first time in 2014", 60))


prop_2010_map_antic
prop_2014_map_antic
prop_2010_map_antip
prop_2014_map_antip

pdf("prop_antip_2010v2014.pdf", height = 10, width = 8)
par(mfrow = c(2,1))
prop_2010_map_antip
prop_2014_map_antip
dev.off()

pdf("prop_antic_2010v2014.pdf", height = 10, width = 8)
par(mfrow = c(2,1))
prop_2010_map_antic
prop_2014_map_antic
dev.off()



