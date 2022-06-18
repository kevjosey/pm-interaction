## Make table 1 ##

library(survival)
library(dplyr)
library(haven)
library(data.table)
library(lubridate)
library(stringr)
library(ipw)
library(zoo)

ndig <- 2 # significant figures

## read in the full dataset ##
cohort <- read_sas("M:/External Users/Shared/To RachelNet/cohortfinalaug2021_a2c2.sas7bdat")
names(cohort)<-tolower(names(cohort))
setDT(cohort)
medclass <- "oralsteroid"

cohort <- cohort[,enddate:=pmin(deathdate,dtend,ymd("2016-11-30"),
                      get(paste0("last_date_", medclass)), na.rm=T)
       ][,firstmed:=get(paste0('first_date_',medclass))
        ][,lastmed:=get(paste0('last_date_',medclass))
          # and remove people who start taking meds or experience a terminating event before their index date
          ][indexdate<=pmin(firstmed,enddate,na.rm=T)
            # if a person first takes meds after a terminating event, treat as though they never take meds
            ][!is.na(firstmed) & firstmed>enddate,':='(firstmed=NA,lastmed=NA)]

## remove sex=0 (there's only one person with this) ##
cohort <- cohort[!(sex=='0')]

enddate <- cohort$enddate
firstmed <- cohort$firstmed
lastmed <- cohort$lastmed
onMeds <- firstmed <= enddate
onMeds[is.na(onMeds)] <- FALSE

## write a function to add a row for continuous variables

add_cont <- function(x, nm_var, nm_level, onMeds=NULL){
  
  if (!is.null(onMeds)){
    n_x<-c(sum(!is.na(x[!onMeds])), sum(!is.na(x[onMeds])))
    mean_x<-c(round(mean(x[!onMeds], na.rm=T), ndig), round(mean(x[onMeds], na.rm=T), ndig))
    sd_x<-c(round(sd(x[!onMeds], na.rm=T), ndig), round(sd(x[onMeds], na.rm=T), ndig))
    return(c(n_x[1], mean_x[1], sd_x[1], n_x[2], mean_x[2], sd_x[2]))
  } else {
    n_x<-sum(!is.na(x))
    mean_x<-round(mean(x, na.rm=T),ndig)
    sd_x<-round(sd(x, na.rm=T),ndig)
    return(c(nm_var, nm_level, n_x, mean_x, sd_x))
  }

}

## write a function to add a row for a level of a categorical variable
add_cat <- function(x,nm_var,nm_level,onMeds=NULL){
  
  if (!is.null(onMeds)){
    n_x<-c(sum(!is.na(x[!onMeds])), sum(!is.na(x[onMeds])))
    n_x1<-c(sum(x[!onMeds] == 1, na.rm=T), sum(x[onMeds]==1,na.rm=T))
    pct_x1<-c(round(100*n_x1[1]/n_x[1], ndig), round(100*n_x1[2]/n_x[2], ndig))
    return(c(n_x[1], n_x1[1], pct_x1[1], n_x[2], n_x1[2], pct_x1[2]))
  } else {
    n_x<-sum(!is.na(x))
    n_x1<-sum(x==1, na.rm=T)
    pct_x1<-round(100*n_x1/n_x, ndig)
    return(c(nm_var, nm_level, n_x, n_x1, pct_x1))
  }
  
}

# Table 1 -----------------------------------------------------------------

## age ##
table1<-c(add_cont(x=cohort[,age],nm_var='Age at Index',nm_level=''),
          add_cont(x=cohort[,age],nm_var='Age at Index',nm_level='',onMeds=onMeds))

## sex ##
table1<-rbind(table1,c(add_cat(x=as.numeric(cohort[,sex]==1),nm_var='Male',nm_level=''),
                       add_cat(x=as.numeric(cohort[,sex]==1),nm_var='Male',nm_level='',onMeds=onMeds)))

## race ##
table1<-rbind(table1,c(add_cat(x=as.numeric(cohort[,race]==1),nm_var='Race',nm_level='White'),
                       add_cat(x=as.numeric(cohort[,race]==1),nm_var='Race',nm_level='White',onMeds=onMeds)))
table1<-rbind(table1,c(add_cat(x=as.numeric(cohort[,race]==2),nm_var='Race',nm_level='Black'),
                       add_cat(x=as.numeric(cohort[,race]==2),nm_var='Race',nm_level='Black',onMeds=onMeds)))
table1<-rbind(table1,c(add_cat(x=as.numeric(cohort[,race]==5),nm_var='Race',nm_level='Hispanic'),
                       add_cat(x=as.numeric(cohort[,race]==5),nm_var='Race',nm_level='Hispanic',onMeds=onMeds)))
table1<-rbind(table1,c(add_cat(x=as.numeric(cohort[,race]==4),nm_var='Race',nm_level='Asian'),
                       add_cat(x=as.numeric(cohort[,race]==4),nm_var='Race',nm_level='Asian',onMeds=onMeds)))
table1<-rbind(table1,c(add_cat(x=as.numeric(cohort[,race]==6),nm_var='Race',nm_level='North American Native'),
                       add_cat(x=as.numeric(cohort[,race]==6),nm_var='Race',nm_level='North American Native',onMeds=onMeds)))
table1<-rbind(table1,c(add_cat(x=as.numeric(cohort[,race]==3),nm_var='Race',nm_level='Other'),
                       add_cat(x=as.numeric(cohort[,race]==3),nm_var='Race',nm_level='Other',onMeds=onMeds)))

## dual eligible ##
table1<-rbind(table1,c(add_cat(x=as.numeric(cohort[,dualeligible]==1),nm_var='Medicaid Eligible',nm_level=''),
                       add_cat(x=as.numeric(cohort[,dualeligible]==1),nm_var='Medicaid Eligible',nm_level='',onMeds=onMeds)))

## chronic conditions ##
cc<-c("chronic_tha", "chronic_tka","chronic_acs","chronic_cancer","chronic_fib","chronic_hemstroke","chronic_hf",       
      "chronic_iscstroke","chronic_mi","chronic_pvd","chronic_tia","chronic_vte","chronic_mi_acs","chronic_cva","chronic_tja", "chronic_carotid")
# comorbidities and medication/health hx
dx <- names(cohort)[grep('dx_',names(cohort))]
hx <- c("hx_n_hospvisits", "hx_n_ervisits", "hx_n_outpvisits", "hx_n_meds")

for (i in 1:length(cc)){
  table1<-rbind(table1,c(add_cat(x=as.numeric(cohort[,get(cc[i])]==1),nm_var='Chronic Condition',nm_level=cc[i]),
                         add_cat(x=as.numeric(cohort[,get(cc[i])]==1),nm_var='Chronic Condition',nm_level=cc[i],onMeds=onMeds)))
}

for (i in 1:length(dx)){
  table1<-rbind(table1,c(add_cat(x=as.numeric(cohort[,get(dx[i])]==1),nm_var='Comorbidities',nm_level=dx[i]),
                         add_cat(x=as.numeric(cohort[,get(dx[i])]==1),nm_var='Comorbidities',nm_level=dx[i],onMeds=onMeds)))
}

for (i in 1:length(hx)){
  table1<-rbind(table1,c(add_cont(x=as.numeric(cohort[,get(hx[i])]),nm_var='Hospitalization History',nm_level=hx[i]),
                         add_cont(x=as.numeric(cohort[,get(hx[i])]),nm_var='Hospitalization History',nm_level=hx[i],onMeds=onMeds)))
}

table1 <- as.data.frame(table1)
names(table1)<-c('Variable','Level', 'Size - All', 'Mean (or N) - All', 'SD (or %) - All', 
                 'Size - No Meds','Mean (or N) - No Meds','SD (or %) - No Meds',
                 'Size - Meds', 'Mean (or N) Meds','SD (or %) Meds')
write.csv(table1, "M:/External Users/KevinJos/output/table1.csv")

# Table 2 -----------------------------------------------------------------

## compute average PM2.5 across all zipcodes used in the study for all years of study ##
## first bring in moving dataset to get extra zipcodes ##
changezip <- read_sas("M:/External Users/Shared/To RachelNet/morethanonezipcode.sas7bdat")
names(changezip)<-tolower(names(changezip))
setDT(changezip)

## make a vector of all zips in the study ##
all_zips<-c(cohort[,zip],changezip[,zip])

## now merge in pm2.5 data ##
pm <- fread("M:/External Users/RachelNet/data/pollution/pm25_seasonalavg_zipcode.csv")
pm[,zip:=str_pad(ZIP,width = 5,side='left',pad='0')]
pm[,yr_in_study:=as.numeric(year>=2008 & year<=2016)]
pm[,zip_in_study:=as.numeric(zip %in% all_zips)]
pm$pm25 <- rowMeans(pm[,c("pm25_winter","pm25_spring","pm25_summer","pm25_fall")])

pm_sub<-pm[yr_in_study == 1 & zip_in_study == 1]

table2<-add_cont(x = pm_sub[,pm25], nm_var='PM2.5', nm_level = 'Combined')
table2<-rbind(table2, add_cont(x=pm_sub[,pm25_winter],nm_var='PM2.5',nm_level='Winter'))
table2<-rbind(table2, add_cont(x=pm_sub[,pm25_spring],nm_var='PM2.5',nm_level='Spring'))
table2<-rbind(table2, add_cont(x=pm_sub[,pm25_summer],nm_var='PM2.5',nm_level='Summer'))
table2<-rbind(table2, add_cont(x=pm_sub[,pm25_fall],nm_var='PM2.5',nm_level='Fall'))

## use same approach to add the neighborhood level features ##
## read in and clean ##
conf <- fread("M:/External Users/RachelNet/data/confounders/census_interpolated_zips.csv")
## subset to only 2008 and later (dates when medicare data available) ##
conf<-conf[year>=2008
           ## make zipcode into a 5-digit character string ##
           ][,zip:=str_pad(ZIP,width=5,side='left',pad='0')
             ## remove unnecessary columns
             ][,c('V1','zcta','ZIP'):=NULL]
conf[,zip_in_study:=as.numeric(zip %in% all_zips)]

conf_sub<-conf[zip_in_study==1]

neigh<-c('popdensity','medianhousevalue','medhouseholdincome','poverty','pct_owner_occ','hispanic','pct_blk','pct_white','education')
nice_names<-c('Population Density','Median House Value','Median Household Income','% Poverty','% Owner Occupied Housing',
              '% Hispanic','% Black', '% White', '% with Bachelors Degree or Higher')

## population density can be used as-is ##
table2<-rbind(table2, add_cont(x=conf_sub[,get(neigh[1])],nm_var=nice_names[1],nm_level='')) # density
table2<-rbind(table2, add_cont(x=conf_sub[,get(neigh[2])],nm_var=nice_names[2],nm_level='')) # house value
table2<-rbind(table2, add_cont(x=conf_sub[,get(neigh[3])],nm_var=nice_names[3],nm_level='')) # income

## other variables need to be multiplied by 100 to get %s rather than proportions ##
for (i in 4:length(neigh)){
  table2 <- rbind(table2,add_cont(x=conf_sub[,get(neigh[i])*100],nm_var=nice_names[i],nm_level=''))
}

table2 <- as.data.frame(table2)
names(table2)<-c('Variable','Level','Size','Mean (or N)','SD (or %)')
write.csv(table2, "M:/External Users/KevinJos/output/table2.csv")

# Table 3 -----------------------------------------------------------------

## outcomes ##
outvar_all<-c('fib','newhf', 'newvte', 'mi_acs','iscstroke_tia','death')
out_names<-c('Atrial Fibrilation','Heart Failure','Venous Thromboembolism',
              'Myocardial Infarction with Acute Coronary Syndrome',
              'Ischemic Stroke with Transient Ischemic Attack', 'All-Cause Mortality')

table3 <- data.frame()

for (i in 1:length(outvar_all)){
  
  outvar <- outvar_all[i]
  endtime <- ifelse(outvar == "death", "deathdate", paste0('first_',outvar,'_date'))
  print(outvar)
  
  vnames <- unique(c("bene_id","zip","dob","deathdate", "dtstart","dtend","indexdate",
                     paste0('last_date_',medclass), paste0('first_date_',medclass),endtime))
  
  subcohort <- cohort[,..vnames]
  
  # format dates
  if(outvar == "death"){
    
    subcohort[,c('dob','dtstart','dtend','indexdate', 
                 paste0('last_date_',medclass),
                 paste0('first_date_',medclass), endtime):=
                .(ymd(dob),ymd(dtstart),ymd(dtend),ymd(indexdate),
                  ymd(get(paste0('last_date_',medclass))),
                  ymd(get(paste0('first_date_',medclass))), 
                  ymd(get(endtime)))]
    
  } else {
    
    subcohort[,c('dob','deathdate','dtstart','dtend','indexdate', 
                 paste0('last_date_',medclass),
                 paste0('first_date_',medclass), endtime):=
                .(ymd(dob),ymd(deathdate),ymd(dtstart),ymd(dtend),ymd(indexdate),
                  ymd(get(paste0('last_date_',medclass))),
                  ymd(get(paste0('first_date_',medclass))),
                  ymd(get(endtime)))]
    
  }
  
  pt0 <- with(subcohort, pmin(deathdate, dtend, ymd("2016-11-30"),
                           get(paste0("first_date_", medclass)),
                           get(endtime), na.rm = T) - indexdate)
  pt1 <- with(subcohort, pmin(deathdate, dtend, ymd("2016-11-30"),
                           get(paste0("last_date_", medclass)),
                           get(endtime),na.rm = T) - 
                get(paste0("first_date_", medclass)))
  
  pt0[is.na(pt0) | pt0 < 0] <- 0
  pt1[is.na(pt1) | pt1 < 0] <- 0
  
  pt0 <- as.vector(pt0 - I(pt1>0))/365.25
  pt1 <- as.vector(pt1)/365.25
  
  after <- as.numeric(subcohort[,get(endtime)] >= firstmed &
    subcohort[,get(endtime)] <= enddate)
  after[is.na(after)] <- 0
  
  before <- as.numeric((subcohort[,get(endtime)] < firstmed | is.na(firstmed)) &
    subcohort[,get(endtime)] <= enddate)
  before[is.na(before)] <- 0
  
  all <- before + after
  pt <- pt0 + pt1
  
  n <- round(sum(pt),ndig)
  n_x <- sum(all, na.rm = T)
  pct_x <- round(100*n_x/n, ndig)
  
  n0 <- round(sum(pt0),ndig)
  n_x0 <- sum(before==1,na.rm=T)
  pct_x0 <- round(100*n_x0/n0, ndig)
  
  n1<-round(sum(pt1),ndig)
  n_x1<- sum(after==1,na.rm=T)
  pct_x1<-round(100*n_x1/n1, ndig)
  
  tbl <- c(outvar, n, n_x, pct_x,
           n0, n_x0, pct_x0,
           n1, n_x1, pct_x1)
  
  table3 <- rbind(table3, tbl)
  
}

names(table3)<-c('Outcome','Person Years','Events','Events by Person Year',
                 'Person Years - No Meds','Events - No Meds','Events by Person Year (%) - No Meds',
                 'Person Years - Meds','Events - Meds','Events by Person Year (%) - Meds')

table3 <- as.data.frame(table3)
write.csv(table3, "M:/External Users/KevinJos/output/table3.csv")
