library(haven)
library(data.table)
library(lubridate)
library(stringr)
library(ranger)
library(geepack)
library(survival)

source("M:/External Users/KevinJos/code/ipw_fun.R")

# get the full name of the drug in the dataset (medclass)
medclass <- "oralsteroid"
outvar_all <- c('mi_acs','iscstroke_tia','newvte','newhf','fib','death')

# read in the full dataset
cohort <- read_sas("M:/External Users/Shared/To RachelNet/cohortfinalaug2021_a2c2.sas7bdat")
names(cohort) <- tolower(names(cohort))
setDT(cohort)

# remove sex=0 (there's only one person with this)
cohort <- cohort[!(cohort$sex=='0')]
cohort$sex <- ifelse(cohort$sex == '1', 0, 1)

# Create race indicators
cohort <- cohort[,race:=factor(race,levels=as.character(c(1:6,0)))]

# subset
# set.seed(42)
# init <- sample(unique(cohort$bene_id), ceiling(length(unique(cohort$bene_id))/2), replace = FALSE)
# cohort <- cohort[cohort$bene_id %in% init,]
# rm(init)

# comorbidities and medication/health hx
select_dx <- names(cohort)[grep('dx_',names(cohort))]
select_hx <- c("hx_n_hospvisits", "hx_n_ervisits", "hx_n_outpvisits", "hx_n_meds")

# condition strata
cohort[,card := as.numeric((chronic_fib+chronic_hf+chronic_pvd+chronic_mi+chronic_acs+
                              chronic_mi_acs+chronic_hemstroke+chronic_iscstroke+
                              chronic_tia+chronic_carotid+chronic_cva)>=1)]
cohort[,thromb_chron := as.numeric((chronic_cancer+chronic_vte)>=1)]
cohort[,thromb_acute := as.numeric((chronic_tja+chronic_tha+chronic_tka)>=1)]

cohort[,qc_type:=0]
cohort[card == 1 & thromb_chron == 0 & thromb_acute == 0, qc_type:=1]
cohort[card == 0 & thromb_chron == 1 & thromb_acute == 0, qc_type:=2]
cohort[card == 0 & thromb_chron == 0 & thromb_acute == 1, qc_type:=3]
cohort[card == 0 & thromb_chron == 1 & thromb_acute == 1, qc_type:=4]
cohort[card == 1 & (thromb_chron == 1 | thromb_acute == 1), qc_type:=5]
cohort[,qc_type:=factor(qc_type)]

# remove non-chronic users
# cohort <- subset(cohort, qc_type != 3)

for (i in 1:length(outvar_all)){
  
  # Data Cleaning -----------------------------------------------------------
  
  outvar <- outvar_all[i]
  endtime <- ifelse(outvar == "death", "deathdate", paste0('first_',outvar,'_date'))
  print(outvar)
  
  # pick off only the variables needed for the analyses for this drug
  vnames <- unique(# identifiers and individual characteristics
    c("bene_id","zip","dob","deathdate","age","sex","race","dualeligible",
      "dtstart","dtend","enrol_mths","indexdate","qc_type",
      # pre-existing conditions and hospitalizations
      select_dx, select_hx,
      # medication inf
      paste0('last_date_',medclass), paste0('first_date_',medclass),
      endtime, paste0('newuser_',medclass,'_365d')))
  
  subcohort <- cohort[,..vnames]
  
  # format dates
  if(outvar == "death"){
    
    subcohort[,c('dob','dtstart','dtend','indexdate', paste0('last_date_',medclass),
                 paste0('first_date_',medclass), endtime):=
                .(ymd(dob),ymd(dtstart),ymd(dtend),ymd(indexdate),
                  ymd(get(paste0('last_date_',medclass))),
                  ymd(get(paste0('first_date_',medclass))), ymd(get(endtime)))]
    
    subcohort<-subcohort[,age:=as.double(difftime(indexdate,dob,units='days'))/365.25
                         # add a column with the end date for each person
                         ][,enddate:=pmin(deathdate,dtend,ymd("2016-11-30"),
                                          get(paste0("last_date_", medclass)), na.rm=T)
                           # and a column that says what type of event the end event is
                           ][,enddate_type:=0
                             ][enddate==get(endtime),enddate_type:=1]
    
  } else {
    
    subcohort[,c('dob','deathdate','dtstart','dtend','indexdate', paste0('last_date_',medclass),
                 paste0('first_date_',medclass), endtime):=
                .(ymd(dob),ymd(deathdate),ymd(dtstart),ymd(dtend),ymd(indexdate),
                  ymd(get(paste0('last_date_',medclass))),
                  ymd(get(paste0('first_date_',medclass))),
                  ymd(get(endtime)))]
    
    subcohort<-subcohort[,age:=as.double(difftime(indexdate,dob,units='days'))/365.25
                         # add a column with the end date for each person
                         ][,enddate:=pmin(deathdate,get(endtime),dtend, ymd("2016-11-30"),
                                          get(paste0("last_date_", medclass)), na.rm=T)
                           # and a column that says what type of event the end event is
                           ][,enddate_type:=0
                             ][enddate==deathdate,enddate_type:=2
                               ][enddate==get(endtime),enddate_type:=1]
    
  }
  
  ## PUT DATA IN COUNTING PROCESS FORMAT 
  
  # make a datset that just includes the index date, first med date, censoring/event date, and end date type columns for each person
  baseevents<-subcohort[,.(bene_id,zip,indexdate,firstmed=get(paste0('first_date_',medclass)),
                           lastmed=get(paste0('last_date_',medclass)),enddate,enddate_type)
                        # and remove people who start taking meds or experience a terminating event before their index date
                        ][indexdate<pmin(firstmed,enddate,na.rm=T)
                          # if a person first takes meds after a terminating event, treat as though they never take meds
                          ][!is.na(firstmed) & firstmed>enddate,':='(firstmed=NA,lastmed=NA)]
  
  # make a counting-process style dataset for base events
  cp <- melt(baseevents, measure.vars = c('indexdate','firstmed','lastmed','enddate'),
             variable.factor = F,variable.name = 'date_type',value.name = 'date')
  
  # add in the base events info for each person into the long form data
  cp <- baseevents[cp,on = .(bene_id=bene_id,zip=zip,enddate_type=enddate_type)
                   # remove firstmed and lastmed date for people who never took steroids
                   ][!is.na(date)
                     # create a factor variable for the different event types
                     ][,date_type:=factor(date_type,levels=c('indexdate','firstmed','lastmed','enddate'))
                       # order by beneficiary id and then event type
                       ][order(bene_id,date_type)]
  
  ## ADD IN INFO ABOUT ZIPCODE CHANGES 
  
  # read in data with info on zipcode changes
  changezip <- read_sas("M:/External Users/Shared/To RachelNet/morethanonezipcode.sas7bdat")
  names(changezip) <- tolower(names(changezip))
  setDT(changezip)
  setnames(changezip,c('year','zip'),c('yr_czip','zipnew'))
  changezip[,date_czip:=ymd(paste0(yr_czip,'-01-01'))]
  
  # remove zipcode changes that occur before the person's index date
  changezip <- changezip[indexdate>date_czip,date_czip:=indexdate][order(bene_id,-yr_czip)]
  changezip <- changezip[!duplicated(changezip[,.(bene_id,date_czip)])][,yr_czip:=year(date_czip)]
  
  # for each individual with a zipcode change, make a sequence of years spanning the entire period they're in the study
  baseevents_zc <- baseevents[bene_id %in% changezip[,bene_id]][,.(yr_czip=year(indexdate):year(enddate)),by=bene_id]
  
  # merge zipcode change data with start/end date info for the person
  changezip<-changezip[baseevents_zc,on=.(bene_id,yr_czip)
                       # carry forward zips for missings
                       ][order(bene_id,yr_czip)
                         ][,zipnew_ff:=zipnew[1],.(bene_id,cumsum(!is.na(zipnew)))]
  
  # merge this dataset with the counting process dataset to get the right zipcode for each year that these people have an event
  cp[,yr_czip:=year(date)]
  cp <- merge(cp,changezip[,.(bene_id,yr_czip,zipnew_ff)], by=c('bene_id','yr_czip'),all.x=T)[
    !is.na(zipnew_ff),zip:=zipnew_ff][,'zipnew_ff':=NULL]
  
  # now add rows to the counting process data for each zip change (but not for first zip)
  addRows <- changezip[(!is.na(date_czip)) & (!indexdate==date_czip),!c('indexdate','zipnew_ff')][,date_type:='changezip']
  setnames(addRows,c('zipnew','date_czip'),c('zip','date'))
  cpCols <- c(which(names(cp)=='bene_id'),which(!(names(cp) %in% names(addRows))))
  addRows <- merge(addRows,cp[date_type=='indexdate',..cpCols],by='bene_id')
  
  # put cp and changezip rows together
  cp <- rbind(cp,addRows)
  cp <- cp[date_type != "changezip" | date != enddate ]
  
  ## ADD IN SEASON INFORMATION 
  
  getSeason <- function(input_date){
    input_month <- data.frame('id'=1:length(input_date),'inmonth'=month(input_date))
    season <- data.frame('inmonth'=c(12,1:11),'season'=rep(c('winter','spring','summer','fall'),each=3))
    final <- merge(input_month,season,by='inmonth',all.x=T)
    final <- final[order(final$id),]
    return(final$season)
  }
  
  # add season
  cp[,season:=getSeason(date)
     # also add adjusted year for linkage with pm data (december is grouped with winter of following year in pm data)
     ][,yr_ssn:=year(date)
       ][month(date)==12,yr_ssn:=yr_ssn+1]
  
  ## LINK PM VALUES WITH BASE EVENTS AND ZIP CHANGES 
  
  # read in the seasonal average pollutants
  pm <- fread("M:/External Users/RachelNet/data/pollution/pm25_seasonalavg_zipcode.csv")
  pm <- melt(pm,measure.vars = c('pm25_winter','pm25_spring','pm25_summer','pm25_fall'),
             variable.factor = F,variable.name = 'season',value.name = 'pm')
  pm$ssn_num <- with(pm, ifelse(season == 'pm25_winter', 1, ifelse(season == 'pm25_spring', 2,ifelse(season == 'pm25_summer', 3, 4))))
  pm <- pm[order(ZIP, year,ssn_num)]
  pm <- pm[,c('pm.lag1','pm.lag2','pm.lag3', 'pm.lag4') := 
             .(shift(pm, 1, type = "lag"),shift(pm, 2, type = "lag"),
               shift(pm, 3, type = "lag"),shift(pm, 4, type = "lag"))]
  pm$ssn_num <- NULL
  pm <- pm[year>=2008]
  setnames(pm,c('ZIP','year'),c('zip','yr_ssn'))
  pm[,':='(season=substr(season,start = 6,stop=nchar(season)),zip=str_pad(zip,width = 5,side='left',pad='0'))]
  
  # merge with health data
  cp <- pm[cp,on=.(zip,yr_ssn,season)]
  
  ## ADD ROWS FOR EACH SEASON CHANGE WITH NEW PM VALUES 
  
  # create a dataset that changes the PM on the first date of each new season
  changepm <- expand.grid('season'=c('winter','spring','summer','fall'),'yr_ssn'=min(pm$yr_ssn):max(pm$yr_ssn))
  md_ssn <- data.frame('season'=c('winter','spring','summer','fall'),'md'=c('12-01','03-01','06-01','09-01'))
  changepm <- setDT(merge(changepm,md_ssn,by='season'))
  changepm[,yrA:=yr_ssn][season=='winter',yrA:=yr_ssn-1][,changedate:=ymd(paste(yrA,md,sep='-'))][,c('md','yrA'):=NULL]
  changepm <- merge(changepm,pm,by=c('season','yr_ssn'))
  
  # merge by zipcode with the index/end dates for each person and then only keep rows with changedate between index and enddate
  cpMerge <- cp[date_type %in% c("indexdate", "changezip", "enddate"),
                ][,c("bene_id","zip","date", "date_type")
                  ][order(bene_id,date)
                    ][,date2 := shift(date, type='lead'), by = bene_id
                      ][date_type != "enddate",]
  
  cpAir <- merge(changepm, cpMerge, by=c('zip'), allow.cartesian = T)[
    data.table::between(changedate,date,date2)
    ][,':='(date=changedate,date_type='changePM',yr_czip=year(changedate))
      ][,c('date2','changedate'):=NULL]
  
  cpCols <- c(which(names(cp)=='bene_id'),which(!(names(cp) %in% names(cpAir))))
  cpAir <- merge(cpAir,cp[date_type=='indexdate',..cpCols],by='bene_id')
  
  # put cp and cpAir together
  cp <- rbind(cp,cpAir)
  
  rm(pm, changepm, baseevents, cpAir, cpCols, cpMerge, addRows,
     baseevents_zc, changezip, md_ssn, getSeason); gc()
  
  # add nicer event names
  eventMat <- data.table('date_type'=c('indexdate','firstmed','lastmed','enddate','changezip','changePM'),
                         'event'=c('enterStudy','startMed','endMed','leaveStudy','changeZip','changePM'))
  cpAll <- merge(cp,eventMat,by='date_type')
  cpAll[,'date_type':=NULL]
  
  rm(cp, eventMat, vnames);  gc()
  
  ## CLEAN UP COUNTING PROCESS STRUCTURE
  
  # counting process dataset can contain multiple rows corresponding to events on same day for a single person, must remove these, prioritizing base events
  cpAll <- cpAll[,event:=factor(event,levels=c('leaveStudy','startMed','endMed','enterStudy','changeZip','changePM'))
                 ][order(bene_id,event)][,dup:=duplicated(date),by=bene_id][dup==0]
  
  # make a binary "on meds" variable
  cpFit <- cpAll[,onMeds:=as.numeric(date>=firstmed & date<lastmed)
                 ][is.na(onMeds),onMeds:=0
                   # make a binary "failed" and "censored" variable
                   ][,':='(died=as.numeric(event=='leaveStudy' & enddate_type==2),
                           failed=as.numeric(event=='leaveStudy' & enddate_type==1))
                     # add time variables to the dataset for the survival function
                     ][,time0:=as.double(difftime(date,indexdate,units='days'))
                       ][order(bene_id,time0)
                         ][,':='(time1=shift(time0,n=1,type='lead'),
                                 died=shift(died,n=1,type='lead'),
                                 failed=shift(failed,n=1,type='lead')), by=bene_id
                           ][order(bene_id,time0,time1)][date <= enddate]
  
  rm(cpAll); gc()
  
  # Inspection
  # tempdat <- with(cpFit, data.table(bene_id, indexdate, date, enddate, season, zip,
  #                                   onMeds, event, failed, died, time0, time1))[order(bene_id, time0, time1)]
  
  ## CREATE DRUG PANELS
  
  # define ordered season time
  cpFit$ssn_time <- as.numeric(4*(cpFit$yr_ssn - 2008) + (cpFit$season == "spring") +
                                 2*(cpFit$season == "summer") + 3*(cpFit$season == "fall")) + 1
  
  # shifting seasons
  setDT(cpFit)
  cpFit_premed <- subset(cpFit, event == "startMed", select = c(bene_id, yr_ssn, season, ssn_time, date))
  cpFit_last <- subset(cpFit, event == "leaveStudy", select = c(bene_id, time0, enddate_type))
  
  cpFit_premed$shift0 <- with(cpFit_premed, ifelse(season == "winter", difftime(date, as.Date(paste0(yr_ssn-1,"-12-01")), units='days'),
                                                   ifelse(season == "spring", difftime(date, as.Date(paste0(yr_ssn,"-03-01")), units='days'),
                                                          ifelse(season == "summer", difftime(date, as.Date(paste0(yr_ssn,"-06-01")), units='days'),
                                                                 difftime(date, as.Date(paste0(yr_ssn,"-09-01")), units='days')))))
  
  cpFit_premed$shift1 <- with(cpFit_premed, ifelse(season == "winter", difftime(date, as.Date(paste0(yr_ssn,"-03-01")), units='days'),
                                                   ifelse(season == "spring", difftime(date, as.Date(paste0(yr_ssn,"-06-01")), units='days'),
                                                          ifelse(season == "summer", difftime(date, as.Date(paste0(yr_ssn,"-09-01")), units='days'),
                                                                 difftime(date, as.Date(paste0(yr_ssn,"-12-01")), units='days')))))
  
  cpFit_premed$shift <- with(cpFit_premed, ifelse(ssn_time == 1 | shift0 > abs(shift1), shift1, shift0))
  
  # merge in shifts
  setDT(cpFit_premed)
  setDT(cpFit_last)
  
  cpFit <- merge(cpFit, data.frame(bene_id = cpFit_premed$bene_id, 
                                   shift = cpFit_premed$shift),
                 by = "bene_id", all.x = TRUE)
  
  cpFit <- merge(cpFit, data.frame(bene_id = cpFit_last$bene_id, last = cpFit_last$time0,
                                   med_died = as.numeric(cpFit_last$enddate_type==2),
                                   med_failed = as.numeric(cpFit_last$enddate_type==1)),
                 by = "bene_id", all.x = TRUE)
  
  rm(cpFit_last, cpFit_premed); gc()
  
  # partition meds/no meds
  cpFit_drug <- subset(setDT(cpFit), !is.na(shift) & event == "changePM")
  cpFit_meds <- subset(setDT(cpFit), !is.na(shift))
  cpFit_nomeds <- subset(setDT(cpFit), is.na(shift))[!is.na(time1) & time0 < last & time0 != time1 & time0 >= 0]
  cpFit_meds$time1[is.na(cpFit_meds$time1)] <- cpFit_meds$time0[is.na(cpFit_meds$time1)]
  setDT(cpFit_drug); setDT(cpFit_meds); setDT(cpFit_nomeds)
  
  # shift time0 around and update time1
  cpFit_drug$time0 <- with(cpFit_drug, time0 + shift, time0)
  cpFit_drug$event <- "changeDrug"
  cpFit_drug$date <- ymd(cpFit_drug$date + cpFit_drug$shift)
  cpFit_meds <- setDT(rbind(cpFit_drug, cpFit_meds))
  cpFit_meds <- cpFit_meds[order(bene_id, date, -event)
                           ][,time1:=ifelse(event != "leaveStudy", shift(time0, n = 1, type = "lead"), time1), by = bene_id
                             ][!is.na(time1) & time0 < last & time0 != time1 & time0 >= 0]
  cpFit_meds$failed <- with(cpFit_meds, ifelse(time1 == last & med_failed == 1, 1, 0))
  cpFit_meds$died <- with(cpFit_meds, ifelse(time1 == last & med_died == 1, 1, 0))
  
  cpFit_nomeds$drug_time <- cpFit_nomeds$ssn_time
  cpFit_meds$drug_time <- cpFit_meds$ssn_time
  
  falseifNA <- function(x){ ifelse(is.na(x), FALSE, x) }
  ifelse2 <- function(x, a, b){ ifelse(falseifNA(x), a, b) }
  
  # ensure counting process indexes and measurements are consistent
  cpFit_meds <- cpFit_meds[order(bene_id, time0, time1)
                           ][,season:=ifelse(event %in% c("changeDrug", "startMed") & !is.na(shift(season, n = 1, type = "lag")),
                                             shift(season, n = 1, type = "lag"), season), by = bene_id]
  
  cpFit_meds <- cpFit_meds[order(bene_id, time0, time1)
                           ][,yr_ssn:=ifelse(event %in% c("changeDrug", "startMed") & !is.na(shift(yr_ssn, n = 1, type = "lag")),
                                             shift(yr_ssn, n = 1, type = "lag"), yr_ssn), by = bene_id]
  
  cpFit_meds <- cpFit_meds[order(bene_id, time0, time1)
                           ][,pm:=ifelse(event %in% c("changeDrug", "startMed") & !is.na(shift(pm, n = 1, type = "lag")),
                                         shift(pm, n = 1, type = "lag"), pm), by = bene_id]
  
  cpFit_meds <- cpFit_meds[order(bene_id, time0, time1)
                           ][,drug_time:=ifelse(event %in% c("changeDrug", "startMed") & !is.na(shift(drug_time, n = 1, type = "lag")),
                                                shift(drug_time, n = 1, type = "lag") + 1, drug_time), by = bene_id]
  
  cpFit_meds <- cpFit_meds[order(bene_id, time0, time1)
                           ][,drug_time:=ifelse2(event  == "enterStudy" & shift(event, n = 1, type = "lead") == "changePM",
                                                 shift(drug_time, n = 1, type = "lead"), drug_time), by = bene_id]
  
  # define ordered season time
  cpFit_meds$ssn_time <- as.numeric(4*(cpFit_meds$yr_ssn - 2008) + (cpFit_meds$season == "spring") +
                                      2*(cpFit_meds$season == "summer") + 3*(cpFit_meds$season == "fall")) + 1
  
  # Inspection
  # tempdat <- with(cpFit_meds, data.table(bene_id, indexdate, date, enddate, season, zip, drug_time, ssn_time,
  #                                        onMeds, pm, event, failed, died, time0, time1))[order(bene_id, time0, time1)]
  
  # put it all together and reset pm, yr_ssn, and season
  cpFit <- setDT(rbind(cpFit_meds, cpFit_nomeds))
  
  rm(cpFit_meds, cpFit_nomeds, cpFit_drug); gc()
  
  ## ADD IN ZIPCODE LEVEL CONFOUNDERS 
  
  cpFit <- cpFit[order(bene_id, time0, time1)
                 ][,zip:=ifelse(event != "changeZip" & !is.na(shift(zip, n = 1, type = "lag")),
                                shift(zip, n = 1, type = "lag"), zip), by = bene_id
                   ][order(bene_id, time0, time1)
                     ][,drug_time:=ifelse2(event == "changeZip" & shift(event, n = 1, type = "lead") == "changePM",
                                           shift(drug_time, n = 1, type = "lead"), drug_time), by = bene_id
                       ][,drug_time:=ifelse2(event == "changeZip" & shift(event, n = 1, type = "lead") %in% c("changeDrug", "starMed"),
                                             shift(drug_time, n = 1, type = "lead") - 1, drug_time), by = bene_id]
  
  # read in and clean the zipcode level confounders
  conf <- fread("M:/External Users/RachelNet/data/confounders/census_interpolated_zips.csv")
  
  # subset to only 2008 and later (dates when medicare data available)
  conf <- conf[year>=2008
               # make zipcode into a 5-digit character string
               ][,zip:=str_pad(ZIP,width=5,side='left',pad='0')
                 # remove unnecessary columns
                 ][,c('V1','zcta','ZIP'):=NULL]
  setnames(conf,c('year','hispanic'),c('yr_ssn','pct_hisp'))
  
  # merge confounders with health data by zipcode and year (either index year or year of zipcode change)
  cpFit <- merge(cpFit, conf, by=c('zip','yr_ssn'))
  
  ## MERGE THE CONSTANT CHARACTERISTICS BACK IN
  
  vnames <- unique(c("bene_id","age","sex","race","dualeligible","qc_type",
                     names(subcohort)[grep('dx_',names(subcohort))],
                     names(subcohort)[grep('hx_',names(subcohort))]))
  cpFit <- merge(cpFit,subcohort[,..vnames],by=c('bene_id'))
  
  rm(conf, subcohort, vnames); gc()
  
  # person-season
  cpFit$bene_id_ssn <- paste(cpFit$bene_id, cpFit$ssn_time, sep = "-")
  cpFit$zip_ssn <- paste(cpFit$zip, cpFit$ssn_time, sep = "-")
  cpFit$bene_id_drug <- paste(cpFit$bene_id, cpFit$drug_time, sep = "-")
  
  # time-varying age
  cpFit$age_tm <- with(cpFit, as.double(difftime(date,indexdate,units='days')/365.25) + age)
  
  # remove missing pm measurements among other covariates
  cpFit <- cpFit[complete.cases(model.frame(data = setDF(cpFit), na.action = NULL,
                                            formula = formula(paste0('~pm+pm.lag1+pm.lag2+pm.lag3+pm.lag4+
                                                                      onMeds+age+sex+race+dualeligible+
                                                                      season+poverty+popdensity+medianhousevalue+
                                                                      pct_owner_occ+education+medhouseholdincome+
                                                                      pct_hisp+pct_blk+pct_white+',
                                                                     paste0(select_dx,collapse = '+'), ' + ',
                                                                     paste0(select_hx,collapse = '+'))))),]
  
  rm.id <- unique(cpFit[which((cpFit$time1 - cpFit$time0) > 137),]$bene_id)
  cpFit <- cpFit[!(cpFit$bene_id %in% rm.id),]
  
  # inspect/check that code works
  # tempdat <- with(cpFit, data.table(bene_id, indexdate, enddate, zip, date, ssn_time, drug_time, onMeds,
  #                                   pm, event, failed, died, time0, time1, shift, last))[order(bene_id, time0, time1)]
  
  # IPW Analysis ------------------------------------------------------------
  
  ## PM IPW WEIGHTS
  
  setDT(cpFit)
  pmFit <- cpFit[!duplicated(zip_ssn)]
  setDF(pmFit); gc()
  
  fmla.denom.pm<-formula(paste0('~ pm.lag1+pm.lag2+pm.lag3+pm.lag4+
                                    season+yr_ssn+poverty+popdensity+
                                    medianhousevalue+pct_owner_occ+education+
                                    medhouseholdincome+pct_hisp+pct_blk+pct_white'))
  
  pm_w <- ipwtm_ranger(exposure = "pm", denominator = fmla.denom.pm, 
                       numerator = formula('~pm.lag1+pm.lag2+pm.lag3+pm.lag4'), 
                       id = "zip",timevar = "ssn_time", data = pmFit, trunc = 0.05, 
                       continuous = TRUE, num.trees = 200, max.depth = 8, num.threads = 12)
  
  cpFit <- merge(cpFit, data.frame(pmFit[,c("zip","ssn_time")], ipw_pm = pm_w$weights.trunc),
                 by = c("zip","ssn_time"), all.x = TRUE)
  
  rm(fmla.denom.pm, pmFit, pm_w); gc()
  
  ## STEROID IPW WEIGHTS
  
  setDT(cpFit)
  medFit <- cpFit[cpFit[,.I[onMeds == max(onMeds)], by = bene_id_drug]$V1]
  medFit <- medFit[!duplicated(bene_id_drug)]
  medFit <- medFit[order(bene_id_drug)]
  medFit <- medFit[,onMeds.lag:=shift(onMeds, n = 1, type = "lag"), by = bene_id]
  medFit <- subset(medFit, onMeds.lag == 0 | is.na(onMeds.lag))
  setDF(medFit); gc()
  
  fmla.denom.med <- formula(paste0('~ pm+age_tm+sex+race+dualeligible+season+yr_ssn+poverty+
                                    popdensity+medianhousevalue+pct_owner_occ+education+
                                    medhouseholdincome+pct_hisp+pct_blk+pct_white+',
                                   paste0(select_hx, collapse = '+'), '+',
                                   paste0(select_dx, collapse = '+')))
  
  med_w <- ipwtm_ranger(exposure = "onMeds", denominator = fmla.denom.med, numerator = formula("~ pm"),
                        id = "bene_id",timevar = "drug_time", data = medFit, trunc = 0.01, 
                        continuous = FALSE, num.trees = 200, max.depth = 8, num.threads = 12)
  
  cpFit <- merge(cpFit, data.frame(medFit[,c("bene_id","drug_time")], ipw_med = med_w$weights.trunc),
                 by = c("bene_id","drug_time"), all.x = TRUE)
  
  cpFit$ipw_med[is.na(cpFit$ipw_med)] <- 1
  
  rm(med_w, fmla.denom.med, medFit); gc()
  
  ## CENSOR IPW WEIGHTS
  
  if (outvar != "death") {
    
    setDT(cpFit)
    censorFit <- cpFit[cpFit[,.I[died == max(died)], by = bene_id_drug]$V1]
    censorFit <- censorFit[!duplicated(bene_id_drug)]
    setDF(censorFit); gc()
    
    fmla.denom.censor <- formula(paste0('~ pm+onMeds+age_tm+sex+race+dualeligible+season+yr_ssn+poverty+
                                        popdensity+medianhousevalue+pct_owner_occ+education+
                                        medhouseholdincome+pct_hisp+pct_blk+pct_white+',
                                        paste0(select_hx, collapse = '+'), '+',
                                        paste0(select_dx, collapse = '+')))
    
    censor_w <- ipwtm_ranger(exposure = "died", denominator = fmla.denom.censor, numerator = formula("~ onMeds+pm"),
                             id = "bene_id", timevar = "drug_time", data = censorFit, trunc = 0.01,
                             continuous = FALSE, censor = TRUE, num.trees = 200, max.depth = 8, num.threads = 12)
    
    cpFit <- merge(cpFit, data.frame(censorFit[,c("bene_id","drug_time")], ipw_censor = censor_w$weights.trunc),
                   by = c("bene_id","drug_time"), all.x = TRUE)
    
    rm(censor_w, fmla.denom.censor, censorFit); gc()
    
    
  } else {
    
    cpFit$ipw_censor <- 1
    
  }
  
  # put it all together
  dat <- setDT(subset(cpFit, select = c(bene_id, date, zip, season, yr_ssn, ssn_time, drug_time,
                                        qc_type, age, onMeds, pm, event, failed, died,
                                        time0, time1, ipw_pm, ipw_med, ipw_censor)))
  
  dat <- dat[order(bene_id,time0)]
  dat$time0 <- with(dat, time0/365.25 + age)
  dat$time1 <- with(dat, time1/365.25 + age)
  dat$pm_med <- ifelse(dat$onMeds == 1, dat$pm, 8)
  dat$pm_nomed <- ifelse(dat$onMeds == 0, dat$pm, 8)
  
  save(dat, file=paste0("M:/External Users/KevinJos/output/fit_data/", medclass, '_', outvar,'.RData'))
  
  rm(dat, cpFit); gc()
  
}

rm(cohort); gc()
