rm(list=ls())

library(tidyverse)
library(lubridate)
library(haven)
library(epiR)
library(arsenal)
library(anthro)
library(zscorer)

dat <- read_dta("data/Pneumonia_influenza.dta") 

##### data management #####

# drop names
dat$idfathfn <- NULL
dat$idfathln <- NULL
dat$idmothfn <- NULL
dat$idmothln <- NULL
dat$flun <- NULL
dat$map <- NULL
dat$fathvita <- NULL
nam <- c("sib1dob","sib1tw","sib1sex","sib1age","sib2bobuk","sib2dob","sib2tw","sib2sex",
         "sib2age","sib3sex","sib3tw","sib2age3","sib3dob")
dat <- dat[,!names(dat) %in% nam]

dat$spin[dat$spin==""]<-NA

### define groups ### 

dat$hiveu <- NA
dat$hiveu[dat$hivstatus==0 & dat$mhivres==1]<-1
dat$hiveu[dat$hivstatus==1 | dat$mhivres==2]<-0

dat$hivuu<-NA
dat$hivuu[dat$hivstatus==0 & dat$mhivres==2]<-1
dat$hivuu[dat$hivstatus==1 | dat$mhivres==1]<-0

# children not classified
dat$excchild <- NA
dat$excchild[is.na(dat$hivstatus) & (is.na(dat$mhivres) | dat$mhivres==3 | dat$mhivres==9)] <- 1
dat$excchild[dat$hivstatus==0 & (dat$mhivres ==3 | dat$mhivres==9)] <- 1
dat$excchild[dat$hivstatus==1 & (dat$mhivres ==3 | dat$mhivres==9)] <- 1

exdat <- dat[which(dat$excchild==1),]
m <- which(dat$hiveu==0 | dat$hivuu==0)

# create a variable seperating the exposure status 
# for creating contingency tables

#dat$clasif <- NA
j <- which(dat$hivuu==1)
k <- which(dat$hiveu==1)
dat[j,"clasif"] <- 1
dat[k,"clasif"] <- 2

## tidying up some variables
dat$thirst[!is.na(dat$thirst) & dat$thirst==9]<-NA
dat$length[!is.na(dat$length) & dat$length>500]<-NA

## convert string dates to date format ##
dat$DOB <- dmy(dat$dob)
dat$ca_date <- dmy(dat$cadate)
dat$dischargeDate <- dmy(dat$dcdate)
dat$admissionDate <- dmy(dat$admitdate)

## vaccination dates
dat$bcgdate <- dmy(dat$bcgdate)
dat$p0date <- dmy(dat$p0date)
dat$p1date <- dmy(dat$p1date)
dat$dpthh1date <- dmy(dat$dpthh1date)
dat$pcv1date <- dmy(dat$pcv1date)
dat$rv1date <- dmy(dat$rv1date)
dat$p2date <- dmy(dat$p2date)
dat$p3date <- dmy(dat$p3date)
dat$dpthh2date <- dmy(dat$dpthh2date)
dat$dpthh3date <- dmy(dat$dpthh3date)
dat$pcv2date <- dmy(dat$pcv2date)
dat$pcv3date <- dmy(dat$pcv3date)
dat$rv2date <- dmy(dat$rv2date)
dat$measdate <- dmy(dat$measdate)


## calculate time from birth to different vaccines

dat$BCGdateDays <- difftime(dat$bcgdate,dat$DOB,units = "days") # bcg
dat$BCGdateDays[dat$BCGdateDays < 0] <- NA

dat$polio0Days <- difftime(dat$p1date,dat$DOB,units = "days")  # polio 0 at birth
dat$polio0Days[dat$polio0Days < 0] <- NA

# 6 weeks
dat$polio1Days <- difftime(dat$p1date,dat$DOB,"weeks") # polio
dat$polio1Days[dat$polio1Days < 0] <- NA
dat$DPTDays <- difftime(dat$dpthh1date,dat$DOB,units = "weeks") # diptheria
dat$DPTDays[dat$DPTDays < 0] <- NA
dat$PCV1Days <- difftime(dat$pcv1date,dat$DOB,units = "weeks") # pcv
dat$PCV1Days[dat$PCV1Days < 0] <- NA
dat$RV1Days <- difftime(dat$rv1date,dat$DOB,units = "weeks") # rotavirus
dat$RV1Days[dat$RV1Days < 0] <- NA
###  VACCINATIONS AT 10 WEEKS ####

dat$polio2Days <- difftime(dat$p2date,dat$DOB,units = "weeks") # polio 2
dat$polio2Days[dat$polio2Days < 0] <- NA
dat$DPT2Days <- difftime(dat$dpthh2date,dat$DOB,units = "weeks") # diptheria2
dat$DPT2Days[dat$DPT2Days < 0] <- NA

dat$PCV2Days <- difftime(dat$pcv2date,dat$DOB,units = "weeks") # pcv2
dat$PCV2Days[dat$PCV2Days < 0] <- NA
dat$RV2Days <- difftime(dat$rv2date,dat$DOB,units = "weeks") # rotavirus 2
dat$RV2Days[dat$RV2Days < 0] <- NA
#### VACCINATIONS AT 14 WEEKS ####
dat$polio3Days <- difftime(dat$p3date,dat$DOB,units = "weeks") # polio 3
dat$polio3Days[dat$polio3Days < 0] <- NA
dat$DPT3Days <- difftime(dat$dpthh3date,dat$DOB,units = "weeks") # diptheria 3
dat$DPT3Days[dat$DPT3Days < 0] <- NA
dat$PCV3Days <- difftime(dat$pcv3date,dat$DOB,units = "weeks") # pcv3
dat$PCV3Days[dat$PCV3Days < 0] <- NA
### VACCINATIONS AT 9 MONTHS ####
dat$measles1Days <- difftime(dat$measdate,dat$DOB,units = "days")/30.41667 # measles
dat$measles1Days[dat$measles1Days < 0] <- NA

# define age
j <- which(dat$ca_date<dat$DOB) # isolate non-sensible dates
dat[j,"ca_date"] <- NA
dat$Age <- as.numeric(dat$ca_date-dat$DOB)/30.41667
dat$AgeWhol <- floor(dat$Age)

# date of birth and age of mother
dat$mothdob[dat$mothdob==""] <- NA
dat$mothdob <- gsub("-","",dat$mothdob)
dat$moth_dob <- as.Date(dat$mothdob,format = "%d%b%y")
# define severe pneumonia


# create age groups

k <- which(dat$AgeWhol < 0) 
dat[k,"AgeWhol"] <- NA 
dat$age_grp <- NA
dat$age_grp[dat$AgeWhol <= 2] <- 1
dat$age_grp[dat$AgeWhol > 2 & dat$AgeWhol <= 6] <- 2
dat$age_grp[dat$AgeWhol > 6 & dat$AgeWhol <= 12] <- 3
dat$age_grp[dat$AgeWhol > 12 & dat$AgeWhol <= 24] <- 4
dat$age_grp[dat$AgeWhol > 24 & dat$AgeWhol <= 60] <- 5

dat$sex12<-dat$sex # recording sex; 1=male,2=female
dat$sex12[!is.na(dat$sex)]<-ifelse(dat$sex[!is.na(dat$sex)]==1,1,2)

## compute weight for age and weight for height

dat$whoWeightForAgeZ<-NA
idx<-which(!is.na(dat$sex12) & !is.na(dat$weight) & !is.na(dat$Age))
for(i in idx){
  dat$whoWeightForAgeZ[i]<-try(getWGS(index="wfa",sexObserved=dat$sex12[i],firstPart=dat$weight[i],secondPart=dat$AgeWhol[i]),silent = T) 
  #http://www.who.int/childgrowth/standards/weight_for_age/en/,
}

dat$whoWeightForAgeZ<-as.numeric(dat$whoWeightForAgeZ)

# weight for height
dat$height <- dat$length/100 # problem calculation here, need to be revised
#dat$heightWhol <- ceiling(dat$height)
#dat$whoWeightForHeightZ <- NA
#xx <- which(!is.na(dat$sex12) & !is.na(dat$weight) & !is.na(dat$height))
#for(i in xx){
#  dat$whoWeightForHeightZ[i]<-try(getWGS(index="wfa",sexObserved=dat$sex12[i],firstPart=dat$weight[i],secondPart=dat$heightWhol[i])) 
  #http://www.who.int/childgrowth/standards/weight_for_age/en/,
#}

# define weight for height categories

#dat$wforHght <- NA
#dat$wforHght[dat$whoWeightForHeightZ > -2] <- 1
#dat$wforHght[dat$whoWeightForHeightZ > -3 & dat$whoWeightForHeightZ < -2] <- 2
#dat$wforHght[dat$whoWeightForHeightZ < -3] <- 3

## weight for height for calculating predictive diagnistics - dichotomized
#dat$zscoreHght <- NA
#dat$zscoreHght[!is.na(dat$whoWeightForHeightZ) & dat$whoWeightForHeightZ > -3 & dat$whoWeightForHeightZ < -2] <- 1
#dat$zscoreHght[!is.na(dat$whoWeightForHeightZ) & dat$whoWeightForHeightZ < -3] <- 2

# revise computation of z-scores
dat$heightcm <- dat$height*100
newzscores <- anthro_zscores(sex = dat$sex12,age = dat$AgeWhol,is_age_in_month = TRUE,weight = dat$weight,lenhei = dat$heightcm)

# add new scores to the data
dat$zweightAge <- newzscores$zwei
dat$zweightHeight <- newzscores$zwfl

# define weight for height categories
# replace old values
dat$zwforHght <- NA
dat$zwforHght[dat$zweightHeight >= -2] <- 1
dat$zwforHght[dat$zweightHeight >= -3 & dat$zweightHeight < -2] <- 2
dat$zwforHght[dat$zweightHeight < -3] <- 3

# vaccination schedules

dat$birthVacc <- apply(dat[,c("BCGdateDays","polio0Days")],1,median,na.rm=TRUE) # at birth
dat$sixwksVacc <- apply(dat[,c("polio1Days","DPTDays","PCV1Days","RV1Days")],1,median,na.rm=TRUE)
dat$tenwksVacc <- apply(dat[,c("polio2Days","DPT2Days","PCV2Days","RV2Days")],1,median,na.rm=TRUE)
dat$fourteenwksVacc <- apply(dat[,c("polio3Days","DPT3Days","PCV3Days")],1,median,na.rm=TRUE)

# variable to show whether vaccines are up to date

# consult with pui-ying on definition
dat$epiVacHue <- ifelse(!is.na(dat$BCGdateDays) & dat$bcgdate<dat$ca_date,1,
                    ifelse(!is.na(dat$polio0Days) & dat$p0date<dat$ca_date,1,
                    ifelse(!is.na(dat$polio1Days) & dat$p1date<dat$ca_date,1,
                    ifelse(!is.na(dat$DPTDays) & dat$dpthh1date<dat$ca_date,1,
                    ifelse(!is.na(dat$PCV1Days) & dat$pcv1date<dat$ca_date,1,
                    ifelse(!is.na(dat$RV1Days) & dat$rv1date<dat$ca_date,1,
                    ifelse(!is.na(dat$polio2Days) & dat$p2date<dat$ca_date,1,
                    ifelse(!is.na(dat$DPT2Days) & dat$dpthh2date<dat$ca_date,1,
                    ifelse(!is.na(dat$PCV2Days) & dat$pcv2date<dat$ca_date,1,
                    ifelse(!is.na(dat$RV2Days) & dat$rv2date<dat$ca_date,1,
                    ifelse(!is.na(dat$measles1Days) & dat$measdate<dat$ca_date,1,0)))))))))))

# check for the presence of danger signs
dat$dangerSigns <- ifelse(!is.na(dat$thirst) & dat$thirst==3,1,
                        ifelse(!is.na(dat$vom) & dat$vom==1,1,
                        ifelse(!is.na(dat$diarr) & dat$diarr==1,1,
                        ifelse(!is.na(dat$conv) & dat$conv==1,1,
                        ifelse(!is.na(dat$rash) & dat$rash==1,1,
                        ifelse(!is.na(dat$oedema) & dat$oedema==1,1,0))))))

# length of hospitsl stay                  
dat$hosplen <- as.numeric(dat$dischargeDate - dat$admissionDate)
dat$hosplen[dat$hosplen < 0] <- NA
# fever
dat$hfever <- ifelse(dat$temp >=38,1,0) 

# household size
dat$hh_size1 <- NA
dat$hh_size1[!is.na(dat$hhsize)] <- ifelse(dat$hhsize[!is.na(dat$hhsize)] > 4,1,0)

# muac scores
dat$muac.score <- ifelse(!is.na(dat$muac) & dat$muac < 11.5,1,0)
dat$muac.scorev2 <- ifelse(dat$muac < 11.5,1,0)

# low birth weight
dat$LBW <- ifelse(dat$bwt < 2.5,1,0) # low birth weight

# antibiotic use in the last 7 days
dat$antib.use <- ifelse(dat$antib==1,1,0) 

# age adjusted tachepenia - < 12 months
dat$tachpenia <- NA
dat$tachpenia[dat$resps > 50 & dat$AgeWhol <= 12] <- 1
dat$tachpenia[dat$resps > 40 & dat$AgeWhol > 12 & dat$AgeWhol <=60] <- 2

# age - adjusted tachycardia
dat$tachycardia <- NA
dat$tachycardia[dat$pulse > 160 & dat$AgeWhol <= 12] <- 1
dat$tachycardia[dat$pulse > 150 & dat$AgeWhol > 12 & dat$AgeWhol <= 36] <- 2
dat$tachycardia[dat$pulse > 140 & dat$AgeWhol > 36 & dat$AgeWhol <= 60] <- 3

# hypoxemia
dat$hypo <- NA
dat$hypo[!is.na(dat$sats)] <- ifelse(dat$sats[!is.na(dat$sats)] < 90,1,0)

# heptomegaly
dat$hepgly <- dat$hepgly
dat$hepgly[dat$hepgly==2] <- 0
#dat$hepgly <- ifelse(!is.na(dat$hepgly) & dat$hepgly==1,1,0)

## stridor
dat$strid <- dat$strid
dat$strid[dat$strid==2] <- 0
#dat$strid <- ifelse(!is.na(dat$strid) & dat$strid==1,1,0)

## lower chest indrawing
#dat$icid <- ifelse(!is.na(dat$icid) & dat$icid==1,1,0) # revise
dat$icid <- dat$icid
dat$icid[dat$icid==2] <- 0

## oral thrush
dat$thrush <- dat$thrush
dat$thrush[dat$thrush==2] <- 0
#dat$thrush <- ifelse(!is.na(dat$thrush) & dat$thrush==1,1,0)

# create a new column for status
dat$Hstatus <- NA
dat$Hstatus[dat$hiveu==1] <- 1
dat$Hstatus[dat$hivuu==1] <- 2

# nasal flaring
dat$nasalfl <- dat$nasalfl
dat$nasalfl[dat$nasalfl==2] <- 0
#dat$nasalfl <- ifelse(!is.na(dat$nasalfl) & dat$nasalfl==1,1,0)

# poor drinking or feeding
dat$drinking <- NA
dat$drinking[!is.na(dat$thirst)] <- ifelse(dat$thirst[!is.na(dat$thirst)]==3,1,0)

# lethargic
dat$lethargic <- ifelse(!is.na(dat$genapp) & dat$genapp==3,1,0)  

# Blantyre comma score (BSC)
dat$BCS <- NA
dat$BCS[!is.na(dat$bcs)] <- ifelse(dat$bcs[!is.na(dat$bcs)] <=2,1,0)


# radialogically confirmed pneumonia
dat$pnemo <- NA
dat$pnemo[dat$ccc=="Radiological confirmed Pneumonia"] <- 1
dat$pnemo[dat$ccc=="0"] <- 0

# bronchial breathing
dat$brbreathing <- ifelse(!is.na(dat$brbr) & dat$brbr!=1,1,0)

# previous visit to the hospital
dat$prvHosp <- ifelse(!is.na(dat$phc) & dat$phc==0|dat$phc==1,1,0) # compared to a single visit

#dat$Hosp <- ifelse(!is.na(dat$phc) & dat$phc==0,0,1)

# crackles
dat$cRackles <- ifelse(!is.na(dat$craks) & dat$craks !=1,1,0)

# reduced air entry
dat$redAir <- ifelse(!is.na(dat$redae) & dat$redae !=1,1,0)

# skin pinching
dat$skin.pinch <- dat$sknpnch
dat$skin.pinch[dat$sknpnch==3] <- NA

# rash
dat$skinrash <- dat$rash
dat$skinrash[dat$rash==9] <- NA

dat$APN <- dat$apn
dat$APN[dat$apn==9] <- NA

# remove missing values from some variables
dat$vomiting <- dat$vom
dat$vomiting[dat$vomiting==9] <- NA
#dat$vomiting < ifelse(!is.na(dat$vomiting) & dat$vomiting == 1,1,0)

dat$diarrhoea <- dat$diarr
dat$diarrhoea[dat$diarrhoea==9] <- NA

dat$convulsions <- dat$conv
dat$convulsions[dat$convulsions==9] <- NA

# age of mother
dat$Agemother <- as.numeric(dat$mothage)
# age group of mother
dat$mothAgeGroup <- ifelse(!is.na(dat$mothage) & dat$mothage < 20,1,2) # mother age

# define disease severity

# first define outcome variable
dat$mortality <- NA
dat$mortality[dat$dispos=="3"] <- 1
dat$mortality[dat$dispos!="3"] <- 0

#dat$diseaseSeverityA <- 0
#dat$diseaseSeverityA[!is.na(dat$mortality) & dat$mortality==1]<-1
#dat$diseaseSeverityA[(!is.na(dat$cpap) & dat$cpap==1) & ((!is.na(dat$sats) & dat$sats<90) | (!is.na(dat$oxygen) & dat$oxygen==1) | (!is.na(dat$ivfluid) & dat$ivfluid==1) | (!is.na(dat$hosplen) & dat$hosplen>=3))]<-1
#dat$diseaseSeverityA[((!is.na(dat$sats) & dat$sats<90) | (!is.na(dat$oxygen) & dat$oxygen==1)) & ( ((!is.na(dat$ivfluid) & dat$ivfluid==1) & (!is.na(dat$hosplen) & dat$hosplen>=3)) |
#                                                                                                     ((!is.na(dat$ivfluid) & dat$ivfluid==1) & (!is.na(dat$muac) & dat$muac<11.5)) | 
#                                                                                                     ((!is.na(dat$muac) & dat$muac<11.5) & (!is.na(dat$hosplen) & dat$hosplen>=3))  )]<-1 

# MUAC left out for now (?)

dat$diseaseSeverityA <- 0
dat$diseaseSeverityA[!is.na(dat$mortality) & dat$mortality==1]<-1
dat$diseaseSeverityA[(!is.na(dat$cpap) & dat$cpap==1) & ((!is.na(dat$sats) & dat$sats<90) | (!is.na(dat$oxygen) & dat$oxygen==1) | (!is.na(dat$ivfluid) & dat$ivfluid==1) | (!is.na(dat$hosplen) & dat$hosplen>=3))]<-1
dat$diseaseSeverityA[((!is.na(dat$sats) & dat$sats<90) | (!is.na(dat$oxygen) & dat$oxygen==1)) & ( (!is.na(dat$ivfluid) & dat$ivfluid==1) & (!is.na(dat$hosplen) & dat$hosplen>=3) )]<-1 # MUAC left out for now (?)

### compute risk scores

# HIV uninfected
riscHIVn<-rep(0,nrow(dat))
riscHIVn[is.na(dat$sats)]<-NA
riscHIVn[!is.na(dat$sats) & dat$sats<=90]<-riscHIVn[!is.na(dat$sats) & dat$sats<=90]+3
riscHIVn[!is.na(dat$sats) & dat$sats>90 & is.na(dat$icid)]<-NA
riscHIVn[!is.na(dat$sats) & dat$sats>90 & !is.na(dat$icid) & dat$icid==1]<-riscHIVn[!is.na(dat$sats) & dat$sats>90 & !is.na(dat$icid) & dat$icid==1]+2
riscHIVn[is.na(dat$whez)]<-NA
riscHIVn[!is.na(dat$whez) & dat$whez==1]<-riscHIVn[!is.na(dat$whez) & dat$whez==1]-2

# feeding
riscHIVn[is.na(dat$thirst)]<-NA # using poor drinking / inability to drink as proxy for poor 
riscHIVn[!is.na(dat$thirst) & dat$thirst==3]<-riscHIVn[!is.na(dat$thirst) & dat$thirst==3]+1

# growth standards
riscHIVn[!is.na(dat$zweightAge) & dat$zweightAge<=(-3)]<-riscHIVn[!is.na(dat$zweightAge) & dat$zweightAge<=(-3)]+2
riscHIVn[!is.na(dat$zweightAge) & dat$zweightAge>(-3) & dat$zweightAge<=(-2)]<-riscHIVn[!is.na(dat$zweightAge) & dat$zweightAge>(-3) & dat$zweightAge<=(-2)]+1


## add the risk scores to the data

dat$riscHIVn <- riscHIVn

# define categories
dat$riscHIVcat <- NA
dat$riscHIVcat[!is.na(dat$riscHIVn)  & dat$riscHIVn <=3] <- 1
dat$riscHIVcat[!is.na(dat$riscHIVn) & dat$riscHIVn >3] <- 2


# mRisc scores
mRisc <- rep(0,nrow(dat))
mRisc[!is.na(dat$genapp) & dat$genapp==3 & dat$AgeWhol<=59] <- mRisc[!is.na(dat$genapp) & dat$genapp==3 & dat$AgeWhol<=59]+1
mRisc[!is.na(dat$thirst)  & !is.na(dat$AgeWhol) & dat$thirst==3 & dat$AgeWhol<=59] <- mRisc[!is.na(dat$thirst) & !is.na(dat$AgeWhol) & dat$thirst==3 & dat$AgeWhol<=59]+1
mRisc[!is.na(dat$icid) & !is.na(dat$AgeWhol)& dat$icid==1 & dat$AgeWhol<=59] <- mRisc[!is.na(dat$icid)& !is.na(dat$AgeWhol) & dat$icid==1 & dat$AgeWhol<=59]+1
mRisc[!is.na(dat$genapp) & !is.na(dat$AgeWhol)& dat$genapp==1 & dat$AgeWhol<=59] <- mRisc[!is.na(dat$genapp) & !is.na(dat$AgeWhol) & dat$genapp==1 & dat$AgeWhol<=59]+2
mRisc[!is.na(dat$zweightAge) & !is.na(dat$AgeWhol)& dat$zweightAge <= (-2) & dat$AgeWhol<=59] <- mRisc[!is.na(dat$zweightAge) & !is.na(dat$AgeWhol) & dat$zweightAge <= (-2) & dat$AgeWhol<=59]+1 # use zweightAge

dat$mRisc <- mRisc
dat$mRisc_cat <- NA
dat$mRisc_cat[!is.na(dat$mRisc) & dat$mRisc <=3] <- 1
dat$mRisc_cat[!is.na(dat$mRisc) & dat$mRisc > 3] <- 2

### number that were followed up
sum(!is.na(dat$fudt))

### export cleaned data for further analysis ###
#write.csv(dat,"data/heu_cleandata07102020.csv",row.names = FALSE)
write.csv(dat,"data/heu_data_clean_12022021.csv",row.names = FALSE)


