### code for data analysis ###

rm(list = ls())

library(tidyverse)
library(epiR)
library(arsenal)
#library(sjPlot)
#library(glmnet)
#library(ROCit)
#library(Amelia)

#dat <- read_csv("data/heu_cleandata07102020.csv")
dat <- read_csv("data/heu_data_clean_12022021.csv")
##### create baseline table ####

### table controls ###
tabcontrols  <- tableby.control(test=TRUE, total=TRUE,
                               numeric.test="kwt", cat.test="fe",  
                               numeric.stats=c("N","medianq1q3"),
                               cat.stats=c("N","countpct"),
                               digits = 2,
                               digits.pct = 2,
                               digits.p = 2,
                               stats.labels=list(
                                 N='Count',
                                 meansd="Mean(SD)",
                                 medianq1q3='Median(Q1,Q3)'),
                               simulate.p.value = T)
# tab controls for z scores
### table controls ###
#ctrls  <- tableby.control(test=TRUE, total=TRUE,
#                                numeric.test="kwt", cat.test="fe",  
#                                numeric.stats=c("N","mean","sd"),
#                                cat.stats=c("N","countpct"),
#                                digits = 2,
#                                digits.pct = 2,
#                                digits.p = 2,
#                                stats.labels=list(
#                                  N='Count',
#                                  meansd="Mean(SD)",
#                                  medianq1q3='Median(Q1,Q3)'),
#                                simulate.p.value = T)


tablabels <- list(
  AgeWhol = "Age",
  age_grp = "Age group",
  sex="Sex",
  bwt="Birth weight(Kg)",
  LBW="Low birth weight",
  mothAgeGroup="Mother's age group",
  hh_size1="Household size", # based on the WHO cut-off point
  antib.use="Antibiotic use",
  birthVacc="Birth vaccinations",
  sixwksVacc="6 week vaccinations",
  tenwksVacc="10 week vaccinations",
  fourteenwksVacc="14 week vaccination",
  wforHght="Weight for Height",
  muac="MUAC",
  muac.score="MUAC (<11.5)",
  tachycardia="Age-adjusted tachycardia",
  hfever="Fever",
  tachpenia="Age-adjusted tachypnea",
  sats="Oxygen saturation",
  hypo="Hypoxemia (<90%)",
  drinking="Ability to drink",
  vomiting="Vomiting",
  diarrhoea="Diarrhoea",
  convulsions="Convulsion",
  skinrash="Rash",
  oedema="Oedema",
  dangerSigns="Danger signs",
  whez="Wheezing",
  nasalfl="nasal flaring",
  strid="Stridor",
  cRackles="Crackles",
  icid="Lower chest wall indrawing",
  redAir="Reduced air entry",
  thrush="Oral thrush/sores",
  club="Finger clubbing",
  bcs="Blantyre Comma Score (BCS)",
  temp="Body temperature",
  hosplen="Hospital stay duration",
  mortality="Mortality",
  brbreathing="Bronchial breathing",
  skin.pinch="SKin pinching",
  snkeye="Sunken eyes",
  apn="Apnoea",
  pallor="Pallor",
  nckstiff="Neck stiffness",
  lymphnd="Lymphadenopathy",
  hepgly="Hepatomegaly",
  splgly="Splenomegaly",
  cardabn="Cardiovascular abnormality",
  earabn="Ear abnormality",
  neuroabn="Neurplogical abnormality",
  cxrdne="Chest x-ray done",
  bcdne="Blood culture done",
  csfdne="Lumbar puncture done",
  mpsdne="mps done",
  ivfluid="Received IVF",
  oxygen="Received oxygen",
  cpap="Placed on CPAP",
  mRisc_cat="mRIsC scores",
  riscHIVncat = 'RISC',
  diseaseSeverityA="Severe disease"
)

# replace odd z-scores with new values

tab1 <- tableby(Hstatus ~ AgeWhol
                + factor(age_grp)
                + factor(sex) 
                + bwt
                + factor(LBW) +
                + Agemother
                + factor(mothAgeGroup)
                + factor(hh_size1)
                + factor(antib.use) + as.factor(epiVacHue)
                + zweightHeight
                + factor(zwforHght) # replaced with a new value
                + muac
                + factor(muac.scorev2)
                + factor(tachpenia)
                + factor(tachycardia)
                + temp
                + factor(hfever)
                + sats
                + factor(hypo) 
                + factor(drinking) 
                + factor(vomiting) 
                + factor(diarrhoea) 
                + as.factor(convulsions) 
                + factor(skinrash) 
                + factor(oedema) + factor(dangerSigns)
                + bcs + factor(bcs) 
                + factor(skin.pinch) + factor(snkeye) + factor(APN) + factor(pallor) + factor(club) + factor(jaund) 
                + factor(oedema) + factor(nasalfl) + factor(icid) + factor(strid) + factor(cRackles) + factor(whez) 
                + factor(brbreathing) + as.factor(redAir) + factor(thrush) + factor(nckstiff) + factor(lymphnd) + factor(hepgly)
                + factor(splgly) + factor(cardabn) + factor(earabn) + factor(neuroabn)
                + factor(cxrdne) + factor(bcdne) + factor(csfdne) + factor(mpsdne) + factor(ivfluid) + factor(oxygen) + factor(cpap)
                + hosplen + factor(mortality) 
                + factor(mRisc_cat) + factor(riscHIVcat) +
                + factor(diseaseSeverityA),data=dat, control = tabcontrols)

summary(tab1,text = TRUE,labelTranslations = tablabels,pfootnote = TRUE)
#write2word(tab1,file = "baselineTable11202020.doc")
write2word(tab1,file = "baselineTable23012021v3.doc")


## risk scores plots ##
tiff("images/risk_scores_revised.tif",width = (25*0.39),height = (25*0.39),units = "in",res = 350,compression = "lzw")
par(mfrow=c(1,1),cex.lab=1.5,cex.axis=1.5,mar=c(4,5,2,1))
barplot(table(dat$riscHIVn),col = "#f5b30c",ylab = "Frequency",xlab = "RISC score")
#abline(v=6.70,lty=2,col=2,lwd=2)
#barplot(table(dat$mRisc),col = "#f5b30c",ylab = "Frequency",xlab = "mRISC score")
#abline(v=4.33,lty=2,col=2,lwd=2)
#mtext("A",side = 2,line = 0,at=0,400,las=2,cex = 2)
dev.off()

## data subset for diagnostic tests ##
epiHeu <- dat[which(dat$hiveu==1),]
epiHuu <- dat[which(dat$hivuu==1),]


## function to calculate diagnostic tests for the HUU group

calculate_diagnostics <- function(outcome,exposure){
  xtab <- table(factor(epiHuu[[exposure]],levels = c("1","0")),factor(epiHuu[[outcome]],levels = c("1","0")))
  xtabTest <- epiR::epi.tests(xtab)
  res <- data.frame(Est=c(xtabTest$elements$ppv,
                          xtabTest$elements$npv,
                          xtabTest$elements$lrpos,
                          youden=xtabTest$elements$youden$est),
                    Lower=c(xtabTest$elements$ppv.low,
                            xtabTest$elements$npv.low,
                            xtabTest$elements$lrpos.low,
                            xtabTest$elements$youden$lower),
                    Upper=c(xtabTest$elements$ppv.up,
                            xtabTest$elements$npv.up,
                            xtabTest$elements$lrpos.up,
                            xtabTest$elements$youden$upper))
  rownames(res) <- c("ppv","npv","lr","youden")
  return(res)
}

calculate_diagnostics("diseaseSeverityA","hfever")
calculate_diagnostics("diseaseseverityA","LBW")
calculate_diagnostics("diseaseSeverityA","dangerSigns")
calculate_diagnostics("diseaseSeverityA","rash")
#calculate_diagnostics("diseaseSeverityA","lethargic")
calculate_diagnostics("diseaseSeverityA","nasalfl")
calculate_diagnostics("diseaseSeverityA","whez")
calculate_diagnostics("diseaseSeverityA","icid")
calculate_diagnostics("diseaseSeverityA","thrush")
calculate_diagnostics("diseaseSeverityA","strid")
#calculate_diagnostics("diseaseSeverityA","snkeye")
#calculate_diagnostics("diseaseSeverityA","hhsize")
calculate_diagnostics("diseaseSeverityA","antib.use") 
calculate_diagnostics("diseaseSeverityA","epiVacHue")
calculate_diagnostics("diseaseSeverityA","tachpenia")
#calculate_diagnostics("diseaseSeverityA","tachpenia2")

calculate_diagnostics("diseaseSeverityA","muac.score")
calculate_diagnostics("diseaseSeverityA","prvHosp")
calculate_diagnostics("diseaseSeverityA","pnemo")
calculate_diagnostics("diseaseSeverityA","brbreathing")
calculate_diagnostics("diseaseSeverityA","vom")
calculate_diagnostics('diseaseSeverityA',"diarr")
calculate_diagnostics("diseaseSeverityA","conv")
calculate_diagnostics("diseaseSeverityA","cRackles")
calculate_diagnostics("diseaseSeverityA","redAir")
calculate_diagnostics("diseaseSeverityA","hepgly")
calculate_diagnostics("diseaseSeverityA","drinking")
calculate_diagnostics("diseaseSeverityA","hypo")
calculate_diagnostics("diseaseSeverityA","cpap")
calculate_diagnostics("diseaseSeverityA","oxygen")
calculate_diagnostics("diseaseSeverityA","hh_size1")
calculate_diagnostics("diseaseSeverityA","sex")
calculate_diagnostics("diseaseSeverityA","stridor")
calculate_diagnostics("diseaseSeverityA","mRisc_cat2")

epiHuu$tachpenia2 <- ifelse(epiHuu$tachpenia==1,1,0)
epiHuu$mRisc_cat2 <- ifelse(epiHuu$mRisc_cat==2,1,0)

# comments from reviewers