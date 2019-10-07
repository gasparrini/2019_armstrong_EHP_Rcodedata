###############################################################################
# Updated version of the R code for the analysis in:
#
#   "The role of humidity in associations of high temperature with mortality:
#       a multicountry multicity study"
#   Ben Armstrong, ..., Antonio Gasparrini
#   Environmental Health Perspectives - 2019
#   http://www.ag-myresearch.com/2019_armstrong_ehp.html
#
# Update: 07 October 2019
# * an updated version of this code, compatible with future versions of the
#   software, is available at:
#   https://github.com/gasparrini/2019_armstrong_EHP_Rcodedata
###############################################################################

# Core R code (may be demonstrated using data_UK_SE.Rdata) - BEN ARMSTRONG
# CORE CODE FOR ANALYSIS OF ASSOCIATION OF RH ILLUSTRATED BY APPLICATION TO DATASET FOR UK SE REGION. 
# NOTES:
#  1/ THIS CODE IS SIMPLIFIED FROM THAT USED IN THE MULTI-LOCATION ANALYSIS TO MAKE MORE TRANSPARENT,
#       BUT THE MODELS FITTED ARE THE SAME.
#  2/ THE DATASET INCLUDES ALL DAYS FROM 29 MAY TO 30 SEPT TO ALLOW CALCULATION OF LAGGED VALUES UP TO LAG 3 
#       WHILE INCLUDING ALL OF JUNE-SEPT DAYS IN MODELS.
#       MODELS WITH CROSS-BASES EXCLUDE 29-31 MAY AUTOMATICALLY; FOR THE BASELINE MODEL THESE ARE EXCLUDED 
#       USING A SUBSET FUNCTION.
#  3/ THIS LOCATION WAS CHOSEN AS ONE FOR WHICH THERE WAS PERMISION TO MAKE PUBLIC, AND IN WHICH MOST RESULTS 
#       CONCERNING RH ARE BROADLY SIMILAR TO THE OVERALL MEANS NOTED IN THE PAPER

load( "data_UK_SER.Rdata" )
summary(data)
library(dlnm) ; library(splines);library(weathermetrics);library(humidity)

# INITIALISE DATA FRAME TO HOLD MODEL FIT STATISTICS TO FACILATE COMPARISON AND DEFINE FUNCTION TO EXTRACT THEM
setNames(modelfitstats <- data.frame(matrix(NA,0,6)),c("model","obs","dfmodel","dispersion","deviance","AIC" ))

extractmodelfit <- function(model,name) {
  modelfit <- data.frame(name,
    summary(model)$df.null+1,                               
    summary(model)$df[1],
    summary(model)$dispersio,       
    summary(model)$deviance,summary(model)$deviance + 2* summary(model)$df[1] * summary(model)$dispersion) 
  names(modelfit) <- c("model","obs","dfmodel","dispersion","deviance","AIC" )
  return(modelfit)
}

# BASELINE MODEL (1 IN PAPER) WITH TIME TERMS ONLY
baseformula <- "death ~  dayofweek + ns(date,df=round(length(unique(year))/10))+   ns(dayofsummer,df=4):factor(year) "
model <- glm(formula=baseformula,family=quasipoisson,data=subset(data,month>5 ), na.action="na.exclude") 
summary(model)
modelfitstats <- rbind(modelfitstats,extractmodelfit(model,"baseline(#1)" )  )

# CONSTRUCT CROSS-BASES FOR TMEAN AND RH
arglagboth <- list(fun="strata", breaks=c(0,1,2), intercept=F) 
cbtmean <- crossbasis(data$tmean,   lag=c(0,3),argvar=list(fun="ns", df=4 ), arglag=arglagboth, group=data$year )
cbRH   <-  crossbasis(data$RH,lag=c(0,3),argvar=list(fun="lin"),       arglag=arglagboth, group=data$year )

# PRIMARY MODEL INCLUDING TMEAN SPLINE AND RH LINEAR TERMS (MODEL 6 IN PAPER)
model <- glm(formula=paste(baseformula,"+ cbtmean + cbRH"),family=quasipoisson,data=data, na.action="na.exclude") 
summary(model)
modelfitstats <- rbind(modelfitstats,extractmodelfit(model,"primary(#6)" )  )

# DISPLAY ESTIMATE OF RR FOR BENCHMARK RH INCREMENT OF 23.4  
crosspred(cbRH, model, at=23.4)[c("allRRfit","allRRlow","allRRhigh")]

# ALTHOUGH PLOTTING IS UNNECESSARY FOR LINEAR MODELS, WE DO SO FOR ILLUSTRATION, SHOWING CUMULATIVE RR RELATIVE TO 
#   THE MEAN RH
plot( crosspred(cbRH, model, cen=mean(data$RH)),"overall")

# PLOT THE ESTIMATE LAG STRUCTURE FOR BENCHMARK 23.4% RH INCREMENT
plot(crosspred(cbRH, model, at=23.4), "slices", var=23.4, ci="bars", type="p", col=2, pch=19,
  main="Lag-response for 23.4%-unit increase in RH")

# GRAPH CUMULATIVE TEMPERATURE-MORTALITY CUMULTIVE ASSCOIATION (AFTER ADJUSTING FOR RH) 
plot( crosspred(cbtmean, model), "overall" , xlab="Temperature", ylab="RR" ,
  main="RR by daily mean temperature adjusting for RH")

# ADDING TERM FOR RHXTMEAN MODIFICATION: 
#  1: LINEAR-LINEAR (MODEL 12 IN PAPER)
# CONSTRUCT TMEANXRH INTERACTION TERM (LINXLIN)
data$tmeanXRH <- (data$tmean-mean(data$tmean))*(data$RH-mean(data$RH)) # centering tmean and RH to reduce co-linearity
cbtmeanXRH <-   crossbasis(data$tmeanXRH,lag=c(0,3),argvar=list(fun="lin"), arglag=arglagboth, group=data$year )
model <- glm(formula=paste(baseformula,"+ cbtmean + cbRH + cbtmeanXRH"),
  ,family=quasipoisson,data=data, na.action="na.exclude")
summary(model)
modelfitstats <- rbind(modelfitstats,extractmodelfit(model,"linXlin tempXRH int (#12)"  ) )

# ESTIMATED MODIFICATION OF RR/3.6  OVER 10 DEG 
crosspred(cbtmeanXRH, model, at=23.4*10)[c("allRRfit","allRRlow","allRRhigh")]

# ESTIMATE OF RR FOR RH INCREMENT OF 23.4% IF INTERACTION TERM IS INCLUDED
crosspred(cbRH, model, at=23.4)[c("allRRfit","allRRlow","allRRhigh")]

# 2: TMEANXRH INTERACTION TERM BY FITTING LINEAR RH EFFECT STRATIFIED B INTO 3 TMEAN GROUPS, HERE ILLUSTRATED 
#   FOR CUT-POINTS AT 22.6 DEG AND 25.2 DEG (MODEL 17 IN PAPER)
tgroupbreaks <- c(0,22.6,25.2)  
tmean3g <- dlnm:::strata(data$tmean,breaks=tgroupbreaks,intercept=F)
tmean3gXRH <- data$RH * tmean3g 
cbtmeang1XRH <- crossbasis(tmean3gXRH[,1],lag=c(0,3),argvar=list(fun="lin"), arglag=arglagboth, group=data$year )
cbtmeang2XRH <- crossbasis(tmean3gXRH[,2],lag=c(0,3),argvar=list(fun="lin"), arglag=arglagboth, group=data$year )
cbtmeang3XRH <- crossbasis(tmean3gXRH[,3],lag=c(0,3),argvar=list(fun="lin"), arglag=arglagboth, group=data$year )
model <- glm(formula=paste(baseformula,"+ cbtmean  + cbtmeang1XRH + cbtmeang2XRH + cbtmeang3XRH ")
  ,family=quasipoisson,data=data, na.action="na.exclude")
summary(model)
modelfitstats <- rbind(modelfitstats,extractmodelfit(model,"3-groupXlin tempXRH int (#17)"  ) )

# OBTAIN CUMULATIVE LAG0-3 RR FOR BENCHMARK RH INCREMENT BY TEMPERATURE GROUP
crosspred(cbtmeang1XRH, model, at=23.4)[c("allRRfit","allRRlow","allRRhigh")] 
crosspred(cbtmeang2XRH, model, at=23.4)[c("allRRfit","allRRlow","allRRhigh")] 
crosspred(cbtmeang3XRH, model, at=23.4)[c("allRRfit","allRRlow","allRRhigh")] 


# ILLUSTRATING REPEAT OF MAIN MODEL  USING SPECIFIC SPECIFIC NOT RELATIVE HUMIDITY (SENSITIVITY ANALYSIS IN PAPER)

# COMPUTE SPECIFIC HUMIDITY FROM RH AND TEMPERATURE USING PACKAGES WEATHERMETRICS AND HUMIDITY
data$DP <- humidity.to.dewpoint(rh=data$RH, t=data$tmean,temperature.metric="celsius")
data$SH <- SH(e=WVP1(data$DP, isK = F)*100 )*1000
summary(data$SH)

# CONSTRUCT CROSS-BASIS FOR SH
cbSH   <-  crossbasis(data$SH,lag=c(0,3),argvar=list(fun="lin"),  arglag=arglagboth, group=data$year )

# MODEL INCLUDING TMEAN SPLINE AND SH LINEAR TERMS (MODEL 6 IN PAPER, WITH SH REPLACING RH)
model <- glm(formula=paste(baseformula,"+ cbtmean + cbSH"),family=quasipoisson,data=data, na.action="na.exclude") 
modelfitstats <- rbind(modelfitstats,extractmodelfit(model,"primary with SH" )  )

# DISPLAY ESTIMATE OF RR FOR BENCHMARK SH INCREMENT OF 3.6 (99TH CENTILE ANOMALY, AS USED FOR RH)  
crosspred(cbSH, model, at=3.6)[c("allRRfit","allRRlow","allRRhigh")]

# PLOT THE ESTIMATED LAG STRUCTURE FOR BENCHMARK 3.6 UNITS SH INCREMENT
plot(crosspred(cbSH, model, at=3.6), "slices", var=3.6, ci="bars", type="p", col=2, pch=19,
  main="Lag-response for 3.6-unit increase in SH")

print(modelfitstats, digits=5)

# COMPUTE AND DISPLAY THE PARTIAL CORRELATION OF RH AND SH, PARTIALLING 
#      OUT THE TEMPERATURE SPLINE USING IN THE MORTALITY MODELS (TABLE S2)
library(ppcor)
pcor(cbind(ns(data$tmean,df=4),  data[c("RH","SH")]))$estimate[5:6,5:6]

  
