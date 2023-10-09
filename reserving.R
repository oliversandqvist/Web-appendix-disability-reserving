library(tidyverse)

LECDK19 <- load(url(URLencode("https://github.com/oliversandqvist/MultistateDelayAdjudication/Data application/raw/main/LECDK19.Rdata")))

date <- as.POSIXct("2019-07-31 23:59:59 UTC", tz="UTC")
studyStart <- "2015-01-31" 
dayYear <- 365.25
time0 <- as.numeric(difftime(date,studyStart, units="days"))/dayYear

LECDK19df <- LECDK19
LECDK19df <- LECDK19df[-c(1,2,3)]
LECDK19df$disab <- LECDK19df$disab %>% filter(dateStart == date)
LECDK19df$reac <- LECDK19df$reac %>% filter(dateStart == date)

set.seed(98765)

n <- 100
RBNSrdfOrig <- LECDK19df$reac[sample(1:nrow(LECDK19df$reac),n),]
RBNSidfOrig <- LECDK19df$disab %>% filter(!is.na(durDisab) & adjState != 3)
RBNSidfOrig <- RBNSidfOrig[sample(1:nrow(RBNSidfOrig),n,replace=TRUE),]
CBNRdfOrig <- LECDK19df$disab %>% filter(is.na(durDisab)) 
CBNRdfOrig <- CBNRdfOrig[sample(1:nrow(CBNRdfOrig),n),]


##ibnr
delayCoef <- c(2.2253524712,1.0531730783,0.3556644822,0.1132067522,0.0008557776)

I <- function(s,t,gender,age){
  1-(1-exp(-(delayCoef[1]*(t-s))^delayCoef[2]))^(delayCoef[3]*(gender=="M")+delayCoef[4]*(gender=="F")+delayCoef[5]*age)
}

##adjudication
tmax=8
nstep=100
rk4 <- function(df, a, b, f0, n) {
  
  h = (b-a)/n
  
  f = f0
  for (i in 1:n) {
    s = a + h * (i-1)
    k1 = h * df(s,f)
    k2 = h * df(s+0.5*h, f+0.5*k1)
    k3 = h * df(s+0.5*h, f+0.5*k2)
    k4 = h * df(s+h, f+k3)
    
    f = f+1/6*(k1+2*k2+2*k3+k4)
  }
  return(f)
}

diffAdjReac = function(t, x, mu12, mu21, mu13, mu14, mu24) {
  d1 = -x[1]*(mu12(t)+mu13(t)+mu14(t))+x[2]*mu21(t)
  d2 = -x[2]*(mu21(t)+mu24(t))+x[1]*mu12(t)
  d3 = x[1]*mu13(t)
  d4= x[1]*mu14(5)+x[2]*mu24(t)
  return( c(d1, d2, d3, d4) )
}

#note: model space is extended so "rejectedBefore" is captured in statespace - it is now state 4 and dead is state 5
diffAdjDisab = function(t, x, mu13, mu12, mu24, mu42, mu43, mu15, mu45, mu25) {
  d1 = -x[1]*(mu13(t)+mu12(t)+mu15(t))
  d2 = -x[2]*(mu24(t)+mu25(t))+x[4]*mu42(t)
  d3 = x[1]*mu13(t)+x[4]*mu43(t)
  d4 = -x[4]*(mu43(t)+mu42(t)+mu45(t))+x[2]*mu24(t)
  d5 = x[4]*mu45(t)+x[1]*mu15(t)+x[2]*mu25(t)
  return( c(d1, d2, d3, d4, d5) )
}


reacRejectCoef <- c(0.005007513,1.183425789,0.972683981,-0.041822943,-0.065379706)
reacReapplyCoef <- c(6.838564e-05,-5.328340e-02,-2.082540e-01,-1.245490e+00,2.055961e-01)
reacAwardCoef <- c(0.01546058,0.77122201,1.29111478,-0.48998518,-0.19835554 )
reacDeadIBNRCoef <- c(0.02319754,-21.09900823,-5.94234612,0.14767608,-0.46433496)

absProbReac <- function(y0,age0,gender,durDisab0,durReac0){
  mu12 <- function(t){exp(reacRejectCoef %*% c(age0+t,(gender=="F"),(gender=="M"),durReac0+t,durDisab0+t))}
  mu21 <- function(t){exp(reacReapplyCoef %*% c(age0+t,(gender=="F"),(gender=="M"),durReac0+t,durDisab0+t))}
  mu13 <- function(t){exp(reacAwardCoef %*% c(age0+t,(gender=="F"),(gender=="M"),durReac0+t,durDisab0+t))}
  mu14 <- function(t){0}
  mu24 <- function(t){exp(reacDeadIBNRCoef %*% c(age0+t,(gender=="F"),(gender=="M"),durReac0+t,durDisab0+t))}
  
  distrib <- rk4(function(t, x) diffAdjReac(t, x, mu12=mu12, mu21=mu21, mu13=mu13, mu14=mu14, mu24=mu24), a=0, b=tmax, f0=y0, nstep)
  return(distrib[3])
} 


disabAwardCoef <- c(0.0007125963,1.1868871662,1.0592085304,0.4697592786,-0.2344138034,-0.1175834324)
disabRejecCoef <- c(-0.004817674,-0.806449160,-0.799717698,-0.249845173,0.077125150,0.982297988)
disabReapplyCoef<- c(0.002596526,1.084541659,1.159223225,-0.090002484,-0.255215380)
deadRBNSCoef <- c(0.01270125,-5.28526166,-4.84594645,-0.53219385,0.16933611,-14.67544111)

absProbDisab <- function(y0,age0,gender,durDisab0,durReport0){
  mu13 <- function(t){exp(disabAwardCoef %*% c(age0+t,(gender=="F"),(gender=="M"),durReport0+t,durDisab0+t,0))}
  mu12 <- function(t){exp(disabRejecCoef %*% c(age0+t,(gender=="F"),(gender=="M"),durReport0+t,durDisab0+t,0))}
  mu24 <- function(t){exp(disabReapplyCoef %*% c(age0+t,(gender=="F"),(gender=="M"),durReport0+t,durDisab0+t))}
  mu42 <- function(t){exp(disabRejecCoef %*% c(age0+t,(gender=="F"),(gender=="M"),durReport0+t,durDisab0+t,1))}
  mu43 <- function(t){exp(disabAwardCoef %*% c(age0+t,(gender=="F"),(gender=="M"),durReport0+t,durDisab0+t,1))} 
  mu15 <- function(t){exp(deadRBNSCoef %*% c(age0+t,(gender=="F"),(gender=="M"),durReport0+t,durDisab0+t,0))}
  mu45 <-function(t){exp(deadRBNSCoef %*% c(age0+t,(gender=="F"),(gender=="M"),durReport0+t,durDisab0+t,1))}
  mu25 <- function(t){0}
  
  distrib <- rk4(function(t, x) diffAdjDisab(t, x, mu13=mu13, mu12=mu12, mu24=mu24, mu42=mu42, mu43=mu43, mu15=mu15, mu45=mu45, mu25=mu25), a=0, b=tmax, f0=y0, nstep)
  return(distrib[3])
} 

adjStateDisab <- function(adjState,rejectedBefore){c((adjState==1 & rejectedBefore==0),(adjState==2),(adjState==3),(adjState==1 & rejectedBefore==1),(adjState==5))}
adjStateReac <- function(adjState) c((adjState==1),(adjState==2),(adjState==3),(adjState==4))

##hazards

#note: time argument refers to time elapsed after time0 (the latter being the time where the reserves are calculated)
disabCoef <- c(0.02276197,-8.65900843,-7.46187930,0.30431135)
disabMu <- function(age0,gender,time){
  return( exp(disabCoef %*% c(age0+time,(gender=="F"),(gender=="M"),time0)) )
}

reacCoef <- c(-0.01217323,0.75640213,0.33436334,-0.12466354,-1.04464466)
reacMu <- function(age0,gender,time,duration0){
  return( exp(reacCoef %*% c(age0+time,(gender=="F"),(gender=="M"),time0,duration0+time)) )
}

durCutoff <- 5
disabDeadCoef <- c(0.09,-6.8,-6.4,-0.25)
disabDeadMu <- function(age0,gender,time, duration0){
  return( exp(disabDeadCoef %*% c(age0+time,(gender=="F"),(gender=="M"),min(duration0+time,durCutoff))) )
}

activeDeadCoef <- c(0.09,-9.8,-9.5)
activeDeadMu <- function(age0,gender,time){
  return(exp(activeDeadCoef %*% c(age0+time,(gender=="F"),(gender=="M"))))
}

reacDeadMu <- function(age0,gender,time){
  return(activeDeadMu(age0,gender,time))
}


##reserve parameters
r <- 0.02
ageRetire <- 67
coveragePeriod <- 3
epsilon <- 10^(-7)


##reserve functions
diffValidTimeModel = function(t, x, mu12, mu14, mu23, mu24, mu34) {
  d1 = -x[1]*(mu12(t)+mu14(t))
  d2 = -x[2]*(mu23(t)+mu24(t))+x[1]*mu12(t)
  d3 = -x[3]*mu34(t)+x[2]*mu23(t)
  d4 = x[1]*mu14(t)+x[2]*mu24(t)+x[3]*mu34(t)
  return( c(d1, d2, d3, d4) )
}

probFromDisabOrReac <- function(y0,fromTime,toTime,toState,age0,gender,duration0){
  mu12 <- function(t){disabMu(age0,gender,t)}
  mu14 <- function(t){activeDeadMu(age0,gender,t)}
  mu23 <- function(t){reacMu(age0,gender,t,duration0)}
  mu24 <- function(t){disabDeadMu(age0,gender,t,duration0)}
  mu34 <- function(t){reacDeadMu(age0,gender,t)}
  
  distrib <- rk4(function(t, x) diffValidTimeModel(t, x, mu12=mu12, mu14=mu14, mu23=mu23, mu24=mu24, mu34=mu34), a=fromTime, b=toTime, f0=y0, nstep)
  return(distrib[toState])
} 

diffThieleViGt <- function(t,x,age0,mu23,mu24){
  d1 = (r+mu23(t)+mu24(t))*x[1]-(age0+t <= ageRetire)
  return(c(d1))
}
  
thieleViGt <- function(age0,gender,duration0,tNoDeath){
  mu23 <- function(t){
    factor <- ifelse(t >= tNoDeath,1,(1-probFromDisabOrReac(c(0,0,1,0),t,tNoDeath,4,age0,gender,0))/(1-probFromDisabOrReac(c(0,1,0,0),t,tNoDeath,4,age0,gender,duration0)) )
    reacMu(age0,gender,t,duration0)*factor
    } 
  mu24 <- function(t){
    factor <- ifelse(t >= tNoDeath,1, 0)
    disabDeadMu(age0,gender,t,duration0)*factor
    } 
  fromTime <- ageRetire-age0
  
  V <- rk4(function(t, x) diffThieleViGt(t, x, age0=age0, mu23=mu23, mu24=mu24), a=fromTime, b=0+epsilon, f0=0, nstep)
  return(V[1])
} 


diffThieleVa <- function(t,x,age0,gender,mu12,mu14){
  d1 = (r+mu12(t)+mu14(t))*x[1]-ifelse(t<=coveragePeriod,mu12(t)*thieleViGt(age0+t,gender,0,0),0)
  return(c(d1))
}

thieleVa <- function(age0,gender){
  mu12 <- function(t){disabMu(age0,gender,t)} 
  mu14 <- function(t){activeDeadMu(age0,gender,t)} 
  fromTime <- min(coveragePeriod,ageRetire-age0)
  
  V <- rk4(function(t, x) diffThieleVa(t, x, age0=age0, gender=gender, mu12=mu12, mu14=mu14), a=fromTime, b=0+epsilon, f0=0, nstep)
  return(V[1])
} 

activeDeadMuInt <- function(age0,gender,time){
  num <- activeDeadMu(age0,gender,time) 
  denom <- activeDeadCoef[1]
  return(num/denom)    
}

disabMuInt <- function(age0,gender,time){
  num <- disabMu(age0,gender,time)
  denom <- disabCoef[1]
  return(num/denom)
}

paa <- function(age0,gender,time){
  disabTerm <- disabMuInt(age0,gender,time)-disabMuInt(age0,gender,0)
  activeDeadTerm <- activeDeadMuInt(age0,gender,time)-activeDeadMuInt(age0,gender,0)
  return( exp(-(disabTerm+activeDeadTerm)) )
}
paaVec <- Vectorize(paa)
 
activeDeadMuVec <- Vectorize(activeDeadMu)
disabMuVec <- Vectorize(disabMu)
IVec <- Vectorize(I)
paaVec <- Vectorize(paa)
pCBNR <- function(age0,gender){
  beforeTerm <- integrate(f = function(time){IVec(time,0,gender,age0+time)*paaVec(age0-time0,gender,time0+time)*disabMuVec(age0-time0,gender,time0+time)}, lower=-time0,upper=0)$value
  afterTerm <- integrate(f = function(time){paaVec(age0-time0,gender,time0+time)*disabMuVec(age0-time0,gender,time0+time)}, lower=0,upper=Inf)$value
  infTerm <- integrate(f = function(time){paaVec(age0-time0,gender,time0+time)*activeDeadMuVec(age0-time0,gender,time0+time)}, lower=0,upper=1000)$value
  
  return(beforeTerm+afterTerm+infTerm)
}

VaCBNRFactor <- function(age0,gender){
  integ <- 1/(pCBNR(age0,gender))*integrate(f = function(time){IVec(time,0,gender,age0+time)*disabMuVec(age0-time0,gender,time0+time)*paaVec(age0-time0,gender,time0+time)}, lower=-time0,upper=0)$value
  return(1-integ)
}

thieleViGtVec <- Vectorize(thieleViGt)
ViCBNRTerm <- function(age0,gender){
  integ <- 1/(pCBNR(age0,gender))*integrate(f = function(time){exp(r*(-time))*thieleViGtVec(age0+time,gender,0,-time)*IVec(time,0,gender,age0+time)*disabMuVec(age0-time0,gender,time0+time)*paaVec(age0-time0,gender,time0+time)}, lower=-time0,upper=0)$value
  return(integ)
}

##CBNR reserve (those who have never reported a disability)

CBNRdf <- CBNRdfOrig 

CBNRdf <- CBNRdf %>%
  rowwise() %>%
  mutate(VCBNRVaTerm = ifelse(age < 67,thieleVa(age,gender),0)*VaCBNRFactor(age,gender),
         VCBNRViTerm = ViCBNRTerm(age,gender),
         VCBNRNaive = ifelse(age < 67,thieleVa(age,gender),0)  ) %>% 
  ungroup()

sum(CBNRdf$VCBNRVaTerm+CBNRdf$VCBNRViTerm) #8.32
sum(CBNRdf$VCBNRNaive) #7.73

##RBNSi reserve (those who have reported a disability, but have not recieved any benefits yet)
RBNSidf <- RBNSidfOrig %>% 
  mutate(tMinusGt = durDisab,
         ageGt = age-tMinusGt,
         Wt = 0) %>%
  rowwise() %>%
  mutate(awardProb = absProbDisab(adjStateDisab(adjState,rejectedBefore),ageGt,gender,durDisab,durDisabReport)) %>%
  ungroup() %>%
  mutate(awardProb = ifelse(is.na(awardProb),1,awardProb))

RBNSidf <- RBNSidf %>%
  rowwise() %>%
  mutate(VRBNSi = exp(tMinusGt*r)*thieleViGt(ageGt,gender,Wt,tMinusGt)*awardProb,
         VRBNSiNaive = ifelse(age < 67,thieleVa(age,gender),0)  ) %>% 
  ungroup()


sum(RBNSidf$VRBNSi) #501.69
sum(RBNSidf$VRBNSiNaive) #5.23

##RBNSr reserve (those currently receiving running payments or have been reactivated)

RBNSrdf <- RBNSrdfOrig %>% 
  mutate(tMinusGt = ifelse(is.na(durReac),0,durReac),
                              ageGt = age-tMinusGt,
                              Wt = durDisab-tMinusGt) %>%
  rowwise() %>%
  mutate(awardProb = absProbReac(adjStateReac(adjState),ageGt,gender,durDisab,durReac)) %>%
  ungroup() %>%
  mutate(awardProb = ifelse(is.na(awardProb),1,awardProb))

RBNSrdf <- RBNSrdf %>%
  rowwise() %>%
  mutate(VRBNSr = exp(tMinusGt*r)*thieleViGt(ageGt,gender,Wt,tMinusGt)*awardProb,
         VRBNSrNaive =thieleViGt(age,gender,durDisab,0)*ifelse(awardProb < 1,0,1) ) %>%
  ungroup()

sum(RBNSrdf$VRBNSr) #961.68
sum(RBNSrdf$VRBNSrNaive) #945.35

##save results:
#dataResult <- list()
#dataResult$CBNR <- CBNRdf
#dataResult$RBNSi <- RBNSidf
#dataResult$RBNSr <- RBNSrdf
#save(dataResult, file = "Results/dataResult.Rda") 


##Full population extrapolation
RBNSrdfFull <- LECDK19df$reac
RBNSidfFull <- LECDK19df$disab %>% filter(!is.na(durDisab) & adjState != 3)
CBNRdfFull <- LECDK19df$disab %>% filter(is.na(durDisab)) 

CBNRAvg <- sum(CBNRdf$VCBNRVaTerm+CBNRdf$VCBNRViTerm)/n 
CBNRNaiveAvg <- sum(CBNRdf$VCBNRNaive)/n
RBNSiAvg <- sum(RBNSidf$VRBNSi)/n
RBNSiNaiveAvg <- sum(RBNSidf$VRBNSiNaive)/n
RBNSrAvg <- sum(RBNSrdf$VRBNSr)/n 
RBNSrNaiveAvg <- sum(RBNSrdf$VRBNSrNaive)/n

nCBNR <- nrow(CBNRdfFull)
nRBNSi <- nrow(RBNSidfFull)
nRBNSr <- nrow(RBNSrdfFull)

dfPortfolio <- data.frame(Reserve = c(CBNRAvg*nCBNR,CBNRNaiveAvg*nCBNR,RBNSiAvg*nRBNSi,RBNSiNaiveAvg*nRBNSi,RBNSrAvg*nRBNSr,RBNSrNaiveAvg*nRBNSr) 
                          ,Method = c("Proposed","Naive","Proposed","Naive","Proposed","Naive") 
                          ,Category =c("CBNR","CBNR","RBNSi","RBNSi","RBNSr","RBNSr") )


png("Figures/Barchart.png", width = 4.5, height = 4, units = 'in', res = 500) 

ggplot(dfPortfolio, aes(x = Method, y = Reserve, fill = Category, label = round(Reserve,0) )) +
  geom_bar(stat = "identity") +
  geom_text(size = 4.5, position = position_stack(vjust = 0.5),color="white") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank())  +
  scale_fill_grey(start=0.7, end=0.3)

dev.off()