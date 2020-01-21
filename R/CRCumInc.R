#' Estimation of Cumulative Incidence of competing events
#'
#' @description Estimation of Cumulative Incidence Functions (CIF) of competing events.
#' @description This function is based on the CIF estimated by the survival package.
#' @param df A data frame containing, at a minimum, exit, event, and exposure.
#' @param time Name of the column in df containing times of event or censoring.
#' @param event Name of the column in df containing codes for censoring (0) and event types (1-4). Analysis of more than 4 competing events is not supported by this function.
#' @param exposed Name of the column in df containing a binary (0/1) exposure variable for stratification.
#' @param entry Name of the column in df containing late entry times.
#' @param weights Name of the column in df containing user-supplied weights. If ipwvars is utilized, this argument is ignored.
#' @param ipwvars A vector of names of columns in `df` containing predictor variables for building a propensity score model for exposure and creating standardized inverse probability weights using this model. Overrides the weights argument.
#' @param print.attr A logical indicator for whether results should be returned in console.
#' @importFrom "stats" "approx" "as.formula" "glm" "predict" "quantile" "sd" "stepfun"
#' @importFrom "survival" "survfit" "Surv"
#' @return A data frame with the following columns:
#' \describe{
#'   \item{event}{Type of event that occurs at the given time.}
#'   \item{exposure}{Exposure group in which the event happens.}
#'   \item{time}{Time of the event.}
#'   \item{CIoinc_comp}{Value of the unexposed (denoted by “o”) composite cumulative incidence at the given time.}
#'   \item{CIxinc_comp}{Value of the exposed (denoted by “x”) composite cumulative incidence at the given time.}
#'   \item{CIoinc_1}{Value of the unexposed cumulative incidence of event 1 at the given time.}
#'   \item{CIxinc_1}{Value of the exposed cumulative incidence of event 1 at the given time.}
#'   \item{R_1}{Sub-hazard ratio/Cause-specific hazard ratio for event 1.}
#'   \item{R_2}{Sub-hazard ratio/Cause-specific hazard ratio for event 2.}
#' }
#' @examples
#' #data from the package
#' data <- hrcomprisk::dat_ckid
#' #Estimate the Cumulative Incidence Functions and Ratios of sHR and csHR
#' mydat.CIF<-CRCumInc(df=data, time=exit, event=event, exposed=b1nb0, print.attr=TRUE)
#' @export

CRCumInc <- function(df,time,event,exposed,entry=NULL,weights=NULL,ipwvars=NULL,print.attr=T){
  time <- round(df[[deparse(substitute(time))]], 4)
  event <- df[[deparse(substitute(event))]]
  if(is.factor(event)) event <- as.numeric(event)-1
  exposed <- df[[deparse(substitute(exposed))]]
  if(is.factor(exposed)) exposed <- as.numeric(exposed)-1
  nevents <- max(event)
  composite <- 0 + (event > 0) #0 is added to convert logical values to numeric values
  entry <- df[[deparse(substitute(entry))]]
  if(is.null(entry)) entry <- rep(0,length(time))
  entry <- round(entry, 4)

  # User-provided weights
  weights <- df[[deparse(substitute(weights))]]
  if(is.null(weights)) weights <- rep(1,length(time))

  # Model-generated weights
  if(!is.null(ipwvars)){
    ipwdf <- as.data.frame(cbind(exposed,df[,ipwvars]))
    ipwmodel <- glm(as.formula(ipwdf), family="binomial", data=ipwdf)
    #fweights <- formula(paste(deparse(substitute(exposed)) , " ~ ", paste(ipwvars, collapse=" + ")))
    #ipwmodel <- glm(fweights, family="binomial", df)
    #print(summary(ipwmodel))
    ipwpred <- predict(ipwmodel,type="response")
    sipw <- as.numeric()
    sipw[exposed==1] <- 1/ipwpred[exposed==1] * (sum(exposed==1)/length(exposed))
    sipw[exposed==0] <- 1/(1-ipwpred[exposed==0]) * (sum(exposed==0)/length(exposed))
    tts <- sd(log(sipw))
    zipw <- log(sipw)/tts
    # do we want to allow the user to provide the value for extreme tts?
    sipw[zipw > 4] <- exp(4*tts)
    sipw[zipw < -4] <- exp(-4*tts)
    weights <- sipw
    #print(quantile(weights))
  }

  # Adjust the times in order to distinguish separate events when times are tied, with censoring after events
  time <- time + 1e-005 * event #this will fail if nevents>=5
  time[event == 0] <- time[event == 0] + 4.9e-005
  # so that times still round to original values
  # unexposed
  t0 <- time[exposed == 0]
  c0 <- composite[exposed == 0]
  w0 <- entry[exposed == 0]
  wt0 <- weights[exposed == 0]
  CIcomposite0 <- 1 - survfit(Surv(w0,t0,c0)~1, type = "kaplan-meier", weights=wt0)$surv
  CIcomposite0 <- c(0, CIcomposite0)
  ttd <- diff(CIcomposite0)
  # diff returns the change in the CI at each time whether uncensored or censored
  jumps0 <- ttd[ttd != 0]
  #values of the "jumps" in the CI curve
  utimes0 <- sort(unique(t0[c0 == 1]))
  # ordered distinct event times, where "jumps" occur
  #exposed
  t1 <- time[exposed == 1]
  c1 <- composite[exposed == 1]
  w1 <- entry[exposed == 1]
  wt1 <- weights[exposed == 1]
  CIcomposite1 <- 1 - survfit(Surv(w1,t1,c1)~1, type = "kaplan-meier", weights=wt1)$surv
  CIcomposite1 <- c(0, CIcomposite1)
  ttd <- diff(CIcomposite1)
  jumps1 <- ttd[ttd != 0]
  utimes1 <- sort(unique(t1[c1 == 1]))
  #recombine groups
  utimesall <- c(utimes0, utimes1)
  exposedall <- c(rep(0, length(utimes0)), rep(1, length(utimes1)))
  jumpsall <- c(jumps0, jumps1)
  #separate event types using fifth decimal place
  evtype <- rep(0, length(utimesall))
  evtype[round(utimesall - round(utimesall, 4), 5) == 1e-005] <- 1
  evtype[round(utimesall - round(utimesall, 4), 5) == 2e-005] <- 2
  evtype[round(utimesall - round(utimesall, 4), 5) == 3e-005] <- 3
  evtype[round(utimesall - round(utimesall, 4), 5) == 4e-005] <- 4
  #separate event- and exposure-specific cumulative incidences
  ret <- NULL
  uorder <- order(utimesall)
  ret$alltimes <- utimesall[uorder]
  ret$CI0inc <- NULL
  ret$CI1inc <- NULL
  for(i in 1:nevents){
    ret$CI0inc[[i]] <- c(0,approx(c(0, utimesall[evtype == i & exposedall == 0]),c(0, cumsum(jumpsall[evtype == i & exposedall == 0])),ret$alltimes,method="constant",rule=2)$y)
    ret$CI1inc[[i]] <- c(0,approx(c(0, utimesall[evtype == i & exposedall == 1]),c(0, cumsum(jumpsall[evtype == i & exposedall == 1])),ret$alltimes,method="constant",rule=2)$y)
  }
  ret$alltimes <- c(0,ret$alltimes)
  ret$evtype<-c(NA, evtype[uorder])
  ret$exposure<-c(NA, exposedall[uorder])
  #Creating a dataframe
  ##Composite dataframe
  ci0comp <- cbind(cumsum(jumps0), utimes0)
  ci1comp <- cbind(cumsum(jumps1), utimes1)
  utimesall1<-as.data.frame(utimesall)
  comp.dat<-merge(utimesall1, ci0comp, by.x="utimesall", by.y = "utimes0", all = TRUE, incomparables=0)
  comp.dat<-merge(comp.dat, ci1comp, by.x="utimesall", by.y = "utimes1", all = TRUE, incomparables=0)
  v <- !is.na(comp.dat$V1.x)
  comp.dat$V1.x<-c(0, comp.dat$V1.x[v])[cumsum(v)+1]
  v <- !is.na(comp.dat$V1.y)
  comp.dat$V1.y<-c(0, comp.dat$V1.y[v])[cumsum(v)+1]
  comp.dat <- rbind(0, comp.dat)
  #Final data
  myCIF<-as.data.frame(cbind(ret$evtype, ret$exposure, ret$alltimes, comp.dat$V1.x, comp.dat$V1.y))
  colnames(myCIF)<-c("event", "exposure", "time", "CIoinc_comp", "CIxinc_comp"  )
  for(i in 1:nevents){
    myCIF<-as.data.frame(cbind(myCIF, ret$CI0inc[[i]], ret$CI1inc[[i]]))
    l<-dim(myCIF)
    l<-l[2]
    j<-l-1
    names(myCIF)[l]<-paste0("CIxinc_", i) #name follows Cumulative Incidence among the EXPOSED for EVENT 'i'
    names(myCIF)[j]<-paste0("CIoinc_", i) #name follows Cumulative Incidence among the UNEXPOSED for EVENT 'i'
  }
  CIx3 <- myCIF$CIxinc_3
  if(is.null(CIx3)) CIx3 <- 0*myCIF$CIxinc_1
  CIo3 <- myCIF$CIoinc_3
  if(is.null(CIo3)) CIo3 <- 0*myCIF$CIxinc_1
  CIx4 <- myCIF$CIxinc_4
  if(is.null(CIx4)) CIx4 <- 0*myCIF$CIxinc_1
  CIo4 <- myCIF$CIoinc_4
  if(is.null(CIo4)) CIo4 <- 0*myCIF$CIxinc_1
  myCIF$R1 <- (1 - (myCIF$CIxinc_2+CIx3+CIx4)/(1 - myCIF$CIxinc_1))/(1 - (myCIF$CIoinc_2+CIo3+CIo4)/(1 - myCIF$CIoinc_1))
  myCIF$R1[!is.finite(myCIF$R1)] <- 99 #just in case
  myCIF$R2 <- (1 - (myCIF$CIxinc_1+CIx3+CIx4)/(1 - myCIF$CIxinc_2))/(1 - (myCIF$CIoinc_1+CIo3+CIo4)/(1 - myCIF$CIoinc_2))
  myCIF$R2[!is.finite(myCIF$R2)] <- 99
  if(nevents>=3){
    myCIF$R3 <- (1 - (myCIF$CIxinc_1+myCIF$CIxinc_2+CIx4)/(1 - CIx3))/(1 - (myCIF$CIoinc_1+myCIF$CIoinc_2+CIo4)/(1 - CIo3))
    myCIF$R3[!is.finite(myCIF$R3)] <- 99
  }
  if(nevents==4){
    myCIF$R4 <- (1 - (myCIF$CIxinc_1+myCIF$CIxinc_2+CIx3)/(1 - CIx4))/(1 - (myCIF$CIoinc_1+myCIF$CIoinc_2+CIo3)/(1 - CIo4))
    myCIF$R4[!is.finite(myCIF$R4)] <- 99
  }
  if(print.attr) print(attributes(myCIF)[-2])
  return(myCIF)
}
