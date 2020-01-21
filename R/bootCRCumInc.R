#' Bootstrap for Ratios of Hazard Ratios
#'
#' @description Bootstrap 95\% Confidence Intervals limits for estimated Ratios of sHR/csHR.
#' @param df A data frame containing, at a minimum, exit, event, and exposure.
#' @param exit Name of the column in df containing times of event or censoring.
#' @param event Name of the column in df containing codes for censoring (0) and event types (1-4). Analysis of more than 4 competing events is not supported by this function.
#' @param exposure Name of the column in df containing a binary (0/1) exposure variable for stratification.
#' @param entry Name of the column in df containing late entry times.
#' @param weights Name of the column in df containing user-supplied weights. If ipwvars is utilized, this argument is ignored.
#' @param ipwvars A vector of names of columns in `df` containing predictor variables for building a propensity score model for exposure and creating standardized inverse probability weights using this model. Overrides the weights argument.
#' @param rep Number of replicates for bootstrapping if confidence intervals for the sHR/csHR estimate are desired. See more details on bootstrapping below.
#' @param print.attr A logical indicator for whether results should be returned in console.
#' @param seed A seed number start for the bootstrap estimation.
#' @importFrom "stats" "approx" "as.formula" "glm" "predict" "quantile" "sd" "stepfun"
#' @return A data frame with the 95\% confidence interval limits (upper and lower) for
#' Sub-hazard ratio/Cause-specific hazard ratio for each event:
#' \describe{
#'   \item{R1.lower}{Lower limit of the 95\%CI of the Sub-hazard ratio/Cause-specific hazard ratio for event 1 at time \code{t}}
#'   \item{R1.upper}{Upper limit of the 95\%CI of the Sub-hazard ratio/Cause-specific hazard ratio for event 1 at time \code{t}}
#'   \item{R2.lower}{Lower limit of the 95\%CI of the Sub-hazard ratio/Cause-specific hazard ratio for event 2 at time \code{t}}
#'   \item{R2.upper}{Upper limit of the 95\%CI of the Sub-hazard ratio/Cause-specific hazard ratio for event 2 at time \code{t}}
#' }
#' @examples
#' #data from the package
#' data <- hrcomprisk::dat_ckid
#' #Obtain the 95%CI by bootstraping
#' ciCIF<-bootCRCumInc(df=data, exit=exit, event=event, exposure=b1nb0, rep=10, print.attr=TRUE)
#' @export

bootCRCumInc <- function (df, exit, event, exposure, entry = NULL, weights = NULL, ipwvars=NULL, rep = 0, print.attr=T, seed=54321)
{
  df$exit <- df[[deparse(substitute(exit))]]
  df$event <- df[[deparse(substitute(event))]]
  if(is.factor(df$event)) df$event <- as.numeric(df$event)-1
  df$exposure <- df[[deparse(substitute(exposure))]]
  if(is.factor(df$exposure)) df$exposure <- as.numeric(df$exposure)-1
  df$entry <- df[[deparse(substitute(entry))]]
  df$weights <- df[[deparse(substitute(weights))]]
  set.seed(seed)
  e3 <- any(df$event==3)
  e4 <- any(df$event==4)
  if (rep == 0) {
    stop("\n", "`n` boostrapping repetitions missing", call. = FALSE)
  }
  else {
    nboot <- rep
    bI1o.all <- NULL
    bI1x.all <- NULL
    bI2o.all <- NULL
    bI2x.all <- NULL
    R1.all <- NULL
    R2.all <- NULL
    bI3o.all <- NULL
    bI3x.all <- NULL
    R3.all <- NULL
    bI4o.all <- NULL
    bI4x.all <- NULL
    R4.all <- NULL
    ttx.all <- NULL
    datx <- df[which(df$exposure == 1), ]
    dato <- df[which(df$exposure == 0), ]
    for (b in 1:nboot) {
      b.data <- rbind(datx[sample(1:nrow(datx),replace=T),],dato[sample(1:nrow(dato),replace=T),])
      b.data$exit <- jitter(b.data$exit)

      dat_boot <- do.call(CRCumInc,list(b.data,substitute(exit), substitute(event), substitute(exposure), substitute(entry), substitute(weights), substitute(ipwvars), F))
      ttx.all[[b]] <- dat_boot$time
      R1.all[[b]] <- dat_boot$R1
      R2.all[[b]] <- dat_boot$R2
      if(e3) R3.all[[b]] <- dat_boot$R3
      if(e4) R4.all[[b]] <- dat_boot$R4
      rm(dat_boot)
    }
    ttx <- sort(c(0, df$exit[df$event > 0]))
    R1.lower <- NULL
    R1.upper <- NULL
    R2.lower <- NULL
    R2.upper <- NULL
    if(e3){
      R3.lower <- NULL
      R3.upper <- NULL
    }
    if(e4){
      R4.lower <- NULL
      R4.upper <- NULL
    }
    for (i in 1:length(ttx)) {
      tty1 <- NULL
      tty2 <- NULL
      tty3 <- NULL
      tty4 <- NULL
      for (b in 1:nboot) {
        tty1[b] <- approx(ttx.all[[b]], R1.all[[b]], xout = ttx[i], method = "constant", f = 0)$y
        tty2[b] <- approx(ttx.all[[b]], R2.all[[b]], xout = ttx[i], method = "constant", f = 0)$y
        if(e3) tty3[b] <- approx(ttx.all[[b]], R3.all[[b]], xout = ttx[i], method = "constant", f = 0)$y
        if(e4) tty4[b] <- approx(ttx.all[[b]], R4.all[[b]], xout = ttx[i], method = "constant", f = 0)$y
      }
      ttq <- quantile(tty1, probs = c(0.025, 0.975), na.rm = T)
      R1.lower[i] <- ttq[1]
      R1.upper[i] <- ttq[2]
      ttq <- quantile(tty2, probs = c(0.025, 0.975), na.rm = T)
      R2.lower[i] <- ttq[1]
      R2.upper[i] <- ttq[2]
      if(e3){
        ttq <- quantile(tty3, probs = c(0.025, 0.975), na.rm = T)
        R3.lower[i] <- ttq[1]
        R3.upper[i] <- ttq[2]
      }
      if(e4){
        ttq <- quantile(tty4, probs = c(0.025, 0.975), na.rm = T)
        R4.lower[i] <- ttq[1]
        R4.upper[i] <- ttq[2]
      }
    }
    CRCumInc_boot <- as.data.frame(cbind(R1.lower, R1.upper, R2.lower, R2.upper))
    if(e3) CRCumInc_boot <- as.data.frame(cbind(CRCumInc_boot,R3.lower,R3.upper))
    if(e4) CRCumInc_boot <- as.data.frame(cbind(CRCumInc_boot,R4.lower,R4.upper))
    if(print.attr) print(attributes(CRCumInc_boot)[1:2])
    invisible(CRCumInc_boot)
  }
}
