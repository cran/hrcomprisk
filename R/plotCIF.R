#' Plot Incidence and Ratio of sHR/csHR
#'
#' @description Plot Cumulative Incidence and Ratio of sHR/csHR.
#' @param cifobj A dataframe containing the Cumulative Incidence of each competing event by exposure group.
#' @param maxtime Largest time to display on the x-axis of all output plots.
#' @param ci A dataframe containing the 95\% CI for each ratio of csHR/sHR.
#' @param eoi Event number for the event of interest, useful when more than two events exist.
#' @return A large list containing 2 figures:
#' \describe{
#'   \item{$plot1 }{  Plot the cumulative incidence of the composite event and of each event by exposure group.}
#'   \item{$plot2 }{  Plot the ratio of Sub-hazard ratio and Cause-specific hazard ratio for each event i (Ri).}
#' }
#' @examples
#' #data from the package
#' data <- hrcomprisk::dat_ckid
#' #Estimate the Cumulative Incidence Functions and Ratios of sHR and csHR
#' mydat.CIF<-CRCumInc(df=data, time=exit, event=event, exposed=b1nb0, print.attr=FALSE)
#' #Plot the CIs and Ratios estimated
#' plots<-plotCIF(cifobj=mydat.CIF, maxtime = 20, eoi = 1)
#' @importFrom "grDevices" "gray" "recordPlot"
#' @importFrom "graphics" "abline" "axis" "box" "lines" "mtext" "par" "plot" "text"
#' @export

plotCIF <- function (cifobj, maxtime = Inf, ci = NULL, eoi=-1)
{
  time <- cifobj$time[cifobj$time <= maxtime]
  I1o <- cifobj$CIoinc_1[cifobj$time <= maxtime]
  I1x <- cifobj$CIxinc_1[cifobj$time <= maxtime]
  I2o <- cifobj$CIoinc_2[cifobj$time <= maxtime]
  I2x <- cifobj$CIxinc_2[cifobj$time <= maxtime]
  if(ncol(cifobj)>12){
    I3o <- cifobj$CIoinc_3[cifobj$time <= maxtime]
    I3x <- cifobj$CIxinc_3[cifobj$time <= maxtime]
  }
  else{
    I3o <- 0*I1o
    I3x <- 0*I1x
  }
  if(ncol(cifobj)>15){
    I4o <- cifobj$CIoinc_4[cifobj$time <= maxtime]
    I4x <- cifobj$CIxinc_4[cifobj$time <= maxtime]
  }
  else{
    I4o <- 0*I1o
    I4x <- 0*I1x
  }
  Io <- cifobj$CIoinc_comp[cifobj$time <= maxtime]
  Ix <- cifobj$CIxinc_comp[cifobj$time <= maxtime]
  if(eoi<=1) {J1o <- I1o; J1x <- I1x; J2o <- I2o+I3o+I4o; J2x <- I2x+I3x+I4x}
  if(eoi==2) {J1o <- I2o; J1x <- I2x; J2o <- I1o+I3o+I4o; J2x <- I1x+I3x+I4x}
  if(eoi==3) {J1o <- I3o; J1x <- I3x; J2o <- I2o+I1o+I4o; J2x <- I2x+I1x+I4x}
  if(eoi==4) {J1o <- I4o; J1x <- I4x; J2o <- I2o+I3o+I1o; J2x <- I2x+I3x+I1x}
  cifplots <- NULL

  #exit call
  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))

  par(mgp = c(2.25, 0.5, 0), las = 1, mfrow = c(1, 2), mar = c(5.1, 4.1, 2.1, 2.1), cex.lab = 1, cex.axis = 1)
  plot(0, 0, pch = 16, xlim = c(0, max(time)), ylim = c(0, 1), main = "", ylab = "", xlab = "Time", cex.lab = 1, font.main = 12, type = "n", yaxt = "n", xaxs = "i")
  mtext("Cumulative incidence of composite event", at = 0.5, side = 2, line = 1.5, las = 0, col = gray(0.25))
  axis(2, at = seq(0, 1, 0.2), col.axis = gray(0.25), labels = seq(0, 100, 20))
  lines(stepfun(time, c(0, Ix), f = 0), lty = 1, cex = 0, lwd = 3, col = gray(0.25))
  lines(stepfun(time, c(0, Io), f = 0), cex = 0, lwd = 1, col = gray(0.25))
  par(mar = c(5.1, 4.1, 2.1, 4.1))

  if(ncol(cifobj)==11 | eoi>0){ # funnel plot
    plot(0, 0, pch = 16, xlim = c(0, max(time)), ylim = c(0, 1), main = "", ylab = "", xlab = "Time", cex.lab = 1, font.main = 12, type = "n", yaxt = "n", xaxs = "i")
    ymax <- round(5*max(c(J1o,J1x)))/5
    ymax2 <- round(5*max(c(J2o,J2x)))/5
    mtext("Cumulative incidence of event of interest", side = 2, at = ymax/2, las = 0, line = 1.5, col = gray(0.5))
    corners <- par("usr")
    if(ncol(cifobj)==11) {text(corners[2] + (3/22)*diff(corners[1:2]), 1-(ymax/2), "Cumulative incidence of competing event", srt = -90, xpd = T)}
    else {text(corners[2] + (3/22)*diff(corners[1:2]), 1-(ymax2/2), "Cumulative incidence of all competing events", srt = -90, xpd = T)}
    axis(2, at = seq(0, ymax, 0.2), col.axis = gray(0.5), labels = 100*seq(0, ymax, 0.2))
    axis(4, at = seq(1, 1-ymax2, -0.2), labels = 100*seq(0, ymax2, 0.2))
    box()
    lines(stepfun(time, c(0, J1x)), lty = 1, cex = 0, lwd = 3, col = gray(0.5))
    lines(stepfun(time, c(1, 1-J2x)), cex = 0, lwd = 3, col = gray(0))
    lines(stepfun(time, c(0, J1o)), cex = 0, lwd = 1, col = gray(0.5))
    lines(stepfun(time, c(1, 1-J2o)), cex = 0, lwd = 1, col = gray(0))
  } #endif2

  if(ncol(cifobj)==14 & eoi<0){ # stack plot
    plot(0, 0, pch = 16, xlim = c(0, max(time)), ylim = c(0, 1), main = "", ylab = "", xlab = "Time", cex.lab = 1, font.main = 12, type = "n", yaxt = "n", xaxs = "i")
    mtext("Cumulative incidence of events", side = 2, at = 0.5, las = 0, line = 1.5, col = gray(0))
    axis(2, at = seq(0, 1, 0.2), col.axis = gray(0), labels = 100*seq(0, 1, 0.2))
    box()
    lines(stepfun(time, c(0, I1o)), cex = 0, lwd = 1, col = gray(0.6))
    lines(stepfun(time, max(I1o)+c(0, I2o)), cex = 0, lwd = 1, col = gray(0.3))
    lines(stepfun(time, max(I1o+I2o)+c(0, I3o)), cex = 0, lwd = 1, col = gray(0))
    lines(stepfun(time, c(0, I1x)), cex = 0, lwd = 3, col = gray(0.6))
    lines(stepfun(time, max(I1x)+c(0, I2x)), cex = 0, lwd = 3, col = gray(0.3))
    lines(stepfun(time, max(I1x+I2x)+c(0, I3x)), cex = 0, lwd = 3, col = gray(0))
  } #endif3

  if(ncol(cifobj)==17 & eoi<0){ # stack plot
    plot(0, 0, pch = 16, xlim = c(0, max(time)), ylim = c(0, 1), main = "", ylab = "", xlab = "Time", cex.lab = 1, font.main = 12, type = "n", yaxt = "n", xaxs = "i")
    mtext("Cumulative incidence of events", side = 2, at = 0.5, las = 0, line = 1.5, col = gray(0))
    axis(2, at = seq(0, 1, 0.2), col.axis = gray(0), labels = 100*seq(0, 1, 0.2))
    box()
    lines(stepfun(time, c(0, I1o)), cex = 0, lwd = 1, col = gray(0.6))
    lines(stepfun(time, max(I1o)+c(0, I2o)), cex = 0, lwd = 1, col = gray(0.4))
    lines(stepfun(time, max(I1o+I2o)+c(0, I3o)), cex = 0, lwd = 1, col = gray(0.2))
    lines(stepfun(time, max(I1o+I2o+I3o)+c(0, I4o)), cex = 0, lwd = 1, col = gray(0))
    lines(stepfun(time, c(0, I1x)), cex = 0, lwd = 3, col = gray(0.6))
    lines(stepfun(time, max(I1x)+c(0, I2x)), cex = 0, lwd = 3, col = gray(0.4))
    lines(stepfun(time, max(I1x+I2x)+c(0, I3x)), cex = 0, lwd = 3, col = gray(0.2))
    lines(stepfun(time, max(I1x+I2x+I3x)+c(0, I4x)), cex = 0, lwd = 3, col = gray(0))
  } #endif4
  cifplots$plot1 <- recordPlot()

  #ttx <- sort(data$time)
  #ttx <- ttx[ttx < maxtime]
  ntimes <- length(time)
  R1 <- cifobj$R1[cifobj$time <= maxtime]
  R2 <- cifobj$R2[cifobj$time <= maxtime]
  R3 <- cifobj$R3[cifobj$time <= maxtime]
  if(is.null(R3)) R3 <- rep(1,ntimes)
  R4 <- cifobj$R4[cifobj$time <= maxtime]
  if(is.null(R4)) R4 <- rep(1,ntimes)
  yrange <- range(list(R1,R2,R3,R4,ci))
  par(mar = c(5.1, 4.1, 2.1, 2.1), mfrow = c(1, 2), pch = -1)
  plot(1, 1, type = "n", xlim = c(0, max(time)), ylim = yrange, xlab = "Time", ylab = expression(paste("Sub-hazard ratio / Cause-specific hazard ratio ( ", R[1], ")")), log = "y", xaxs = "i")
  xmid <- mean(par("usr")[1:2])
  abline(h = 1, lty = 2, col = gray(0.5))
  box()
  lines(stepfun(time, c(1, R1)), col = gray(0), lwd = 3, cex = 0)
  #text(xmid, quantile(yrange,0.95), "Event 1", cex = 1.2, col = gray(0))
  mtext("Event 1",at=xmid,side=3,line=1,cex=1.2)
  if (!is.null(ci)) {
    lines(stepfun(time, c(1, ci$R1.lower[1:ntimes])), col = gray(0), lwd = 2, cex = 0)
    lines(stepfun(time, c(1, ci$R1.upper[1:ntimes])), col = gray(0), lwd = 2, cex = 0)
  }
  plot(1, 1, type = "n", xlim = c(0, max(time)), ylim = yrange, xlab = "Time", ylab = expression(paste("Sub-hazard ratio / Cause-specific hazard ratio ( ", R[2], ")")), log = "y", xaxs = "i")
  abline(h = 1, lty = 2)
  box()
  lines(stepfun(time, c(1, R2)), col = gray(0), lwd = 3, cex = 0)
  if (!is.null(ci)) {
    lines(stepfun(time, c(1, ci$R2.lower[1:ntimes])), col = gray(0), lwd = 2, cex = 0)
    lines(stepfun(time, c(1, ci$R2.upper[1:ntimes])), col = gray(0), lwd = 2, cex = 0)
  }
  #text(xmid, quantile(yrange,0.95), "Event 2", cex = 1.2, col = gray(0))
  mtext("Event 2",at=xmid,side=3,line=1,cex=1.2)
  cifplots$plot2 <- recordPlot()

  if(ncol(cifobj)>12){
    par(mar = c(5.1, 4.1, 2.1, 2.1), mfrow = c(1, 2), pch = -1)
    plot(1, 1, type = "n", xlim = c(0, max(time)), ylim = yrange, xlab = "Time", ylab = expression(paste("Sub-hazard ratio / Cause-specific hazard ratio ( ", R[3], ")")), log = "y", xaxs = "i")
    abline(h = 1, lty = 2)
    box()
    lines(stepfun(time, c(1, R3)), col = gray(0), lwd = 3, cex = 0)
    if (!is.null(ci)) {
      lines(stepfun(time, c(1, ci$R3.lower[1:ntimes])), col = gray(0), lwd = 2, cex = 0)
      lines(stepfun(time, c(1, ci$R3.upper[1:ntimes])), col = gray(0), lwd = 2, cex = 0)
    }
    #text(xmid, quantile(yrange,0.95), "Event 3", cex = 1.2, col = gray(0))
    mtext("Event 3",at=xmid,side=3,line=1,cex=1.2)
    if(ncol(cifobj)>15){
      plot(1, 1, type = "n", xlim = c(0, max(time)), ylim = yrange, xlab = "Time", ylab = expression(paste("Sub-hazard ratio / Cause-specific hazard ratio ( ", R[4], ")")), log = "y", xaxs = "i")
      abline(h = 1, lty = 2)
      box()
      lines(stepfun(time, c(1, R4)), col = gray(0), lwd = 3, cex = 0)
      if (!is.null(ci)) {
        lines(stepfun(time, c(1, ci$R4.lower[1:ntimes])), col = gray(0), lwd = 2, cex = 0)
        lines(stepfun(time, c(1, ci$R4.upper[1:ntimes])), col = gray(0), lwd = 2, cex = 0)
      }
      #text(xmid, quantile(yrange,0.95), "Event 4", cex = 1.2, col = gray(0))
      mtext("Event 4",at=xmid,side=3,line=1,cex=1.2)
    } #endif4
    cifplots$plot3 <- recordPlot()
  } #endif3
  invisible(cifplots)
}
