## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%",
  fig.width = 10,
  fig.height = 8
)

## ----example------------------------------------------------------------------
library(hrcomprisk)
data <- hrcomprisk::dat_ckid
dim(data) #dimensions
names(data) #varible names

## ----CRCumInc_example---------------------------------------------------------
mydat.CIF<-CRCumInc(df=data, time=exit, event=event, exposed=b1nb0, print.attr=T)

## ----plotCIF------------------------------------------------------------------
plots<-plotCIF(cifobj=mydat.CIF, maxtime = 20, eoi = 1)

## ----boot_example-------------------------------------------------------------
ciCIF<-bootCRCumInc(df=data, exit=exit, event=event, exposure=b1nb0, rep=100, print.attr=T)

## ----plot_ci------------------------------------------------------------------
plotCIF(cifobj=mydat.CIF, maxtime= 20, ci=ciCIF)

## ----npcrest------------------------------------------------------------------
npcrest(df=data, exit=exit, event=event, exposure=b1nb0,rep=100, maxtime=20, print.attr=T)

