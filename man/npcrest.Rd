% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/npcrest.R
\name{npcrest}
\alias{npcrest}
\title{Nonparametric Estimation of the Ratios of Hazard Ratios}
\usage{
npcrest(
  df,
  exit,
  event,
  exposure,
  entry = NULL,
  weights = NULL,
  ipwvars = NULL,
  maxtime = Inf,
  rep = NULL,
  eoi = -1,
  print.attr = T
)
}
\arguments{
\item{df}{A data frame containing, at a minimum, exit, event, and exposure.}

\item{exit}{Name of the column in df containing times of event or censoring.}

\item{event}{Name of the column in df containing codes for censoring (0) and event types (1-4). Analysis of more than 4 competing events is not supported by this function.}

\item{exposure}{Name of the column in df containing a binary (0/1) exposure variable for stratification.}

\item{entry}{Name of the column in df containing late entry times.}

\item{weights}{Name of the column in df containing user-supplied weights. If ipwvars is utilized, this argument is ignored.}

\item{ipwvars}{A vector of names of columns in `df` containing predictor variables for building a propensity score model for exposure and creating standardized inverse probability weights using this model. Overrides the weights argument.}

\item{maxtime}{Largest time to display on the x-axis of all output plots. As data can become sparse and thus more widely variable at times get large, this argument may be used to restrict plots to a range of the data that is discerned to be more accurate and reliable.}

\item{rep}{Number of replicates for bootstrapping if confidence intervals for the sHR/csHR estimate are desired. See more details on bootstrapping below.}

\item{eoi}{Event number for the event of interest, useful when more than two events exist. If utilized, only two cumulative incidence curves will be plotted: one for the event of interest, and one for the composite of all competing events. Each event will still have its sHR/csHR ratio plotted.}

\item{print.attr}{A logical indicator for whether results should be returned in console.}
}
\value{
An object containing the plotted figures (\code{$plots}) and a data frame (\code{$cuminc}) with the following columns:
\describe{
  \item{event}{Type of event that occurs at the given time.}
  \item{exposure}{Exposure group in which the event happens.}
  \item{time}{Time of the event.}
  \item{CIoinc_comp}{Value of the unexposed (denoted by “o”) composite cumulative incidence at the given time.}
  \item{CIxinc_comp}{Value of the exposed (denoted by “x”) composite cumulative incidence at the given time.}
  \item{CIoinc_1}{Value of the unexposed cumulative incidence of event 1 at the given time.}
  \item{CIxinc_1}{Value of the exposed cumulative incidence of event 1 at the given time.}
  \item{R_1}{Sub-hazard ratio/Cause-specific hazard ratio for event 1.}
  \item{R_2}{Sub-hazard ratio/Cause-specific hazard ratio for event 2.}
}
}
\description{
A comprehensive wrapper function for implementing the competing risks diagnostic of Ng, Antiporta, Matheson and Muñoz (2019) to compare sub-hazard ratio (a la Fine and Gray; sHR) and cause-specific hazard ratio (csHR) approaches. While doing either analysis individually involves parametric or semi-parametric estimation, because of the fact that the derivatives of the cumulative incidences involved in the quantities cancel out when their ratio is considered, this ratio can be characterized completely using nonparametric estimates of the event-specific cumulative incidences. This provides a useful diagnostic when both analyses are performed as well as a method for estimating the sub-hazard ratios in a way that is valid and free of tethering assumptions characterized by Muñoz et al. (2013).

This function calls datcheck, CRCumInc, bootCRCumInc (if confidence intervals are requested), and plotCIF, all of which are also included in the hrcomprisk package. These functions should generally not be utilized directly.
}
\details{
\strong{1.	Bootstrapped confidence intervals for the sHR/csHR quantities}

If a positive number of bootstrap replicates is requested via the rep argument, the program will calculate and provide pointwise percentile-based bootstrap confidence intervals for the sHR/csHR ratios. The bootstrapping process uses two loops. In the first loop, rep bootstrap samples are taken stratified by exposure (so each sample has the same exposure prevalence as the original data) and all event-specific cumulative incidences are calculated and stored for each of them. In the second loop, for each event time (i.e., each change in any one of the cumulative incidence functions), the 2.5th and 97.5th percentiles of the bootstrap estimates of the rep sHR/csHR ratios are stored as the lower and upper confidence limits. These are not directly returned to the user, but are used in the plotting of the sHR/csHR ratios.

\strong{2.	User-supplied vs. program-generated weights}

If confidence intervals are not desired, the user can supply a column of weights (e.g, inverse probability weights from a model predicting exposure) which will be used in the estimation of the cumulative incidences and the sHR/csHR ratios derived from them. However, the use of bootstrapping for calculation of confidence intervals as described above necessitates that such a model be refit for each bootstrap sample, generating new weights for new estimates of all quantities. If this is desired, the user should include all the predictor variables as columns in the data frame so that the appropriate model can be fit automatically using the ipwvars argument. If this method is used, the program uses a logistic regression model to calculate probability of exposure, stabilizes the resulting weights to the sample size, and winsorizes weights that fall outside ±4 standard deviations on the log scale. Using this method can increase computation time as the model must be refit on each bootstrap replicate.

\strong{3.	Use of the nonparametric sHR/csHR ratio for calculation of sHR estimates}

Simultaneous estimation of all subhazard ratios and cause-specific hazard ratios is often fraught with problematic results due to incompatible modeling assumptions (e.g., not all ratios can be proportional) and the tethering inherent in subhazard analysis of multiple events as described by Muñoz et al. (2013). Cause-specific hazard ratios are not subject to such tethering constraints, and will admissible regardless of whether the model is misspecified. As such, the output of this function – valid nonparametric estimates of the sHR/csHR ratio – can be combined with (i.e., multiplied by) cause-specific hazard ratio estimates (e.g., from a proportional or loglinear cause-specific hazards model) to produce subhazard ratio estimates which do not violate the principles of tethering. See below for an example of how to implement this.
}
\examples{
#data from the package
data <- hrcomprisk::dat_ckid
#Using the wrapper function
npcrest(df=data, exit=exit, event=event, exposure=b1nb0,rep=10, maxtime=20, print.attr=TRUE)
}
\references{
1. Ng D, Antiporta DA, Matheson M, Munoz A. Nonparametric assessment of differences between competing risks hazard ratios: application to racial differences in pediatric chronic kidney disease progression. (Clinical Epidemiology, 2019-in print)

2. Muñoz A, Abraham AG, Matheson M, Wada N. In: Risk Assessment and Evaluation of Predictions. Lee MLT, Gail M, Pfeiffer R, Satten G, Cai T, Gandy A, editor. New York: Springer; 2013. Non-proportionality of hazards in the competing risks framework; pp. 3–22. [Google Scholar](https://link.springer.com/chapter/10.1007/978-1-4614-8981-8_1)
}
