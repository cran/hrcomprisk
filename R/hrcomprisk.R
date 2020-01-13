#' Nonparametric Assessment Between Competing Risks Hazards
#'
#' Estimate nonparametric cumulative-incidence based estimation
#' of the ratios of sub-hazard ratios to cause-specific hazard ratios from Ng, Antiporta,
#' Matheson and Muñoz (2019)[1] to compare sub-hazard ratio (a la Fine and Gray; sHR) and cause-specific hazard ratio
#' (csHR) approaches.
#'
#' While doing either analysis individually involves parametric or semi-parametric estimation,
#' because of the fact that the derivatives of the cumulative incidences involved in the quantities cancel out when
#' their ratio is considered, this ratio can be characterized completely using nonparametric estimates of the event-specific
#' cumulative incidences. This provides a useful diagnostic when both analyses are performed as well as a method for estimating
#' the sub-hazard ratios in a way that is valid and free of tethering assumptions characterized by Muñoz et al. [2].
#'
#' @details \strong{1.	Bootstrapped confidence intervals for the sHR/csHR quantities}
#' @details If a positive number of bootstrap replicates is requested via the rep argument, the program will
#' calculate and provide pointwise percentile-based bootstrap confidence intervals for the sHR/csHR ratios.
#' The bootstrapping process uses two loops. In the first loop, rep bootstrap samples are taken stratified by
#' exposure (so each sample has the same exposure prevalence as the original data) and all event-specific cumulative
#' incidences are calculated and stored for each of them. In the second loop, for each event time (i.e., each change
#' in any one of the cumulative incidence functions), the 2.5th and 97.5th percentiles of the bootstrap estimates of
#' the rep sHR/csHR ratios are stored as the lower and upper confidence limits. These are not directly returned to the
#' user, but are used in the plotting of the sHR/csHR ratios.
#'
#' @details \strong{2.	User-supplied vs. program-generated weights}
#' @details If confidence intervals are not desired, the user can supply a column of weights (e.g, inverse probability weights from a model predicting exposure) which will be used in the estimation of the cumulative incidences and the sHR/csHR ratios derived from them. However, the use of bootstrapping for calculation of confidence intervals as described above necessitates that such a model be refit for each bootstrap sample, generating new weights for new estimates of all quantities. If this is desired, the user should include all the predictor variables as columns in the data frame so that the appropriate model can be fit automatically using the ipwvars argument. If this method is used, the program uses a logistic regression model to calculate probability of exposure, stabilizes the resulting weights to the sample size, and winsorizes weights that fall outside ±4 standard deviations on the log scale. Using this method can increase computation time as the model must be refit on each bootstrap replicate.
#'
#' @details \strong{3.	Use of the nonparametric sHR/csHR ratio for calculation of sHR estimates}
#' @details Simultaneous estimation of all subhazard ratios and cause-specific hazard ratios is often fraught with problematic results due to incompatible modeling assumptions (e.g., not all ratios can be proportional) and the tethering inherent in subhazard analysis of multiple events as described by Muñoz et al. [2]. Cause-specific hazard ratios are not subject to such tethering constraints, and will admissible regardless of whether the model is misspecified. As such, the output of this function – valid nonparametric estimates of the sHR/csHR ratio – can be combined with (i.e., multiplied by) cause-specific hazard ratio estimates (e.g., from a proportional or loglinear cause-specific hazards model) to produce subhazard ratio estimates which do not violate the principles of tethering.
#'
#' @section Functions:
#'
#' The hrcomprisk package provides 3 main functions and a wrapper function:
#'
#' - \link{npcrest} : Main wrapper function.
#'
#' - \link{CRCumInc} : Estimation of Cumulative Incidence Functions (CIF) of competing events.
#'
#' - \link{plotCIF} : Plot Cumulative Incidence and Ratio of sHR/csHR.
#'
#' - \link{bootCRCumInc} : Bootstrap 95% Confidence Intervals limits for estimated Ratios of sHR/csHR.
#'
#' @references
#' 1. Ng D, Antiporta DA, Matheson M, Munoz A. Nonparametric assessment of differences
#' between competing risks hazard ratios: application to racial differences in pediatric
#' chronic kidney disease progression. Clinical Epidemiology, 2020 (in press)
#'
#' 2. Muñoz A, Abraham AG, Matheson M, Wada N. In: Risk Assessment and Evaluation of Predictions.
#' Lee MLT, Gail M, Pfeiffer R, Satten G, Cai T, Gandy A, editor. New York: Springer; 2013.
#' Non-proportionality of hazards in the competing risks framework; pp. 3–22. \href{https://link.springer.com/chapter/10.1007/978-1-4614-8981-8_1}{[Google Scholar]}
#'
#' @author
#'
#' \strong{Mantainer:}  Daniel Antiporta <\email{dantiporta@@jhu.edu}>
#'
#' \strong{Authors:}
#'
#' - Daniel Antiporta
#'
#' - Matthew Matheson
#'
#' - Derek Ng
#'
#' - Alvaro Muñoz
#'
#' @seealso
#'
#' Useful links:
#'
#' \url{https://github.com/AntiportaD/hrcomprisk}
#'
#' @examples
#' #data from the package - See fuctions for specific examples
#' data <- hrcomprisk::dat_ckid
#' @docType package
#' @name hrcomprisk
NULL
