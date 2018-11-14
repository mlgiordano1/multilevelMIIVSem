#' MIIV-2SLS estimation for multilevel cfa
#'
#' Description
#'
#' Details:
#'
#' @param withinModel Model statement using Lavaan/MIIVsem syntax
#' @param betweenModel Model statement using Lavaan/MIIVsem syntax
#' @param estimator "muthen" or "goldstein"
#' @param allIndicators character vector of the names of indicators
#' @l1Var string for level 1 variable name
#' @l2Var string for level 2 variable name
#' @df dataframe
#' @var.cov boolean. if you want variance covariance parameters estimated
#'
#' @export
# -----------------------------------------------------------------------------

# return within and between models
mlcfaMIIV <- function(withinModel,
                      betweenModel,
                      estimator = "muthen",
                      allIndicators,
                      l1Var,
                      l2Var,
                      df,
                      var.cov = FALSE) {
  # Program some checks, like is long a DF, fitWith =nlme or lmer,
  # all indicators is charater, etc.

  # -----------------------------------------------------------------
  # Do some universal steps
  # number of subjects
  n <- nrow(df)
  # number of groups
  g <- length(unique(df[[l2Var]]))
  # -----------------------------------------------------------------


  # -----------------------------------------------------------------
  # Process the data
  if (tolower(estimator) == "muthen") {
    # decompose muthen style
    muth <- mcfa.input(l1Var, l2Var, dat = df)
    covMats <- list(within = muth$pw.cov, between=muth$ab.cov)
  } else if (tolower(estimator) == "goldstein") {
    # decompose muthen style
    covMats <- decompGoldstein(allIndicators, l1Var, l2Var, df)
  } else {
    print("incorrect estimator entered")
  }

  # -----------------------------------------------------------------
  # fit with MIIVsem
  # Fit covariance matrices with MIIVsem
  w <- MIIVsem::miive(withinModel,
                      sample.cov = covMats[["within"]],
                      sample.nobs = n,
                      var.cov = var.cov)
  b <- MIIVsem::miive(betweenModel,
                      sample.cov = covMats[["between"]],
                      sample.nobs = g,
                      # overid.degree = 1,
                      # overid.method = "random",
                      var.cov = var.cov)
  # return the list of within and between models
  return(list(within=w, between=b))
} # end function
