#' Decompose two level by Goldstein method
#'
#' Desc:
#'
#' Details to be written
#'
#' @param allIndicators character vector of indicator names
#' @param l1Var string
#' @param l2Var string
#' @param df dataframe
#'
#' @export
# -----------------------------------------------------------------------------

decompGoldstein <- function(allIndicators,
                            l1Var,
                            l2Var,
                            df) {
  # make long (necessary for fitting in nlme)
  long <- reshape2::melt(df, id.vars = c(l1Var, l2Var), variable.name = "item")
  # dummy code the indicators
  for (level in unique(long$item)) {
    long[level] <- ifelse(long$item == level, 1, 0)
  }
  # save combinations of indicators
  combinations <- t(combn(allIndicators, m=2))
  # create the matrices
  wCovMat <- matrix(nrow=length(allIndicators), ncol=length(allIndicators),
                    dimnames = list(allIndicators, allIndicators))
  bCovMat <- matrix(nrow=length(allIndicators), ncol=length(allIndicators),
                    dimnames = list(allIndicators, allIndicators))
  # loop through all bivariate combinations. for each fit and pull out covariances with getVCE()
  for (i in seq(nrow(combinations))) {
    # indicators in question
    i1 <- combinations[i,1]
    i2 <- combinations[i,2]
    # Create the model and random statements
    form <- paste(i1, i2 , sep = "+" )
    model <- as.formula(paste0("value~ -1 +", form))
    ranef <- as.formula(paste("~-1 + ", form, "|", l2Var, "/", l1Var))
    # subset the rows we want
    subset <- dplyr::filter(.data = long, item %in% c(i1, i2))
    # Try catch incase models do not converge
    tryCatch({
      #print(i)
      fit <- nlme::lme(fixed   = model,
                       random  = ranef,
                       data    = subset,
                       method  = "REML",
                       control = nlme::lmeControl(opt='optim'))
      covMats <- getVCE(fit, names = c(i1, i2))
      # fill in the larger
      wCovMat[i1, i2] <- covMats$within[i1, i2]
      wCovMat[i1, i1] <- covMats$within[i1, i1]
      wCovMat[i2, i2] <- covMats$within[i2, i2]
      wCovMat[i2, i1] <- covMats$within[i2, i1]
      # fill in the larger between
      bCovMat[i1, i2] <- covMats$between[i1, i2]
      bCovMat[i1, i1] <- covMats$between[i1, i1]
      bCovMat[i2, i2] <- covMats$between[i2, i2]
      bCovMat[i2, i1] <- covMats$between[i2, i1]
    },
    error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  } # end of for loop making covariance matrices
  return(list(within=wCovMat, between=bCovMat))
} # end of the function
