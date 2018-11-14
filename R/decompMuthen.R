#' DEPRECIATED:
#'
#' Desc
#'
#' Details
#'
#' @param allIndicators
#' @param l1Var
#' @param l2Var
#' @param df
#' @param n
#' @param g
#'
#' @export
# -----------------------------------------------------------------------------

decompMuthen <- function(allIndicators,
                         l1Var,
                         l2Var,
                         df,
                         n,
                         g) {
  #-------------------------------------------------------------------
  # some general setup
  # ------------------------------------------------------------------
  grpMnNames <- paste(allIndicators, "_gmn", sep = "")
  mnNames    <- paste(allIndicators, "_mn" , sep = "")
  grpSizeNames <- paste(allIndicators, "_n" , sep = "")
  # creating an aggregated data set (with group means and group sizes)
  agg <- cbind(aggregate(df[allIndicators],
                         by = list(cluster = df[[l2Var]]),
                         FUN = mean),
               aggregate(df[l2Var],
                         by = list(cluster = df[[l2Var]]),
                         FUN = length)
  )
  # dropping the repeat "cluster" var
  agg <- agg[,-(ncol(agg)-1)]
  # adding names
  names(agg) <- c(l2Var, grpMnNames, "grpSize")

  #-------------------------------------------------------------------
  # within covariance matrices
  # ------------------------------------------------------------------
  # merge agg with full df
  df <- merge(df, agg, by = l2Var)
  # Create difference matrix
  d <- as.matrix(df[allIndicators]-df[grpMnNames])
  # Create covariance matrix
  wCovMat <- ((1/(n-g-1)) * t(d) %*% d)

  #-------------------------------------------------------------------
  # between covariance matrices
  # ------------------------------------------------------------------
  # average cluster size
  nMn <- mean(agg[['grpSize']])
  # save overall average for each indicator
  for (i in seq(allIndicators)) {
    agg <- cbind(agg, mean(df[[allIndicators[i]]]))
  }
  names(agg) <- c(l2Var, grpMnNames, "grpSize", mnNames)
  # difference matrix
  d <- as.matrix(agg[mnNames]-agg[grpMnNames])
  # covariance matrix
  # bCovMat <- ((nMn/(g-1)) * t(d) %*% d)
  bCovMat <- ((g-1) * t(nMn*d) %*% d)
  # names
  colnames(bCovMat)  <- allIndicators
  rownames(bCovMat)  <- allIndicators
  # return
  return(list(within=wCovMat, between=bCovMat))

}
