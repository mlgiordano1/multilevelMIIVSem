#' Extract the variance parameter from nlme object
#'
#' Desc:
#'
#' Details to be written
#'
#' @param twoLvlnlmeObj nlme model obj
#' @param names
#'
#' @export
# -----------------------------------------------------------------------------

getVCE <- function(twoLvlnlmeObj, names) {
  # --------------------------------
  # do some checks?
  # --------------------------------
  # Save varcorr output
  v               <- nlme::VarCorr(twoLvlnlmeObj)
  # pull the between variances and save them on the diag
  b.cov           <- diag(v[2:3,1])
  # compute the covariance
  b.cov[2,1]      <- sqrt(b.cov[1,1])*sqrt(b.cov[2,2])*as.numeric(v[3,3])
  b.cov[1,2]      <- sqrt(b.cov[1,1])*sqrt(b.cov[2,2])*as.numeric(v[3,3])
  colnames(b.cov) <- names
  rownames(b.cov) <- names
  # now for the within
  w.cov           <- diag(v[5:6,1])
  w.cov[2,1]      <- sqrt(w.cov[1,1])*sqrt(w.cov[2,2])*as.numeric(v[6,3])
  w.cov[1,2]      <- sqrt(w.cov[1,1])*sqrt(w.cov[2,2])*as.numeric(v[6,3])
  colnames(w.cov) <- names
  rownames(w.cov) <- names
  # compile two matices in a list
  return(list(within=w.cov, between=b.cov))
}
