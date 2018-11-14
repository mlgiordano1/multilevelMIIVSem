#' Title
#'
#' Desc:
#'
#' Details to be written
#'
#' @param dat
#'
#' @export
# -----------------------------------------------------------------------------

alpha<-function(dat){
  covar<-dat[lower.tri(dat)] #get unique covariances
  n<-ncol(dat) #number of items in the scale
  al<-((sum(covar)/length(covar))*n^2/sum(dat))
  cat("alpha:",return(al),"\n")
}
