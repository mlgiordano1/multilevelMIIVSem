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
  # ---------------------------------------------------------------------------
  # there has to be a better way to parse through the model for variables

  # ---------------------------------------------------------------------------
  # ---------------------------------------------------------------------------
  # check if all indicators matches the df
  if(!all(allIndicators %in% names(df))) {
    stop("allIndicators do no match variable names in DF")
  }
  if(!(l1Var %in% names(df))){
    stop("l1Var not in DF")
  }
  if(!(l2Var %in% names(df))){
    stop("l2Var not in DF")
  }
  if(!grepl(estimator, "muthen|goldstein", ignore.case = TRUE)){
    stop("estimator argument not recognized")
  }

  # -----------------------------------------------------------------
  # # this is lifted straight from MIIVsem for parsing the model...
  # # I need to go through and see if it works with my current setup
  # #-------------------------------------------------------#
  # # Check class of model.
  # #-------------------------------------------------------#
  # if ( "miivs" == class(model) ){
  #
  #   d  <- model$eqns
  #   pt <- model$pt
  #
  # } else {
  #
  #   res <- miivs(model)
  #   d   <- res$eqns
  #   pt  <- res$pt
  #
  # }
  #
  # #-------------------------------------------------------#
  # # parseInstrumentSyntax
  # #-------------------------------------------------------#
  # d  <- parseInstrumentSyntax(d, instruments, miiv.check)
  #
  # #-------------------------------------------------------#
  # # remove equations where there are not sufficient MIIVs
  # #-------------------------------------------------------#
  # underid   <- unlist(lapply(d, function(eq) {
  #   length(eq$MIIVs) < length(eq$IVobs)
  # }))
  # d.un <- d[underid]; d   <- d[!underid]
  #
  # #-------------------------------------------------------#
  # # Remove variables from data that are not in model
  # # syntax and preserve the original column ordering.
  # #-------------------------------------------------------#
  # if(!is.null(data)){
  #
  #   obs.vars <- unique(unlist(lapply(d,"[", c("DVobs", "IVobs", "MIIVs"))))
  #
  #   if (any(!obs.vars %in% colnames(data))){
  #    stop(paste(
  #      "miive: model syntax contains variables not in data.")
  #    )
  #   }
  #
  #   data <- data[,colnames(data) %in% obs.vars]
  #   data <- as.data.frame(data)
  #
  #   # convert any variables listed in
  #   # ordered to categorical
  #   data[,ordered]  <- lapply(
  #     data[,ordered, drop = FALSE], ordered
  #   )
  # }


  # -----------------------------------------------------------------
  # Trimming the df, saving the number of subjects and groups
  # -----------------------------------------------------------------
  df <- df[c(l2Var, l1Var, allIndicators)]
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
