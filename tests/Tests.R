try(setwd("C:/Users/mgiordan/Google Drive/PSYC/RESEARCH/masters/Empirical Examples"))
try(setwd("C:/Users/Michael/Google Drive/PSYC/RESEARCH/masters/Empirical Examples"))

library(multiLevelMIIVSem)

dat <- read.table("famIQData.dat")
names(dat) <- c("fam", "id", paste0("y", 1:6))

wModel <- '
l1 =~ y1 + y2 + y3
l2 =~ y4 + y5 + y6
l1~~l2
'

bModel <- '
Lb =~ y1 + y2 + y3 + y4 + y5 +y6
'

set.seed(20181108)
fit <- mlcfaMIIV(withinModel   = wModel,
                 betweenModel  = bModel,
                 estimator     = "Muthen",
                 allIndicators = paste0("y", 1:6),
                 l1Var         = "id",
                 l2Var         = "fam",
                 var.cov = TRUE,
                 df            = dat)
