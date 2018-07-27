rm(list=ls())
library(mstate)
load("../Models/ERM.RData")
res <- events(res)
res$Frequencies
rm(list=ls())
load("../Models/FourGroupsM.RData")
res <- events(res)
res$Frequencies
load("../Models/Pam50M.RData")
res <- events(res)
res$Frequencies
load("../Models/ICM.RData")
res <- events(res)
res$Frequencies



