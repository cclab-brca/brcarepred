id <- as.numeric(Sys.getenv("LSB_JOBINDEX"))
library(mstate)
library(brcarepred)
load("ICM.RData")
fitted <- msfit(m, newdata=newdata, trans=tra)

## Predictions

timepoints <- seq(from=0, to=20, by=0.25)
pt <- list()
pt[['DR']] <- getProbsDR(m, group=id, newdata, timepoints=timepoints)

pt[['LR']] <- getProbsLR(m, group=id, newdata, timepoints=timepoints, compact=F)

timepoints <- c(seq(from=0, to=15, by=0.25), 16:20)
system.time(pt[['S']] <- getProbsS(m, group=id, newdata, timepoints=timepoints, compact=FALSE))

save(pt, file=paste0("TP_IntCLUST_", id, ".RData"))
