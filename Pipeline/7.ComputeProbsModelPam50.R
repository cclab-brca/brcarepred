rm(list=ls())
id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(id)
library(mstate)
library(brcarepred)
load("../Models/Pam50M.RData")

## Predictions

timepoints <- seq(from=0, to=20, by=0.25)
pt <- list()
newdata.DR <- newdata
newdata.DR$AGE.DR[which(newdata.DR$strata %in% seq(from=9,by=9, length=5))] <-
    newdata.DR$AGE.DR[which(newdata.DR$strata %in% seq(from=9,by=9, length=5))] +
        newdata.DR$TLastSurgery[which(newdata.DR$strata %in%
                                       seq(from=8, by=9, length=5))]

pt[['DR']] <- getProbsDR(m, group=id, newdata.DR, timepoints=timepoints)
newdata.LR <- newdata
newdata.LR$AGE.LR[which(newdata.LR$strata %in% seq(from=7, by=9, length=5))] <-
newdata.LR$AGE.LR[which(newdata.LR$strata %in% seq(from=7, by=9, length=5))] +
    newdata.LR$TLastSurgery[which(newdata.LR$strata %in%
                                   seq(from=6, by=9, length=5))]
newdata.LR$AGE.DR[which(newdata.LR$strata %in% seq(from=9,by=9, length=5))] <-
newdata.LR$AGE.LR[which(newdata.LR$strata %in% seq(from=7, by=9, length=5))]

pt[['LR']] <- getProbsLR(m, group=id, newdata.LR, timepoints=timepoints, compact=F)

timepoints <- c(seq(from=0, to=15, by=0.25), 16:20)
system.time(pt[['S']] <- getProbsS(m, group=id, newdata, timepoints=timepoints, compact=FALSE))

save(pt, file=paste0("TP_PAM50_", id, ".RData"))
