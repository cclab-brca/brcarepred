rm(list=ls())
library(mstate)
load("./Models/Pam50M.RData")
nat.death <- as.numeric((tra[,grep("NaturalDeath", colnames(tra))]))
nat.death <- unique(na.omit(nat.death))
i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(i)


missing <- unique(res[m$na.action,'id'])
Clinical <- Clinical[-which(Clinical$METABRIC.ID %in% missing),]
IC <- levels(Clinical$Group)
newdata <- NULL
for (ic in 1:length(IC)) {
    newdata <- rbind(newdata,
                     data.frame(strata=(9*(ic-1)+1):(9*(ic-1)+9),
                          AGE=Clinical$AGE[i],
                          LN=Clinical$LN[i],
                          GRADE=Clinical$GRADE[i],
                          SIZE=Clinical$SIZE[i],
                          TLastSurgery=c(rep(0, 4),
                          rep(Clinical$TLR[i], 3), rep(Clinical$TDR[i], 2))))
}

rel.can <- c(1:4, 7, 9, 10:13, 16, 18, 19:22, 25, 27, 28:31, 34, 36, 37:40,
             43, 45, 46:49, 52, 54, 55:58, 61, 63, 64:67, 70, 72, 73:76,
             79, 81, 82:85, 88, 90, 91:94, 97, 99)
newdata$AGE[which(!newdata$strata %in% nat.death)] <- 0
newdata$TLastSurgery[which(newdata$strata %in% rel.can)] <- 0

newdata$AGE.PS <- newdata$AGE * (newdata$strata %in%
                                         as.vector(tra[grep("Post", rownames(tra)),
                                                       grep("Nat", colnames(tra))]))
newdata$AGE.LR <- newdata$AGE * (newdata$strata %in%
                                         as.vector(tra[grep("Loco", rownames(tra)),
                                                       grep("Nat", colnames(tra))]))
newdata$AGE.DR <- newdata$AGE * (newdata$strata %in%
                                         as.vector(tra[grep("Distant", rownames(tra)),
                                                       grep("Nat", colnames(tra))]))

newdata$GRADE.PS <- newdata$GRADE * (newdata$strata %in%
                                         as.vector(tra[grep("Post", rownames(tra)),
                                                       -grep("Nat", colnames(tra))]))
newdata$GRADE.R <- newdata$GRADE * (newdata$strata %in%
                                        as.vector(tra[c(grep("Loco", rownames(tra)),
                                                        grep("Distant", rownames(tra))),
                                                       -grep("Nat", colnames(tra))]))

newdata$SIZE.PS <- newdata$SIZE * (newdata$strata %in%
                                         as.vector(tra[grep("Post", rownames(tra)),
                                                       -grep("Nat", colnames(tra))]))
newdata$SIZE.LR <- newdata$SIZE * (newdata$strata %in%
                                         as.vector(tra[grep("Loco", rownames(tra)),
                                                       -grep("Nat", colnames(tra))]))

newdata$LN.PS <- newdata$LN * (newdata$strata %in%
                                         as.vector(tra[grep("Post", rownames(tra)),
                                                       -grep("Nat", colnames(tra))]))
newdata$LN.R <- newdata$LN * (newdata$strata %in%
                                  as.vector(tra[c(grep("Loco", rownames(tra)),
                                                  grep("Distant", rownames(tra))),
                                                       -grep("Nat", colnames(tra))]))



class(newdata) <- c("msdata", "data.frame")
attr(newdata, "trans") <- tra
fitted <- msfit(m, newdata=newdata, trans=tra)
library(brcarepred)
tt <- seq(from=6, to=20, by=0.25)
pt <- list()
newdata.DR <- newdata
newdata.DR$AGE.DR[which(newdata.DR$strata %in% seq(from=9,by=9, length=5))] <-
    newdata.DR$AGE.DR[which(newdata.DR$strata %in% seq(from=9,by=9, length=5))] +
        newdata.DR$TLastSurgery[which(newdata.DR$strata %in%
                                       seq(from=8, by=9, length=5))]
system.time(pt[['DR']] <- getProbsDR(m, group=as.numeric(Clinical$Group[i]), newdata.DR, timepoints=tt))

newdata.LR <- newdata
newdata.LR$AGE.LR[which(newdata.LR$strata %in% seq(from=7, by=9, length=5))] <-
newdata.LR$AGE.LR[which(newdata.LR$strata %in% seq(from=7, by=9, length=5))] +
    newdata.LR$TLastSurgery[which(newdata.LR$strata %in%
                                   seq(from=6, by=9, length=5))]
newdata.LR$AGE.DR[which(newdata.LR$strata %in% seq(from=9,by=9, length=5))] <-
newdata.LR$AGE.LR[which(newdata.LR$strata %in% seq(from=7, by=9, length=5))]

system.time(pt[['LR']] <- getProbsLR(m, group=as.numeric(Clinical$Group[i]), newdata.LR, timepoints=tt, compact=FALSE))
tt <- c(seq(from=6, to=15, by=0.5), 16:20)
system.time(pt[['S']] <- getProbsS(m, group=as.numeric(Clinical$Group[i]), newdata, timepoints=tt, compact=FALSE))

save(pt, file=paste0("./IndivProbs/FROM5/Pam50TP_Patient", i, ".RData"))


