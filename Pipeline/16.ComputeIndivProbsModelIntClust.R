rm(list=ls())

library(mstate)
load("ICM.RData")
nat.death <- as.numeric((tra[,grep("NaturalDeath", colnames(tra))]))
nat.death <- unique(na.omit(nat.death))
i <- as.numeric(Sys.getenv("LSB_JOBINDEX"))

if(!file.exists(paste0("./IndivProbs/ICTP_Patient", i, ".RData"))) {

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
tt <- c(seq(from=0, by=0.25), c(16:20))
pt <- list()
system.time(pt[['DR']] <- getProbsDR(m, group=as.numeric(Clinical$Group[i]), newdata, timepoints=tt))
system.time(pt[['LR']] <- getProbsLR(m, group=as.numeric(Clinical$Group[i]), newdata, timepoints=tt, compact=FALSE))
tt <- c(2, 5, 10, 15, 20)
system.time(pt[['S']] <- getProbsS(m, group=as.numeric(Clinical$Group[i]), newdata, timepoints=tt, compact=FALSE))

save(pt, file=paste0("./IndivProbs/ICTP_Patient", i, ".RData"))


}

