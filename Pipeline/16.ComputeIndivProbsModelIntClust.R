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

load(file="ICM.RData")
missing <- unique(res[m$na.action,'id'])
Clinical <- Clinical[-which(Clinical$METABRIC.ID %in% missing),]


ps <- list()
pslc <- list()
pslo <- list()
psld <- list()
psldc <- list()
psldo <- list()
psl <- list()
psdc <- list()
psdo <- list()
psd <- list()
psc <- list()
pso <- list()

pl <- list()
pld <- list()
plc <- list()
pldc <- list()
pldo <- list()
plo <- list()

pd <- list()
pdc <- list()
pdo <- list()


for (ic in levels(Clinical$Group)) {
ps[[ic]] <- NULL
pslc[[ic]] <- NULL
pslo[[ic]] <- NULL
psld[[ic]] <- NULL
psldc[[ic]] <- NULL
psldo[[ic]] <- NULL
psl[[ic]] <- NULL
psdc[[ic]] <- NULL
psdo[[ic]] <- NULL
psd[[ic]] <- NULL
psc[[ic]] <- NULL
pso[[ic]] <- NULL
pl[[ic]] <- NULL
pld[[ic]] <- NULL
plc[[ic]] <- NULL
pldc[[ic]] <- NULL
pldo[[ic]] <- NULL
plo[[ic]] <- NULL
pd[[ic]] <- NULL
pdc[[ic]] <- NULL
pdo[[ic]] <- NULL

x <- which(Clinical$Group==ic)
for (i in x) {
load(paste0("./IndivProbs/ICTP_Patient", i, ".RData"))

ps[[ic]] <- cbind(ps[[ic]], pt$S[['p.s']])
pslc[[ic]] <- cbind(pslc[[ic]], pt$S[['p.s.l.c']])
pslo[[ic]] <- cbind(pslo[[ic]], pt$S[['p.s.l.o']])
psld[[ic]] <- cbind(psld[[ic]], pt$S[['p.s.l.d']])
psldc[[ic]] <- cbind(psldc[[ic]], pt$S[['p.s.l.d.c']])
psldo[[ic]] <- cbind(psldo[[ic]], pt$S[['p.s.l.d.o']])
psl[[ic]] <- cbind(psl[[ic]], pt$S[['p.s.l']])
psdc[[ic]] <- cbind(psdc[[ic]], pt$S[['p.s.d.c']])
psdo[[ic]] <- cbind(psdo[[ic]], pt$S[['p.s.d.o']])
psd[[ic]] <- cbind(psd[[ic]], pt$S[['p.s.d']])
psc[[ic]] <- cbind(psc[[ic]], pt$S[['p.s.c']])
pso[[ic]] <- cbind(pso[[ic]], pt$S[['p.s.o']])

if (Clinical$DR[i]==1) {
pd[[ic]] <- cbind(pd[[ic]], pt$DR[['p.d']])
pdc[[ic]] <- cbind(pdc[[ic]], pt$DR[['p.d.c']])
pdo[[ic]] <- cbind(pdo[[ic]], pt$DR[['p.d.o']])
}

if (Clinical$LR[i]==1) {
pl[[ic]] <- cbind(pl[[ic]], pt$LR[['p.l']])
pld[[ic]] <- cbind(pld[[ic]], pt$LR[['p.l.d']])
plc[[ic]] <- cbind(plc[[ic]], pt$LR[['p.l.c']])
pldc[[ic]] <- cbind(pldc[[ic]], pt$LR[['p.l.d.c']])
pldo[[ic]] <- cbind(pldo[[ic]], pt$LR[['p.l.d.o']])
plo[[ic]] <- cbind(plo[[ic]], pt$LR[['p.l.o']])
}
}

ps[[ic]] <- cbind(pt$S[['Times']], ps[[ic]])
pslc[[ic]] <- cbind(pt$S[['Times']], pslc[[ic]])
pslo[[ic]] <- cbind(pt$S[['Times']], pslo[[ic]])
psld[[ic]] <- cbind(pt$S[['Times']], psld[[ic]])
psldc[[ic]] <- cbind(pt$S[['Times']], psldc[[ic]])
psldo[[ic]] <- cbind(pt$S[['Times']], psldo[[ic]])
psl[[ic]] <- cbind(pt$S[['Times']], psl[[ic]])
psdc[[ic]] <- cbind(pt$S[['Times']], psdc[[ic]])
psdo[[ic]] <- cbind(pt$S[['Times']], psdo[[ic]])
psd[[ic]] <- cbind(pt$S[['Times']], psd[[ic]])
psc[[ic]] <- cbind(pt$S[['Times']], psc[[ic]])
pso[[ic]] <- cbind(pt$S[['Times']], pso[[ic]])

pd[[ic]] <- cbind(pt$S[['Times']], pd[[ic]])
pdc[[ic]] <- cbind(pt$S[['Times']], pdc[[ic]])
pdo[[ic]] <- cbind(pt$S[['Times']], pdo[[ic]])

pl[[ic]] <- cbind(pt$S[['Times']], pl[[ic]])
pld[[ic]] <- cbind(pt$S[['Times']], pld[[ic]])
pldc[[ic]] <- cbind(pt$S[['Times']], pldc[[ic]])
pldo[[ic]] <- cbind(pt$S[['Times']], pldo[[ic]])
plc[[ic]] <- cbind(pt$S[['Times']], plc[[ic]])
plo[[ic]] <- cbind(pt$S[['Times']], plo[[ic]])

colnames(ps[[ic]]) <- c("Times", Clinical$METABRIC.ID[x])
colnames(pslc[[ic]]) <- c("Times", Clinical$METABRIC.ID[x])
colnames(pslo[[ic]]) <- c("Times", Clinical$METABRIC.ID[x])
colnames(psld[[ic]]) <- c("Times", Clinical$METABRIC.ID[x])
colnames(psldc[[ic]]) <- c("Times", Clinical$METABRIC.ID[x])
colnames(psldo[[ic]]) <- c("Times", Clinical$METABRIC.ID[x])
colnames(psl[[ic]]) <- c("Times", Clinical$METABRIC.ID[x])
colnames(psdc[[ic]]) <- c("Times", Clinical$METABRIC.ID[x])
colnames(psdo[[ic]]) <- c("Times", Clinical$METABRIC.ID[x])
colnames(psd[[ic]]) <- c("Times", Clinical$METABRIC.ID[x])
colnames(psc[[ic]]) <- c("Times", Clinical$METABRIC.ID[x])
colnames(pso[[ic]]) <- c("Times", Clinical$METABRIC.ID[x])

colnames(pd[[ic]]) <- c("Times", Clinical$METABRIC.ID[x][which(Clinical$DR[x]==1)])
colnames(pdc[[ic]]) <- c("Times", Clinical$METABRIC.ID[x][which(Clinical$DR[x]==1)])
colnames(pdo[[ic]]) <- c("Times", Clinical$METABRIC.ID[x][which(Clinical$DR[x]==1)])

colnames(pl[[ic]]) <- c("Times", Clinical$METABRIC.ID[x][which(Clinical$LR[x]==1)])
colnames(pld[[ic]]) <- c("Times", Clinical$METABRIC.ID[x][which(Clinical$LR[x]==1)])
colnames(plc[[ic]]) <- c("Times", Clinical$METABRIC.ID[x][which(Clinical$LR[x]==1)])
colnames(pldc[[ic]]) <- c("Times", Clinical$METABRIC.ID[x][which(Clinical$LR[x]==1)])
colnames(pldo[[ic]]) <- c("Times", Clinical$METABRIC.ID[x][which(Clinical$LR[x]==1)])
colnames(plo[[ic]]) <- c("Times", Clinical$METABRIC.ID[x][which(Clinical$LR[x]==1)])

}
save(ps, pslc, pslo, psld, psldc, psldo, psl, psdc, psdo, psd, psc, pso, pd, pdc, pdo, pl, pld, plc, pldc, pldo, plo, file="IntClust_AllProbs.RData")
