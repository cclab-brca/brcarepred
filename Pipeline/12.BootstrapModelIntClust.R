rm(list=ls())
## Change group for the other IntClust

load(file="ICM.RData")
Oldnewdata <- newdata
source("mssampleOscar.R")
library(mstate)
tmp <- msfit(m, newdata=Oldnewdata, trans=tra)
tmp <- tmp$Haz
tmp <- split(tmp, tmp$trans)
Times <- tmp[[1]]$time
x <- Clinical
nat.death <- as.numeric((tra[,grep("NaturalDeath", colnames(tra))]))
nat.death <- unique(na.omit(nat.death))
id <- as.numeric(Sys.getenv("LSB_JOBINDEX"))
set.seed(3 * id)
if (length(m$na.action) > 0) {
        good.ids <- as.character(unique(res[-(m$na.action),'id']))
        x <- x[which(x$METABRIC.ID %in% good.ids),]
    }
## Bootstrap datasets
xb <- NULL
## indexes <- sample(1:nrow(x), nrow(x), replace=T)
indexes <- 1:nrow(x)
indexes <- table(indexes)
indexes <- data.frame(indexes)
counter <- 0
for (i in as.numeric(as.character(indexes$indexes))) {
    newdata <- data.frame(strata=1:99,
                          AGE=x$AGE[i],
                          LN=x$LN[i],
                          GRADE=x$GRADE[i],
                          SIZE=x$SIZE[i],
                          TLastSurgery=0)
    tstate <- c(0,0,0,0,0)
    if (x$LR[i]==1) {
        newdata$TLastSurgery[sort(c(seq(from=5, length=11, by=9),
                                 seq(from=6, length=11, by=9)))] <- x$TLR[i]
        tstate[2] <- x$TLR[i]
    }
    if (x$DR[i]==1) {
        newdata$TLastSurgery[seq(from=8, length=11, by=9)] <- x$TDR[i]
        tstate[3] <- x$TDR[i]
    }
    if (any(is.na(newdata$AGE))) newdata$AGE <-
            mean(x$AGE, na.rm=T)
        if (any(is.na(newdata$LN))) newdata$LN <-
            mean(x$LN, na.rm=T)
        if (any(is.na(newdata$GRADE))) {
            newdata$GRADE <-
                names(which.max(table(x$GRADE)))
        }
        newdata$GRADE <- as.numeric(as.character(newdata$GRADE))
        if (any(is.na(newdata$SIZE))) newdata$SIZE <-
            mean(x$SIZE, na.rm=T)
        if (any(is.na(newdata$BS))) newdata$BS <-
            names(which.max(table(x$BS)))

        newdata$AGE[which(!newdata$strata %in% nat.death)] <- 0
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


    state <- (as.numeric(x$Group[i])-1) * 5 + 1
    tstate <- c(rep(rep(0, 5), as.numeric(x$Group[i])-1),
                tstate, rep(rep(0, 5), 11-as.numeric(x$Group[i])))
    class(newdata) <- c("msdata", "data.frame")
    attr(newdata, "trans") <- tra
    beta.state <- matrix(0, 55, 99)
    beta.state[2,5:6] <- coef(m)['TLastSurgery']
    beta.state[7,14:15] <- coef(m)['TLastSurgery']
    beta.state[12,23:24] <- coef(m)['TLastSurgery']
    beta.state[17,32:33] <- coef(m)['TLastSurgery']
    beta.state[22,41:42] <- coef(m)['TLastSurgery']
    beta.state[27,50:51] <- coef(m)['TLastSurgery']
    beta.state[32,59:60] <- coef(m)['TLastSurgery']
    beta.state[37,68:69] <- coef(m)['TLastSurgery']
    beta.state[42,77:78] <- coef(m)['TLastSurgery']
    beta.state[47,86:87] <- coef(m)['TLastSurgery']
    beta.state[52,95:96] <- coef(m)['TLastSurgery']

    beta.state[3,8] <- coef(m)['TLastSurgery']
    beta.state[8,17] <- coef(m)['TLastSurgery']
    beta.state[13,26] <- coef(m)['TLastSurgery']
    beta.state[18,35] <- coef(m)['TLastSurgery']
    beta.state[23,44] <- coef(m)['TLastSurgery']
    beta.state[28,53] <- coef(m)['TLastSurgery']
    beta.state[33,62] <- coef(m)['TLastSurgery']
    beta.state[38,71] <- coef(m)['TLastSurgery']
    beta.state[43,80] <- coef(m)['TLastSurgery']
    beta.state[48,89] <- coef(m)['TLastSurgery']
    beta.state[53,98] <- coef(m)['TLastSurgery']


    beta.state[seq(from=3, length=11, by=5),seq(from=8, length=11, by=9)] <-
        coef(m)['TLastSurgery']
    fitted.Haz <- msfit(m, newdata=newdata, trans=tra)
    tmp <- try(oscar.mssample(fitted.Haz$Haz, trans=tra, clock="reset",
                              output="data",
                              history=list(state=state,
                                  time=0, tstate=tstate),
                              beta.state=beta.state,
                              M=1, cens=cens,
                              tvec=seq(from=0, to=20, length=100)))
    if (class(tmp)!="try-error") {
            class(tmp) <- c("msdata", "data.frame")
            tmp$id <- x$METABRIC.ID[i]
            xb <- rbind(xb, tmp)
        }
    }

    xb <- merge(xb, x[,c('METABRIC.ID', 'AGE', 'GRADE',
                         'SIZE', 'LN', 'Group')], by.x="id", by.y="METABRIC.ID")

nat.death <- as.numeric((tra[,grep("NaturalDeath", colnames(tra))]))
    nat.death <- unique(na.omit(nat.death))

    xb$AGE[which(!xb$trans %in% nat.death)] <- 0
    xb$GRADE[which(xb$trans %in% nat.death)] <- 0
    xb$SIZE[which(xb$trans %in% nat.death)] <- 0
    xb$LN[which(xb$trans %in% nat.death)] <- 0
xb$TLastSurgery <- 0
xb$TLastSurgery[which(xb$from %in% seq(from=2, by=5, length=11) &
                          xb$to %in% sort(c(seq(from=3, by=5, length=11),
                                            seq(from=4, by=5, length=11))))] <-
                                                xb$Tstart[which(xb$from %in% seq(from=2, by=5, length=11) &
                          xb$to %in% sort(c(seq(from=3, by=5, length=11),
                                            seq(from=4, by=5, length=11))))]
xb$TLastSurgery[which(xb$from %in% seq(from=3, by=5, length=11) &
                          xb$to %in% seq(from=4, by=5, length=11))] <-
                              xb$Tstart[which(xb$from %in% seq(from=3, by=5, length=11) &
                          xb$to %in% seq(from=4, by=5, length=11))]

xb$GRADE.PS <- xb$GRADE * (xb$trans %in%
                               as.vector(tra[grep("Post", rownames(tra)),
                                             -grep("Nat", colnames(tra))]))
xb$GRADE.R <- xb$GRADE * (xb$trans %in%
                              as.vector(tra[c(grep("Loco", rownames(tra)),
                                              grep("Distant", rownames(tra))),
                                            -grep("Nat", colnames(tra))]))

xb$SIZE.PS <- xb$SIZE * (xb$trans %in%
                             as.vector(tra[grep("Post", rownames(tra)),
                                           -grep("Nat", colnames(tra))]))
xb$SIZE.LR <- xb$SIZE * (xb$trans %in%
                             as.vector(tra[grep("Loco", rownames(tra)),
                                           -grep("Nat", colnames(tra))]))

xb$LN.PS <- xb$LN * (xb$trans %in%
                         as.vector(tra[grep("Post", rownames(tra)),
                                       -grep("Nat", colnames(tra))]))
xb$LN.R <- xb$LN * (xb$trans %in%
                        as.vector(tra[c(grep("Loco", rownames(tra)),
                                        grep("Distant", rownames(tra))),
                                      -grep("Nat", colnames(tra))]))
    colnames(xb)[4] <- "time"

class(xb) <- c("mstate", "data.frame")
attr(xb, "trans") <- tra
xb$Tstop[which(!is.finite(xb$Tstop))] <- max(x$T, na.rm=T)
xb$time <- xb$Tstop - xb$Tstart


mb <- coxph(Surv(time, status) ~ strata(trans) + AGE + GRADE.PS + GRADE.R +                SIZE.PS + SIZE.LR +
                LN.PS + LN.R + TLastSurgery, data=xb)

library(brcarepred)

timepoints <- seq(from=0, to=20, by=0.25)
pt.boot <- list()
pt.boot[['DR']] <- getProbsDR(mb, group=1, Oldnewdata, timepoints=timepoints)

pt.boot[['LR']] <- getProbsLR(mb, group=1, Oldnewdata, timepoints=timepoints, compact=F)

timepoints <- c(seq(from=0, to=15, by=0.25), 16:20)
system.time(pt.boot[['S']] <- getProbsS(mb, group=1, Oldnewdata, timepoints=timepoints, compact=FALSE))


save(mb, pt.boot, file=paste("./Bootstraps/BootsPredsINTCLUSTMODEL_", id, "_GROUP1.RData", sep=""))
