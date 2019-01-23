rm(list=ls())
## Change group for the other Pam50 groups

load(file="./Models/Pam50M.RData")
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
id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
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
    newdata <- data.frame(strata=1:45,
                          AGE=x$AGE[i],
                          LN=x$LN[i],
                          GRADE=x$GRADE[i],
                          SIZE=x$SIZE[i],
                          TLastSurgery=0)
    tstate <- c(0,1,1,0,0)
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

        newdata$AGE[which(!newdata$strata %in% nat.death)] <- 0
newdata$AGE.PS <- newdata$AGE * (newdata$strata %in%
                                             as.vector(tra[grep("Post", rownames(tra)),
                                                           grep("Nat", colnames(tra))]))
newdata$AGE.LR <- newdata$AGE * (newdata$strata %in%
                                        as.vector(tra[c(grep("Loco", rownames(tra))),
                                                       grep("Nat", colnames(tra))]))
newdata$AGE.DR <- newdata$AGE * (newdata$strata %in%
                                        as.vector(tra[c(grep("Distant", rownames(tra))),
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

    state <- (as.numeric(x$Group[i])-1) * 5 + 1
    tstate <- c(rep(rep(0, 5), as.numeric(x$Group[i])-1),
                tstate, rep(rep(0, 5), 5-as.numeric(x$Group[i])))

        class(newdata) <- c("msdata", "data.frame")
        attr(newdata, "trans") <- tra
        beta.state <- matrix(0, 25, 45)
    beta.state[2,5:6] <- coef(m)['TLastSurgery']
    beta.state[7,14:15] <- coef(m)['TLastSurgery']
    beta.state[12,23:24] <- coef(m)['TLastSurgery']
    beta.state[17,32:33] <- coef(m)['TLastSurgery']
    beta.state[22,41:42] <- coef(m)['TLastSurgery']
    beta.state[2,7] <- coef(m)['AGE.LR']
    beta.state[7,16] <- coef(m)['AGE.LR']
    beta.state[12,25] <- coef(m)['AGE.LR']
    beta.state[17,34] <- coef(m)['AGE.LR']
    beta.state[22,43] <- coef(m)['AGE.LR']
    beta.state[3,8] <- coef(m)['TLastSurgery']
    beta.state[8,17] <- coef(m)['TLastSurgery']
    beta.state[13,26] <- coef(m)['TLastSurgery']
    beta.state[18,35] <- coef(m)['TLastSurgery']
    beta.state[23,44] <- coef(m)['TLastSurgery']
    beta.state[3,9] <- coef(m)['AGE.DR']
    beta.state[8,18] <- coef(m)['AGE.DR']
    beta.state[13,27] <- coef(m)['AGE.DR']
    beta.state[18,36] <- coef(m)['AGE.DR']
    beta.state[23,45] <- coef(m)['AGE.DR']

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
xb$TLastSurgery[which(xb$from %in% seq(from=2, by=5, length=5) &
                          xb$to %in% sort(c(seq(from=3, by=5, length=5),
                                            seq(from=4, by=5, length=5))))] <-
                                                xb$Tstart[which(xb$from %in% seq(from=2, by=5, length=5) &
                          xb$to %in% sort(c(seq(from=3, by=5, length=5),
                                            seq(from=4, by=5, length=5))))]
xb$TLastSurgery[which(xb$from %in% seq(from=3, by=5, length=5) &
                          xb$to %in% seq(from=4, by=5, length=5))] <-
                              xb$Tstart[which(xb$from %in% seq(from=3, by=5, length=5) &
                          xb$to %in% seq(from=4, by=5, length=5))]


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
    ## Adjust AGE with relapse
xb$AGE[which(xb$from %in% c(2,3,7,8, 12,13, 17,18,22,23) &
                 xb$to %in% c(5, 10, 15, 20,25))] <-
                     xb$AGE[which(xb$from %in% c(2,3,7,8,12,13,17,18,22,23) &
                                      xb$to %in% c(5, 10,15,20,25))] +
        xb$Tstart[which(xb$from %in% c(2,3,7,8,12,13,17,18,22,23) &
                            xb$to %in% c(5, 10,15,20,25))]

xb$AGE.PS <- xb$AGE * (xb$trans %in%
                               as.vector(tra[grep("Post", rownames(tra)),
                                             grep("Nat", colnames(tra))]))
xb$AGE.LR <- xb$AGE * (xb$trans %in%
                              as.vector(tra[c(grep("Loco", rownames(tra))),
                                            grep("Nat", colnames(tra))]))
xb$AGE.DR <- xb$AGE * (xb$trans %in%
                              as.vector(tra[c(grep("Distant", rownames(tra))),
                                            grep("Nat", colnames(tra))]))
colnames(xb)[4] <- "time"

class(xb) <- c("mstate", "data.frame")
attr(xb, "trans") <- tra
xb$Tstop[which(!is.finite(xb$Tstop))] <- max(x$T, na.rm=T)
xb$time <- xb$Tstop - xb$Tstart


mb <- coxph(Surv(time, status) ~ strata(trans) + AGE.PS + AGE.LR + AGE.DR + GRADE.PS + GRADE.R +                SIZE.PS + SIZE.LR +
                LN.PS + LN.R + TLastSurgery, data=xb)

library(brcarepred)

timepoints <- seq(from=0, to=20, by=0.25)
pt.boot <- list()
Oldnewdata.DR <- Oldnewdata
Oldnewdata.DR$AGE.DR[which(Oldnewdata.DR$strata %in% seq(from=9,by=9, length=11))] <-
    Oldnewdata.DR$AGE.DR[which(Oldnewdata.DR$strata %in% seq(from=9,by=9, length=11))] +
        Oldnewdata.DR$TLastSurgery[which(Oldnewdata.DR$strata %in%
                                       seq(from=8, by=9, length=11))]
pt.boot[['DR']] <- getProbsDR(mb, group=1, Oldnewdata.DR, timepoints=timepoints)

Oldnewdata.LR <- Oldnewdata
Oldnewdata.LR$AGE.LR[which(Oldnewdata.LR$strata %in% seq(from=7, by=9, length=11))] <-
Oldnewdata.LR$AGE.LR[which(Oldnewdata.LR$strata %in% seq(from=7, by=9, length=11))] +
    Oldnewdata.LR$TLastSurgery[which(Oldnewdata.LR$strata %in%
                                   seq(from=6, by=9, length=11))]
Oldnewdata.LR$AGE.DR[which(Oldnewdata.LR$strata %in% seq(from=9,by=9, length=11))] <-
Oldnewdata.LR$AGE.LR[which(Oldnewdata.LR$strata %in% seq(from=7, by=9, length=11))]

pt.boot[['LR']] <- getProbsLR(mb, group=1, Oldnewdata.LR, timepoints=timepoints, compact=F)

timepoints <- c(seq(from=0, to=15, by=0.25), 16:20)
system.time(pt.boot[['S']] <- getProbsS(mb, group=1, Oldnewdata, timepoints=timepoints, compact=FALSE))


save(mb, pt.boot, file=paste("./Bootstraps/BootsPredsPAM50MODEL_", id, "_GROUP1.RData", sep=""))
