## Change group for ER + instead of ER-

load(file="ERM.RData")
Oldnewdata <- newdata
source("mssampleOscar.R")
library(mstate)
tmp <- msfit(fm, newdata=Oldnewdata, trans=tra)
tmp <- tmp$Haz
tmp <- split(tmp, tmp$trans)
Times <- tmp[[1]]$time
m <- fm
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
    newdata <- data.frame(strata=1:18,
                          AGE=x$AGE[i],
                          LN=x$LN[i],
                          GRADE=x$GRADE[i],
                          SIZE=x$SIZE[i],
                          TLastSurgery.LR.NEG=0,
                          TLastSurgery.LR.POS=0,
                          TLastSurgery.DR.NEG=0,
                          TLastSurgery.DR.POS=0)
    tstate <- c(0,0,0,0,0)
    if (x$LR[i]==1) {
        newdata$TLastSurgery.LR.NEG[c(5,6)] <- x$TLR[i]
        newdata$TLastSurgery.LR.POS[c(14,15)] <- x$TLR[i]
        tstate[2] <- x$TLR[i]
    }
    if (x$DR[i]==1) {
        newdata$TLastSurgery.DR.NEG[c(8)] <- x$TDR[i]
        newdata$TLastSurgery.DR.POS[c(17)] <- x$TDR[i]
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

        rel.can <- c(1:4, 10:13)
        newdata$AGE[which(!newdata$strata %in% nat.death)] <- 0

        newdata$GRADE.PS.NEG <- newdata$GRADE * (newdata$strata %in% 1:3)
        newdata$GRADE.PS.POS <- newdata$GRADE * (newdata$strata %in% 10:12)
        newdata$GRADE.LR.NEG <- newdata$GRADE * (newdata$strata %in% 5:6)
        newdata$GRADE.LR.POS <- newdata$GRADE * (newdata$strata %in% 14:15)
        newdata$GRADE.DR.NEG <- newdata$GRADE * (newdata$strata %in% 8)
        newdata$GRADE.DR.POS <- newdata$GRADE * (newdata$strata %in% 17)

        newdata$SIZE.PS.NEG <- newdata$SIZE * (newdata$strata %in% 1:3)
        newdata$SIZE.PS.POS <- newdata$SIZE * (newdata$strata %in% 10:12)
        newdata$SIZE.LR.NEG <- newdata$SIZE * (newdata$strata %in% 5:6)
        newdata$SIZE.LR.POS <- newdata$SIZE * (newdata$strata %in% 14:15)
        newdata$SIZE.DR.NEG <- newdata$SIZE * (newdata$strata %in% 8)
        newdata$SIZE.DR.POS <- newdata$SIZE * (newdata$strata %in% 17)

        newdata$LN.PS.NEG <- newdata$LN * (newdata$strata %in% 1:3)
        newdata$LN.PS.POS <- newdata$LN * (newdata$strata %in% 10:12)
        newdata$LN.LR.NEG <- newdata$LN * (newdata$strata %in% 5:6)
        newdata$LN.LR.POS <- newdata$LN * (newdata$strata %in% 14:15)
        newdata$LN.DR.NEG <- newdata$LN * (newdata$strata %in% 8)
        newdata$LN.DR.POS <- newdata$LN * (newdata$strata %in% 17)


    if (x$ER[i]=="ER-") {
            state <- 1
            tstate <- c(tstate, rep(0, 5))
        } else {
            state <- 6
            tstate <- c(rep(0, 5), tstate)
        }

        class(newdata) <- c("msdata", "data.frame")
        attr(newdata, "trans") <- tra
        beta.state <- matrix(0, 10, 18)
        beta.state[2,5:6] <- coef(m)['TLastSurgery.LR.NEG']
        beta.state[3,8] <- coef(m)['TLastSurgery.DR.NEG']
        beta.state[7,14:15] <- coef(m)['TLastSurgery.LR.POS']
        beta.state[8,17] <- coef(m)['TLastSurgery.DR.POS']
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
                         'SIZE', 'LN', 'ER.Status')], by.x="id", by.y="METABRIC.ID")

    nat.death <- as.numeric((tra[,grep("NaturalDeath", colnames(tra))]))
    nat.death <- unique(na.omit(nat.death))

    xb$AGE[which(!xb$trans %in% nat.death)] <- 0
    xb$GRADE[which(xb$trans %in% nat.death)] <- 0
    xb$SIZE[which(xb$trans %in% nat.death)] <- 0
    xb$LN[which(xb$trans %in% nat.death)] <- 0
xb$TLastSurgery <- 0
xb$TLastSurgery[which(xb$from %in% c(2,7) & xb$to %in% c(3,4,8,9))] <- xb$Tstart[which(xb$from %in% c(2,7) & xb$to %in% c(3,4,8,9))]
xb$TLastSurgery[which(xb$from %in% c(3,8) & xb$to %in% c(4,9))] <- xb$Tstart[which(xb$from %in% c(3,8) & xb$to %in% c(4, 9))]
xb$TLastSurgery.NEG <- xb$TLastSurgery * (xb$ER.Status=="ER-")
xb$TLastSurgery.POS <- xb$TLastSurgery * (xb$ER.Status=="ER+")

    xb$GRADE.NEG <- xb$GRADE * (xb$ER.Status=="ER-")
    xb$GRADE.POS <- xb$GRADE * (xb$ER.Status=="ER+")
    xb$SIZE.NEG <- xb$SIZE * (xb$ER.Status=="ER-")
    xb$SIZE.POS <- xb$SIZE * (xb$ER.Status=="ER+")
    xb$LN.NEG <- xb$LN * (xb$ER.Status=="ER-")
    xb$LN.POS <- xb$LN * (xb$ER.Status=="ER+")

    xb$GRADE.PS <- xb$GRADE * (xb$from %in% grep("PostSurgery", colnames(tra)))
    xb$GRADE.LR <- xb$GRADE * (xb$from %in% grep("LocalRelapse", colnames(tra)))
    xb$GRADE.DR <- xb$GRADE * (xb$from %in% grep("DistantRelapse", colnames(tra)))
    xb$SIZE.PS <- xb$SIZE * (xb$from %in% grep("PostSurgery", colnames(tra)))
    xb$SIZE.LR <- xb$SIZE * (xb$from %in% grep("LocalRelapse", colnames(tra)))
    xb$SIZE.DR <- xb$SIZE * (xb$from %in% grep("DistantRelapse", colnames(tra)))
    xb$LN.PS <- xb$LN * (xb$from %in% grep("PostSurgery", colnames(tra)))
    xb$LN.LR <- xb$LN * (xb$from %in% grep("LocalRelapse", colnames(tra)))
    xb$LN.DR <- xb$LN * (xb$from %in% grep("DistantRelapse", colnames(tra)))
xb$TLastSurgery.LR.NEG <- xb$TLastSurgery.NEG * (xb$from %in% grep("LocalRelapse", colnames(tra)))
xb$TLastSurgery.LR.POS <- xb$TLastSurgery.POS * (xb$from %in% grep("LocalRelapse", colnames(tra)))
xb$TLastSurgery.DR.NEG <- xb$TLastSurgery.NEG * (xb$from %in% grep("DistantRelapse", colnames(tra)))
xb$TLastSurgery.DR.POS <- xb$TLastSurgery.POS * (xb$from %in% grep("DistantRelapse", colnames(tra)))



    xb$GRADE.PS.NEG <- xb$GRADE.NEG * (xb$from %in% grep("PostSurgery", colnames(tra)))
    xb$GRADE.PS.POS <- xb$GRADE.POS * (xb$from %in% grep("PostSurgery", colnames(tra)))
    xb$GRADE.LR.NEG <- xb$GRADE.NEG * (xb$from %in% grep("LocalRelapse", colnames(tra)))
    xb$GRADE.LR.POS <- xb$GRADE.POS * (xb$from %in% grep("LocalRelapse", colnames(tra)))
    xb$GRADE.DR.NEG <- xb$GRADE.NEG * (xb$from %in% grep("DistantRelapse", colnames(tra)))
    xb$GRADE.DR.POS <- xb$GRADE.POS * (xb$from %in% grep("DistantRelapse", colnames(tra)))
    xb$SIZE.PS.NEG <- xb$SIZE.NEG * (xb$from %in% grep("PostSurgery", colnames(tra)))
    xb$SIZE.PS.POS <- xb$SIZE.POS * (xb$from %in% grep("PostSurgery", colnames(tra)))
    xb$SIZE.LR.NEG <- xb$SIZE.NEG * (xb$from %in% grep("LocalRelapse", colnames(tra)))
    xb$SIZE.LR.POS <- xb$SIZE.POS * (xb$from %in% grep("LocalRelapse", colnames(tra)))
    xb$SIZE.DR.NEG <- xb$SIZE.NEG * (xb$from %in% grep("DistantRelapse", colnames(tra)))
    xb$SIZE.DR.POS <- xb$SIZE.POS * (xb$from %in% grep("DistantRelapse", colnames(tra)))
    xb$LN.PS.NEG <- xb$LN.NEG * (xb$from %in% grep("PostSurgery", colnames(tra)))
    xb$LN.PS.POS <- xb$LN.POS * (xb$from %in% grep("PostSurgery", colnames(tra)))
    xb$LN.LR.NEG <- xb$LN.NEG * (xb$from %in% grep("LocalRelapse", colnames(tra)))
    xb$LN.LR.POS <- xb$LN.POS * (xb$from %in% grep("LocalRelapse", colnames(tra)))
    xb$LN.DR.NEG <- xb$LN.NEG * (xb$from %in% grep("DistantRelapse", colnames(tra)))
    xb$LN.DR.POS <- xb$LN.POS * (xb$from %in% grep("DistantRelapse", colnames(tra)))

    colnames(xb)[4] <- "time"

    class(xb) <- c("mstate", "data.frame")
    attr(xb, "trans") <- tra
    xb$Tstop[which(!is.finite(xb$Tstop))] <- max(x$T, na.rm=T)
    xb$time <- xb$Tstop - xb$Tstart


mb <- coxph(Surv(time, status) ~ strata(trans) + AGE + GRADE.PS.NEG + GRADE.LR.NEG + GRADE.DR.NEG + GRADE.PS.POS + GRADE.LR.POS + GRADE.DR.POS +
                SIZE.PS.NEG + SIZE.LR.NEG + SIZE.DR.NEG + SIZE.PS.POS + SIZE.LR.POS + SIZE.DR.POS +
                    LN.PS.NEG + LN.LR.NEG + LN.DR.NEG + LN.PS.POS + LN.LR.POS + LN.DR.POS +
                        TLastSurgery.LR.NEG + TLastSurgery.DR.NEG +
                            TLastSurgery.LR.POS + TLastSurgery.DR.POS,
                          data=xb)
library(brcarepred)

timepoints <- seq(from=0, to=20, by=0.25)
pt.boot <- list()
pt.boot[['DR']] <- getProbsDR(mb, group=1, Oldnewdata, timepoints=timepoints)
beta <- coef(mb)['TLastSurgery.LR.NEG']
x <- Oldnewdata[8,'TLastSurgery.DR.NEG']
LR<- Oldnewdata[6,'TLastSurgery.LR.NEG']

pt.boot[['LR']] <- getProbsLR(mb, group=1, Oldnewdata, timepoints=timepoints, beta=beta, LR=LR, x=x, compact=F)

timepoints <- c(seq(from=0, to=15, by=0.25), 16:20)
beta <- c(coef(mb)['TLastSurgery.LR.NEG'], coef(mb)['TLastSurgery.DR.NEG'])
x <- c(Oldnewdata[6,'TLastSurgery.LR.NEG'], Oldnewdata[8,'TLastSurgery.DR.NEG'])

system.time(pt.boot[['S']] <- getProbsS(mb, group=1, Oldnewdata, timepoints=timepoints, x=x, beta=beta, compact=FALSE))


save(mb, pt.boot, file=paste("./Bootstraps/BootsPredsERMODEL_", id, "_NEG.RData", sep=""))
