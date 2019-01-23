
rm(list=ls())
set.seed(534)
library(rms)

val <- list()
load("../Models/ERM.RData")
units(res$time) <- "Year"
m <- cph(Surv(time, status) ~ strat(trans) + AGE.PS.LR.NEG + AGE.DR.NEG +
             AGE.PS.LR.POS + AGE.DR.POS + GRADE.PS.NEG +
    GRADE.LR.NEG + GRADE.DR.NEG + GRADE.PS.POS + GRADE.LR.POS +
    GRADE.DR.POS + SIZE.PS.NEG + SIZE.LR.NEG + SIZE.DR.NEG +
    SIZE.PS.POS + SIZE.LR.POS + SIZE.DR.POS + LN.PS.NEG + LN.LR.NEG +
    LN.DR.NEG + LN.PS.POS + LN.LR.POS + LN.DR.POS + TLastSurgery.LR.NEG +
    TLastSurgery.DR.NEG + TLastSurgery.LR.POS + TLastSurgery.DR.POS +
        TLastLocal.DR.NEG + TLastLocal.DR.POS + cluster(id), data = res,
         x=T, y=T, surv=T)
val[['ER']] <- list()
for (i in c(5, 10, 15)) {
    val[['ER']][[as.character(i)]] <- validate(m, B=200, dxy=T, pr=F, u=i)
}
## save(val, cal, file="../Models/Internal_validation_ERM.RData")

load("../Models/FourGroupsM.RData")
units(res$time) <- "Year"
m <- cph(Surv(time, status) ~ strat(trans) + AGE.PS.LR + AGE.DR + GRADE.PS +
             GRADE.R + SIZE.PS + SIZE.LR + LN.PS + LN.R + TLastSurgery +
                cluster(id),   data = res,
         x=T, y=T, surv=T)
val[['FOURGROUPS']] <- list()
for (i in c(5, 10, 15)) {
    val[['FOURGROUPS']][[as.character(i)]] <- validate(m, B=200, dxy=T, pr=F, u=i)
}
## save(val, cal, file="../Models/Internal_validation_FOURGROUPSM.RData")

load("../Models/Pam50M.RData")
units(res$time) <- "Year"
m <- cph(Surv(time, status) ~ strat(trans) + AGE.PS.LR + AGE.DR + GRADE.PS +
             GRADE.R + SIZE.PS + SIZE.LR + LN.PS + LN.R + TLastSurgery +
                cluster(id),   data = res,
         x=T, y=T, surv=T)
val[['PAM50']] <- list()
for (i in c(5, 10, 15)) {
    val[['PAM50']][[as.character(i)]] <- validate(m, B=200, dxy=T, pr=F, u=i)
}

load("../Models/ICM.RData")
units(res$time) <- "Year"
m <- cph(Surv(time, status) ~ strat(trans) + AGE.PS.LR + AGE.DR + GRADE.PS +
             GRADE.R + SIZE.PS + SIZE.LR + LN.PS + LN.R + TLastSurgery +
                cluster(id),
         x=T, y=T, surv=T, data = res)
val[['INTCLUST']] <- list()
for (i in c(5, 10, 15)) {
    val[['INTCLUST']][[as.character(i)]] <- validate(m, B=200, dxy=T, pr=F, u=i)
}


save(val, file="../Validation/InternalValidation.RData")

rm(list=ls())
set.seed(43534)
cal <- list()

load("../Models/ERM.RData")
units(res$time) <- "Year"
cal[['ER']] <- list()
for (i in c(5, 10, 15)) {
    m2 <- cph(Surv(time, status) ~ strat(trans) + AGE.PS.LR.NEG + AGE.DR.NEG +
                  AGE.PS.LR.POS + AGE.DR.POS + GRADE.PS.NEG +
    GRADE.LR.NEG + GRADE.DR.NEG + GRADE.PS.POS + GRADE.LR.POS +
    GRADE.DR.POS + SIZE.PS.NEG + SIZE.LR.NEG + SIZE.DR.NEG +
    SIZE.PS.POS + SIZE.LR.POS + SIZE.DR.POS + LN.PS.NEG + LN.LR.NEG +
    LN.DR.NEG + LN.PS.POS + LN.LR.POS + LN.DR.POS + TLastSurgery.LR.NEG +
    TLastSurgery.DR.NEG + TLastSurgery.LR.POS + TLastSurgery.DR.POS +
        TLastLocal.DR.NEG + TLastLocal.DR.POS + cluster(id), data = res,
         x=T, y=T, surv=T,time.inc=i)
    cal[['ER']][[as.character(i)]] <- calibrate(m2, B=200, dxy=T, pr=F, u=i, m=200)
    attributes(cal[['ER']][[as.character(i)]])$predicted <-
        na.omit(attributes(cal[['ER']][[as.character(i)]])$predicted)
}

load("../Models/FourGroupsM.RData")
units(res$time) <- "Year"
cal[['FOURGROUPS']] <- list()
for (i in c(5, 10, 15)) {
    m2 <- cph(Surv(time, status) ~ strat(trans) + AGE.PS.LR + AGE.DR + GRADE.PS +
             GRADE.R + SIZE.PS + SIZE.LR + LN.PS + LN.R + TLastSurgery +
                cluster(id),   data = res,
              x=T, y=T, surv=T, time.inc=i)
    cal[['FOURGROUPS']][[as.character(i)]] <- calibrate(m2, B=200, dxy=T, pr=F, u=i, m=200)
    attributes(cal[['FOURGROUPS']][[as.character(i)]])$predicted <-
        na.omit(attributes(cal[['FOURGROUPS']][[as.character(i)]])$predicted)
}

load("../Models/Pam50M.RData")
units(res$time) <- "Year"
cal[['PAM50']] <- list()
for (i in c(5, 10, 15)) {
    m2 <- cph(Surv(time, status) ~ strat(trans) + AGE.PS.LR + AGE.DR + GRADE.PS +
             GRADE.R + SIZE.PS + SIZE.LR + LN.PS + LN.R + TLastSurgery +
                cluster(id),   data = res,
         x=T, y=T, surv=T, time.inc=i)
    cal[['PAM50']][[as.character(as.character(i))]] <- calibrate(m2, B=200,
                                                                 dxy=T, pr=F, u=i, m=200)
    attributes(cal[['PAM50']][[as.character(i)]])$predicted <-
        na.omit(attributes(cal[['PAM50']][[as.character(i)]])$predicted)
}

load("../Models/ICM.RData")
units(res$time) <- "Year"
cal[['INTCLUST']] <- list()
for (i in c(5, 10, 15)) {
    m2 <- cph(Surv(time, status) ~ strat(trans) + AGE.PS.LR + AGE.DR + GRADE.PS +
             GRADE.R + SIZE.PS + SIZE.LR + LN.PS + LN.R + TLastSurgery +
                cluster(id),
         x=T, y=T, surv=T, time.inc=i, data = res)
    cal[['INTCLUST']][[as.character(i)]] <- calibrate(m2, B=200, dxy=T, pr=F, u=i,  m=200)
    attributes(cal[['INTCLUST']][[as.character(i)]])$predicted <-
        na.omit(attributes(cal[['INTCLUST']][[as.character(i)]])$predicted)
}

save(cal, file="../Validation/InternalCalibration.RData")


rm(list=ls())


load("../Validation/InternalValidation.RData")

ER <- do.call("rbind", val[['ER']])
ER <- data.frame(Model="ER", Year=rep(c(5, 10, 15), c(7, 7, 7)),
                 Statistic=rownames(ER), ER)
FOURGROUPS <- do.call("rbind", val[['FOURGROUPS']])
FOURGROUPS <- data.frame(Model="CLINICAL", Year=rep(c(5, 10, 15), c(7, 7, 7)),
                 Statistic=rownames(FOURGROUPS), FOURGROUPS)
PAM <- do.call("rbind", val[['PAM50']])
PAM <- data.frame(Model="PAM50", Year=rep(c(5, 10, 15), c(7, 7, 7)),
                 Statistic=rownames(PAM), PAM)
INTCLUST <- do.call("rbind", val[['INTCLUST']])
INTCLUST <- data.frame(Model="INTCLUST", Year=rep(c(5, 10, 15), c(7, 7, 7)),
                 Statistic=rownames(INTCLUST), INTCLUST)
VAL <- rbind(ER, FOURGROUPS, PAM, INTCLUST)

write.table(VAL, file="../Validation/InternalValidation.txt", row.names=F, sep="\t", quote=F)


load("../Validation/InternalCalibration.RData")
pdf("../Validation/InternalCalibration.pdf", width=10, height=5)
names(cal)[2] <- "IHC"
for (i in 1:4) {
    par(mfrow=c(1, 3), mar=c(5, 2, 2, 2), oma=c(2, 2, 2, 2))
    for (j in 1:3) {
        plot(cal[[i]][[j]], main=paste("Model", names(cal)[i], names(cal[[i]][j]),
                                       "years"), cex.subtitles=0.5)
    }
}
dev.off()
