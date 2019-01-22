rm(list=ls())
library(mstate)
library(brcarepredDEV)
load("./Models/ERM.RData")
m <- fm
## Predictions ER-

timepoints <- seq(from=0, to=20, by=0.25)
pt <- list()
newdata.DR <- newdata
newdata.DR$AGE.DR.NEG[which(newdata.DR$strata %in% c(9))] <-
    newdata.DR$AGE.DR.NEG[which(newdata.DR$strata %in% c(9))] +
    newdata.DR$TLastSurgery.DR.NEG[which(newdata.DR$strata %in% c(8))]
newdata.DR$AGE.DR.POS[which(newdata.DR$strata %in% c(18))] <-
    newdata.DR$AGE.DR.POS[which(newdata.DR$strata %in% c(18))] +
    newdata.DR$TLastSurgery.DR.POS[which(newdata.DR$strata %in% c(17))]

pt[['DR']] <- getProbsDR(fm, group=1, newdata.DR, timepoints=timepoints)
beta <- coef(m)['TLastSurgery.LR.NEG']
betaAGE <- coef(m)['AGE.DR.NEG']
DR <- newdata[8,'TLastSurgery.DR.NEG']
LR<- newdata[6,'TLastSurgery.LR.NEG']
LRAGE <- 0
DRAGE <-  newdata[9, 'AGE.DR.NEG'] - newdata[7, 'AGE.LR.NEG']
newdata.LR <- newdata
newdata.LR$AGE.PS.LR.NEG[which(newdata.LR$strata %in% c(7))] <-
    newdata.LR$AGE.PS.LR.NEG[which(newdata.LR$strata %in% c(7))] +
    newdata.LR$TLastSurgery.LR.NEG[which(newdata.LR$strata %in% c(6))]
newdata.LR$AGE.PS.LR.POS[which(newdata.LR$strata %in% c(16))] <-
    newdata.LR$AGE.PS.LR.POS[which(newdata.LR$strata %in% c(16))] +
    newdata.LR$TLastSurgery.LR.POS[which(newdata.LR$strata %in% c(15))]
newdata.LR$AGE.DR.NEG[which(newdata.LR$strata %in% c(9))] <-
newdata.LR$AGE.PS.LR.NEG[which(newdata.LR$strata %in% c(7))]
newdata.LR$AGE.DR.POS[which(newdata.LR$strata %in% c(18))] <-
newdata.LR$AGE.PS.LR.POS[which(newdata.LR$strata %in% c(16))]
pt[['LR']] <- getProbsLR(fm, group=1, newdata.LR, timepoints=timepoints,
 beta=beta, LR=LR, DR=DR, betaAGE=betaAGE, LRAGE=LRAGE, DRAGE=DRAGE, compact=FALSE)

timepoints <- c(seq(from=0, to=15, by=0.25), 16:20)
beta <- c(coef(m)['TLastSurgery.LR.NEG'], coef(m)['TLastSurgery.DR.NEG'])
betaAGE <- c(coef(m)['AGE.LR.NEG'], coef(m)['AGE.DR.NEG'])
x <- c(newdata[6,'TLastSurgery.LR.NEG'], newdata[8,'TLastSurgery.DR.NEG'])
xAGE <- c(x=newdata[7, 'AGE.LR.NEG'] - newdata[4, 'AGE.PS.NEG'],
newdata[9, 'AGE.DR.NEG'] - newdata[4, 'AGE.PS.NEG'])
system.time(pt[['S']] <- getProbsS(fm, group=1, newdata, timepoints=timepoints,
x=x, beta=beta, betaAGE=betaAGE, xAGE=xAGE, compact=FALSE))

save(pt, file="ERTP_NEG.RData")

################################################################
################################################################
################################################################
################################################################
## Predictions ER+

timepoints <- seq(from=0, to=20, by=0.25)
pt <- list()
newdata.DR <- newdata
newdata.DR$AGE.DR.NEG[which(newdata.DR$strata %in% c(9))] <-
    newdata.DR$AGE.DR.NEG[which(newdata.DR$strata %in% c(9))] +
    newdata.DR$TLastSurgery.DR.NEG[which(newdata.DR$strata %in% c(8))]
newdata.DR$AGE.DR.POS[which(newdata.DR$strata %in% c(18))] <-
    newdata.DR$AGE.DR.POS[which(newdata.DR$strata %in% c(18))] +
    newdata.DR$TLastSurgery.DR.POS[which(newdata.DR$strata %in% c(17))]

pt[['DR']] <- getProbsDR(fm, group=2, newdata.DR, timepoints=timepoints)

beta <- coef(m)['TLastSurgery.LR.POS']
betaAGE <- coef(m)['AGE.DR.POS']
DR <- newdata[17,'TLastSurgery.DR.POS']
LR<- newdata[15,'TLastSurgery.LR.POS']
LRAGE <- 0
DRAGE <-  newdata[18, 'AGE.DR.POS'] - newdata[16, 'AGE.LR.POS']

newdata.LR <- newdata
newdata.LR$AGE.PS.LR.NEG[which(newdata.LR$strata %in% c(7))] <-
    newdata.LR$AGE.PS.LR.NEG[which(newdata.LR$strata %in% c(7))] +
    newdata.LR$TLastSurgery.LR.NEG[which(newdata.LR$strata %in% c(6))]
newdata.LR$AGE.PS.LR.POS[which(newdata.LR$strata %in% c(16))] <-
    newdata.LR$AGE.PS.LR.POS[which(newdata.LR$strata %in% c(16))] +
    newdata.LR$TLastSurgery.LR.POS[which(newdata.LR$strata %in% c(15))]
newdata.LR$AGE.DR.NEG[which(newdata.LR$strata %in% c(9))] <-
newdata.LR$AGE.PS.LR.NEG[which(newdata.LR$strata %in% c(7))]
newdata.LR$AGE.DR.POS[which(newdata.LR$strata %in% c(18))] <-
newdata.LR$AGE.PS.LR.POS[which(newdata.LR$strata %in% c(16))]

pt[['LR']] <- getProbsLR(fm, group=2, newdata.LR, timepoints=timepoints,
 beta=beta, LR=LR, DR=DR, betaAGE=betaAGE, LRAGE=LRAGE, DRAGE=DRAGE, compact=FALSE)

timepoints <- c(seq(from=0, to=15, by=0.25), 16:20)
beta <- c(coef(m)['TLastSurgery.LR.POS'], coef(m)['TLastSurgery.DR.POS'])
betaAGE <- c(coef(m)['AGE.LR.POS'], coef(m)['AGE.DR.POS'])
x <- c(newdata[15,'TLastSurgery.LR.POS'], newdata[17,'TLastSurgery.DR.POS'])
xAGE <- c(x=newdata[16, 'AGE.LR.POS'] - newdata[13, 'AGE.PS.POS'],
newdata[18, 'AGE.DR.POS'] - newdata[13, 'AGE.PS.POS'])
system.time(pt[['S']] <- getProbsS(fm, group=2, newdata, timepoints=timepoints,
x=x, beta=beta, betaAGE=betaAGE, xAGE=xAGE, compact=FALSE))

save(pt, file="ERTP_POS.RData")
