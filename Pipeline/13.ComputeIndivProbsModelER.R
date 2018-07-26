rm(list=ls())
library(mstate)
load("ERM.RData")
m <- fm
rm(fm)
nat.death <- as.numeric((tra[,grep("NaturalDeath", colnames(tra))]))
nat.death <- unique(na.omit(nat.death))
i <- as.numeric(Sys.getenv("LSB_JOBINDEX"))

if(!file.exists(paste0("./IndivProbs/ERTP_Patient", i, ".RData"))) {

missing <- unique(res[m$na.action,'id'])
Clinical <- Clinical[-which(Clinical$METABRIC.ID %in% missing),]
IC <- levels(Clinical$ER.Status)
newdata <- NULL
for (ic in 1:length(IC)) {
    newdata <- rbind(newdata,
                     data.frame(strata=(9*(ic-1)+1):(9*(ic-1)+9),
                          AGE=Clinical$AGE[i],
                          LN=Clinical$LN[i],
                          GRADE=Clinical$GRADE[i],
                          SIZE=Clinical$SIZE[i],
		          ER=Clinical$ER.Status[i],
 	                  TLastSurgery.LR=c(rep(0, 4),
                          rep(Clinical$TLR[i], 3), rep(0, 2)),
 	                  TLastSurgery.DR=c(rep(0, 4),
                          rep(0, 3), rep(Clinical$TDR[i], 2))
))
}

rel.can <- c(1:4, 7, 9, 10:13, 16, 18, 19:22, 25, 27, 28:31, 34, 36, 37:40,
             43, 45, 46:49, 52, 54, 55:58, 61, 63, 64:67, 70, 72, 73:76,
             79, 81, 82:85, 88, 90, 91:94, 97, 99)

newdata$AGE[which(!newdata$strata %in% nat.death)] <- 0
newdata$TLastSurgery.LR[which(newdata$strata %in% c(1:4, 7:13, 16:18))] <- 0
newdata$TLastSurgery.DR[which(newdata$strata %in% c(1:7, 9:16, 18))] <- 0
newdata$TLastSurgery.LR.NEG <- newdata$TLastSurgery.LR * (newdata$strata %in% 1:8)
newdata$TLastSurgery.LR.POS <- newdata$TLastSurgery.LR * (newdata$strata %in% 10:17)
newdata$TLastSurgery.DR.NEG <- newdata$TLastSurgery.DR * (newdata$strata %in% 1:8)
newdata$TLastSurgery.DR.POS <- newdata$TLastSurgery.DR * (newdata$strata %in% 10:17)

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


class(newdata) <- c("msdata", "data.frame")
attr(newdata, "trans") <- tra
fitted <- msfit(m, newdata=newdata, trans=tra)
library(brcarepred)
tt <- seq(from=0, to=20, by=0.25)
pt <- list()
system.time(pt[['DR']] <- getProbsDR(m, group=as.numeric(Clinical$ER.Status[i]), newdata, timepoints=tt))
if(Clinical$ER.Status[i]=="ER+") {
	beta <- coef(m)['TLastSurgery.LR.POS']
	x <- newdata[17,'TLastSurgery.DR.POS']
	LR<- newdata[15,'TLastSurgery.LR.POS']
}
if(Clinical$ER.Status[i]=="ER-") {
	beta <- coef(m)['TLastSurgery.LR.NEG']
	x <- newdata[8,'TLastSurgery.DR.NEG']
	LR<- newdata[6,'TLastSurgery.LR.NEG']
}
system.time(pt[['LR']] <- getProbsLR(m, group=as.numeric(Clinical$ER.Status[i]), newdata, timepoints=tt,
beta=beta, LR=LR, x=x, compact=FALSE))
tt <- c(2, 5, 10, 15, 20)
if(Clinical$ER.Status[i]=="ER+") {
	beta <- c(coef(m)['TLastSurgery.LR.POS'], coef(m)['TLastSurgery.DR.POS'])
	x <- c(newdata[15,'TLastSurgery.LR.POS'], newdata[17,'TLastSurgery.DR.POS'])
}
if(Clinical$ER.Status[i]=="ER-") {
	beta <- c(coef(m)['TLastSurgery.LR.NEG'], coef(m)['TLastSurgery.DR.NEG'])
	x <- c(newdata[6,'TLastSurgery.LR.NEG'], newdata[8,'TLastSurgery.DR.NEG'])
}
system.time(pt[['S']] <- getProbsS(m, group=as.numeric(Clinical$ER.Status[i]), newdata, timepoints=tt,
x=x, beta=beta, compact=FALSE))

save(pt, file=paste0("./IndivProbs/ERTP_Patient", i, ".RData"))

}



