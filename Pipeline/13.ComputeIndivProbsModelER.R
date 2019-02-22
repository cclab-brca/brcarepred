 rm(list=ls())
 library(mstate)
 load("./Models/ERM.RData")
 m <- fm
 rm(fm)
 nat.death <- as.numeric((tra[,grep("NaturalDeath", colnames(tra))]))
 nat.death <- unique(na.omit(nat.death))
 i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

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

 newdata$AGE.PS <- newdata$AGE * (newdata$strata %in%
                                          as.vector(tra[grep("Post", rownames(tra)),
                                                        grep("Nat", colnames(tra))]))
 newdata$AGE.LR <- newdata$AGE * (newdata$strata %in%
                                          as.vector(tra[grep("Loc", rownames(tra)),
                                                        grep("Nat", colnames(tra))]))
 newdata$AGE.DR <- newdata$AGE * (newdata$strata %in%
                                          as.vector(tra[grep("Distant", rownames(tra)),
                                                        grep("Nat", colnames(tra))]))

 newdata$AGE.PS.NEG <- newdata$AGE.PS * (newdata$strata %in% c(4))
 newdata$AGE.PS.POS <- newdata$AGE.PS * (newdata$strata %in% c(13))
 newdata$AGE.LR.NEG <- newdata$AGE.LR * (newdata$strata %in% c(7))
 newdata$AGE.LR.POS <- newdata$AGE.LR * (newdata$strata %in% c(16))
 newdata$AGE.DR.NEG <- newdata$AGE.DR * (newdata$strata %in% c(9))
 newdata$AGE.DR.POS <- newdata$AGE.DR * (newdata$strata %in% c(18))


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
 newdata.DR <- newdata
 newdata.DR$AGE.DR.NEG[which(newdata.DR$strata==9)] <-
     newdata.DR$AGE.DR.NEG[which(newdata.DR$strata==9)] +
         newdata.DR$TLastSurgery.DR.NEG[which(newdata.DR$strata==8)]
 newdata.DR$AGE.DR.POS[which(newdata.DR$strata==18)] <-
     newdata.DR$AGE.DR.POS[which(newdata.DR$strata==18)] +
         newdata.DR$TLastSurgery.DR.POS[which(newdata.DR$strata==17)]
 system.time(pt[['DR']] <- getProbsDR(m, group=as.numeric(Clinical$ER.Status[i]), newdata.DR, timepoints=tt))
 if(Clinical$ER.Status[i]=="ER+") {
 				   beta <- coef(m)['TLastSurgery.DR.POS']
 				   betaAGE <- coef(m)['AGE.DR.POS']
 				   DR <- newdata[17,'TLastSurgery.DR.POS']
 				   LR<- newdata[15,'TLastSurgery.LR.POS']
 				   LRAGE <- 0
 				   DRAGE <-  newdata[18, 'AGE.DR.POS'] - newdata[16, 'AGE.LR.POS']

 }
 if(Clinical$ER.Status[i]=="ER-") {
 				   beta <- coef(m)['TLastSurgery.DR.NEG']
 				   betaAGE <- coef(m)['AGE.DR.NEG']
 				   DR <- newdata[8,'TLastSurgery.DR.NEG']
 				   LR<- newdata[6,'TLastSurgery.LR.NEG']
 				   LRAGE <- 0
 				   DRAGE <-  newdata[9, 'AGE.DR.NEG'] - newdata[7, 'AGE.LR.NEG']
 }
 newdata.LR <- newdata
 newdata.LR$AGE.LR.NEG[which(newdata.LR$strata==7)] <-
 newdata.LR$AGE.LR.NEG[which(newdata.LR$strata==7)] 	 +
   newdata.LR$TLastSurgery.LR.NEG[which(newdata.LR$strata == 6)]
 newdata.LR$AGE.LR.POS[which(newdata.LR$strata==16)] <-
 newdata.LR$AGE.LR.POS[which(newdata.LR$strata==16)] 	    +
   newdata.LR$TLastSurgery.LR.POS[which(newdata.LR$strata == 15)]

 newdata.LR$AGE.DR.NEG[which(newdata.LR$strata==9)] <-
 newdata.LR$AGE.LR.NEG[which(newdata.LR$strata==7)]
 newdata.LR$AGE.DR.POS[which(newdata.LR$strata==18)] <-
 newdata.LR$AGE.LR.POS[which(newdata.LR$strata==16)]

 system.time(pt[['LR']] <- getProbsLR(m, group=as.numeric(Clinical$ER.Status[i]), newdata.LR, timepoints=tt,
 beta=beta, LR=LR, DR=DR, betaAGE=betaAGE, LRAGE=LRAGE, DRAGE=DRAGE, compact=FALSE))
 tt <- c(2, 5, 10, 15, 20)
 if(Clinical$ER.Status[i]=="ER+") {
 				   beta <- c(coef(m)['TLastSurgery.LR.POS'], coef(m)['TLastSurgery.DR.POS'])
 				   betaAGE <- c(coef(m)['AGE.LR.POS'], coef(m)['AGE.DR.POS'])
 				   x <- c(newdata[15,'TLastSurgery.LR.POS'], newdata[17,'TLastSurgery.DR.POS'])
 				   xAGE <- c(x=newdata[16, 'AGE.LR.POS'] - newdata[13, 'AGE.PS.POS'],
             newdata[18, 'AGE.DR.POS'] - newdata[13, 'AGE.PS.POS'])

 }
 if(Clinical$ER.Status[i]=="ER-") {
 				   beta <- c(coef(m)['TLastSurgery.LR.NEG'], coef(m)['TLastSurgery.DR.NEG'])
 				   betaAGE <- c(coef(m)['AGE.LR.NEG'], coef(m)['AGE.DR.NEG'])
 				   x <- c(newdata[6,'TLastSurgery.LR.NEG'], newdata[8,'TLastSurgery.DR.NEG'])
 				   xAGE <- c(x=newdata[7, 'AGE.LR.NEG'] - newdata[4, 'AGE.PS.NEG'],
             newdata[9, 'AGE.DR.NEG'] - newdata[4, 'AGE.PS.NEG'])

 }
 system.time(pt[['S']] <- getProbsS(m, group=as.numeric(Clinical$ER.Status[i]), newdata, timepoints=tt,
 x=x, beta=beta, betaAGE=betaAGE, xAGE=xAGE, compact=FALSE))

 save(pt, file=paste0("./IndivProbs/ERTP_Patient", i, ".RData"))

