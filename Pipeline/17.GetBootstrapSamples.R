rm(list=ls())
setwd("./Bootstraps")
B <- 100

allB <- NULL
for (i in 1:B) {
load(paste0("BootsPredsERMODEL_", i, "_NEG.RData"))
allB <- cbind(allB, coef(mb))
}
save(allB, file="../MODELSBOOTER.RData")
allB <- NULL
for (i in 1:B) {
load(paste0("BootsPreds4GROUPSMODEL_", i, "_GROUP1.RData"))
allB <- cbind(allB, coef(mb))
}
save(allB, file="../MODELSBOOTFOURGROUPS.RData")
allB <- NULL
for (i in 1:B) {
load(paste0("BootsPredsPAM50MODEL_", i, "_GROUP1.RData"))
allB <- cbind(allB, coef(mb))
}
save(allB, file="../MODELSBOOTPAM50.RData")
allB <- NULL
for (i in 1:B) {
load(paste0("BootsPredsINTCLUSTMODEL_", i, "_GROUP1.RData"))
allB <- cbind(allB, coef(mb))
}
save(allB, file="../MODELSBOOTINTCLUST.RData")


allF <- dir(pattern=".RData")
allF <- allF[grep("ERMODEL", allF)]
ER.probs <- list()
elems <- c("Times", "p.l", "p.l.d", "p.l.c", "p.l.d.c", "p.l.o", "p.l.d.o")
for (j in c("NEG", "POS")) {
	ER.probs[[j]] <- list()
	ER.probs[[j]]$LR <- list()
	for (i in elems) {
	ER.probs[[j]]$LR[[i]] <- NULL
}
}
elems <- c("Times", "p.d", "p.d.c", "p.d.o")
for (j in c("NEG", "POS")) {
	ER.probs[[j]]$DR <- list()
	for (i in elems) {
	ER.probs[[j]]$DR[[i]] <- NULL
}
}
elems <- c("Times", "p.s", "p.s.l.c", "p.s.l.o", "p.s.l.d", "p.s.l.d.c", 
"p.s.l.d.o", "p.s.l", "p.s.d.c", "p.s.d.o", "p.s.d", "p.s.c", "p.s.o")
for (j in c("NEG", "POS")) {
	ER.probs[[j]]$S <- list()
	for (i in elems) {
	ER.probs[[j]]$S[[i]] <- NULL
}
}

library(survival)
coef.models <- list()
for (i in 1:B) {
	ff <- paste0("BootsPredsERMODEL_", i, "_NEG.RData")	
	load(ff)
	coef.models[[i]] <- coef(mb)
	ff <- paste0("BootsPredsERMODEL_", i, "_NEG.RData")	
	load(ff)
	elems <- c("Times", "p.s", "p.s.l.c", "p.s.l.o", "p.s.l.d", "p.s.l.d.c", 
	"p.s.l.d.o", "p.s.l", "p.s.d.c", "p.s.d.o", "p.s.d", "p.s.c", "p.s.o")
	for (j in elems) {
        ER.probs[["NEG"]]$S[[j]] <- cbind(ER.probs[["NEG"]]$S[[j]], pt.boot$S[[j]])
	}
	elems <- c("Times", "p.l", "p.l.d", "p.l.c", "p.l.d.c", "p.l.o", "p.l.d.o")
	for (j in elems) {
        ER.probs[["NEG"]]$LR[[j]] <- cbind(ER.probs[["NEG"]]$LR[[j]], pt.boot$LR[[j]])
	}
	elems <- c("Times", "p.d", "p.d.c", "p.d.o")
	for (j in elems) {
        ER.probs[["NEG"]]$DR[[j]] <- cbind(ER.probs[["NEG"]]$DR[[j]], pt.boot$DR[[j]])
	}
	ff <- paste0("BootsPredsERMODEL_", i, "_POS.RData")	
	load(ff)
	elems <- c("Times", "p.s", "p.s.l.c", "p.s.l.o", "p.s.l.d", "p.s.l.d.c", 
	"p.s.l.d.o", "p.s.l", "p.s.d.c", "p.s.d.o", "p.s.d", "p.s.c", "p.s.o")
	for (j in elems) {
        ER.probs[["POS"]]$S[[j]] <- cbind(ER.probs[["POS"]]$S[[j]], pt.boot$S[[j]])
	}
	elems <- c("Times", "p.l", "p.l.d", "p.l.c", "p.l.d.c", "p.l.o", "p.l.d.o")
	for (j in elems) {
        ER.probs[["POS"]]$LR[[j]] <- cbind(ER.probs[["POS"]]$LR[[j]], pt.boot$LR[[j]])
	}
	elems <- c("Times", "p.d", "p.d.c", "p.d.o")
	for (j in elems) {
        ER.probs[["POS"]]$DR[[j]] <- cbind(ER.probs[["POS"]]$DR[[j]], pt.boot$DR[[j]])
	}
cat("Bootstrap ", i, "done\n")
}
save(coef.models, ER.probs, file="BootstrapModel_ER_Results.RData")

rm(list=ls())
B <- 100
allF <- dir(pattern=".RData")
allF <- allF[grep("4GROUPSMODEL", allF)]
FOURGROUPS.probs <- list()
elems <- c("Times", "p.l", "p.l.d", "p.l.c", "p.l.d.c", "p.l.o", "p.l.d.o")
for (j in 1:4) {
	FOURGROUPS.probs[[j]] <- list()
	FOURGROUPS.probs[[j]]$LR <- list()
	for (i in elems) {
	FOURGROUPS.probs[[j]]$LR[[i]] <- NULL
}
}
elems <- c("Times", "p.d", "p.d.c", "p.d.o")
for (j in 1:4) {
	FOURGROUPS.probs[[j]]$DR <- list()
	for (i in elems) {
	FOURGROUPS.probs[[j]]$DR[[i]] <- NULL
}
}
elems <- c("Times", "p.s", "p.s.l.c", "p.s.l.o", "p.s.l.d", "p.s.l.d.c", 
"p.s.l.d.o", "p.s.l", "p.s.d.c", "p.s.d.o", "p.s.d", "p.s.c", "p.s.o")
for (j in 1:4) {
	FOURGROUPS.probs[[j]]$S <- list()
	for (i in elems) {
	FOURGROUPS.probs[[j]]$S[[i]] <- NULL
}
}
library(survival)
coef.models <- list()
for (i in 1:B) {
	ff <- paste0("BootsPreds4GROUPSMODEL_", i, "_GROUP1.RData")	
	load(ff)
	coef.models[[i]] <- coef(mb)
	for (k in 1:4) {
	ff <- paste0("BootsPreds4GROUPSMODEL_", i, "_GROUP", k, ".RData")	
	load(ff)
	elems <- c("Times", "p.s", "p.s.l.c", "p.s.l.o", "p.s.l.d", "p.s.l.d.c", 
	"p.s.l.d.o", "p.s.l", "p.s.d.c", "p.s.d.o", "p.s.d", "p.s.c", "p.s.o")
	for (j in elems) {
        FOURGROUPS.probs[[k]]$S[[j]] <- cbind(FOURGROUPS.probs[[k]]$S[[j]], pt.boot$S[[j]])
	}
	elems <- c("Times", "p.l", "p.l.d", "p.l.c", "p.l.d.c", "p.l.o", "p.l.d.o")
	for (j in elems) {
        FOURGROUPS.probs[[k]]$LR[[j]] <- cbind(FOURGROUPS.probs[[k]]$LR[[j]], pt.boot$LR[[j]])
	}
	elems <- c("Times", "p.d", "p.d.c", "p.d.o")
	for (j in elems) {
        FOURGROUPS.probs[[k]]$DR[[j]] <- cbind(FOURGROUPS.probs[[k]]$DR[[j]], pt.boot$DR[[j]])
	}
}

cat("Bootstrap ", i, "done\n")
}
save(coef.models, FOURGROUPS.probs, file="BootstrapModel_4GROUPS_Results.RData")

rm(list=ls())
B <- 100
allF <- dir(pattern=".RData")
allF <- allF[grep("PAM50MODEL", allF)]
PAM50.probs <- list()
elems <- c("Times", "p.l", "p.l.d", "p.l.c", "p.l.d.c", "p.l.o", "p.l.d.o")
for (j in 1:5) {
	PAM50.probs[[j]] <- list()
	PAM50.probs[[j]]$LR <- list()
	for (i in elems) {
	PAM50.probs[[j]]$LR[[i]] <- NULL
}
}
elems <- c("Times", "p.d", "p.d.c", "p.d.o")
for (j in 1:5) {
	PAM50.probs[[j]]$DR <- list()
	for (i in elems) {
	PAM50.probs[[j]]$DR[[i]] <- NULL
}
}
elems <- c("Times", "p.s", "p.s.l.c", "p.s.l.o", "p.s.l.d", "p.s.l.d.c", 
"p.s.l.d.o", "p.s.l", "p.s.d.c", "p.s.d.o", "p.s.d", "p.s.c", "p.s.o")
for (j in 1:5) {
	PAM50.probs[[j]]$S <- list()
	for (i in elems) {
	PAM50.probs[[j]]$S[[i]] <- NULL
}
}

library(survival)
coef.models <- list()
for (i in c(1:23, 25:97, 99:B)) {
	ff <- paste0("BootsPredsPAM50MODEL_", i, "_GROUP1.RData")	
	load(ff)
	coef.models[[i]] <- coef(mb)
	for (k in 1:5) {
	ff <- paste0("BootsPredsPAM50MODEL_", i, "_GROUP", k, ".RData")	
	load(ff)
	elems <- c("Times", "p.s", "p.s.l.c", "p.s.l.o", "p.s.l.d", "p.s.l.d.c", 
	"p.s.l.d.o", "p.s.l", "p.s.d.c", "p.s.d.o", "p.s.d", "p.s.c", "p.s.o")
	for (j in elems) {
        PAM50.probs[[k]]$S[[j]] <- cbind(PAM50.probs[[k]]$S[[j]], pt.boot$S[[j]])
	}
	elems <- c("Times", "p.l", "p.l.d", "p.l.c", "p.l.d.c", "p.l.o", "p.l.d.o")
	for (j in elems) {
        PAM50.probs[[k]]$LR[[j]] <- cbind(PAM50.probs[[k]]$LR[[j]], pt.boot$LR[[j]])
	}
	elems <- c("Times", "p.d", "p.d.c", "p.d.o")
	for (j in elems) {
        PAM50.probs[[k]]$DR[[j]] <- cbind(PAM50.probs[[k]]$DR[[j]], pt.boot$DR[[j]])
	}

}
cat("Bootstrap ", i, "done\n")
}
save(coef.models, PAM50.probs, file="BootstrapModel_PAM50_Results.RData")

rm(list=ls())
B <- 100
allF <- dir(pattern=".RData")
allF <- allF[grep("INTCLUSTMODEL", allF)]
INTCLUST.probs <- list()
elems <- c("Times", "p.l", "p.l.d", "p.l.c", "p.l.d.c", "p.l.o", "p.l.d.o")
for (j in 1:11) {
	INTCLUST.probs[[j]] <- list()
	INTCLUST.probs[[j]]$LR <- list()
	for (i in elems) {
	INTCLUST.probs[[j]]$LR[[i]] <- NULL
}
}
elems <- c("Times", "p.d", "p.d.c", "p.d.o")
for (j in 1:11) {
	INTCLUST.probs[[j]]$DR <- list()
	for (i in elems) {
	INTCLUST.probs[[j]]$DR[[i]] <- NULL
}
}
elems <- c("Times", "p.s", "p.s.l.c", "p.s.l.o", "p.s.l.d", "p.s.l.d.c", 
"p.s.l.d.o", "p.s.l", "p.s.d.c", "p.s.d.o", "p.s.d", "p.s.c", "p.s.o")
for (j in 1:11) {
	INTCLUST.probs[[j]]$S <- list()
	for (i in elems) {
	INTCLUST.probs[[j]]$S[[i]] <- NULL
}
}

library(survival)
coef.models <- list()
for (i in 1:B) {
	ff <- paste0("BootsPredsINTCLUSTMODEL_", i, "_GROUP1.RData")	
	load(ff)
	coef.models[[i]] <- coef(mb)
	for (k in 1:11) {
	ff <- paste0("BootsPredsINTCLUSTMODEL_", i, "_GROUP", k, ".RData")	
	load(ff)
	elems <- c("Times", "p.s", "p.s.l.c", "p.s.l.o", "p.s.l.d", "p.s.l.d.c", 
	"p.s.l.d.o", "p.s.l", "p.s.d.c", "p.s.d.o", "p.s.d", "p.s.c", "p.s.o")
	for (j in elems) {
        INTCLUST.probs[[k]]$S[[j]] <- cbind(INTCLUST.probs[[k]]$S[[j]], pt.boot$S[[j]])
	}
	elems <- c("Times", "p.l", "p.l.d", "p.l.c", "p.l.d.c", "p.l.o", "p.l.d.o")
	for (j in elems) {
        INTCLUST.probs[[k]]$LR[[j]] <- cbind(INTCLUST.probs[[k]]$LR[[j]], pt.boot$LR[[j]])
	}
	elems <- c("Times", "p.d", "p.d.c", "p.d.o")
	for (j in elems) {
        INTCLUST.probs[[k]]$DR[[j]] <- cbind(INTCLUST.probs[[k]]$DR[[j]], pt.boot$DR[[j]])
	}

}
cat("Bootstrap ", i, "done\n")
}
save(coef.models, INTCLUST.probs, file="BootstrapModel_INTCLUST_Results.RData")