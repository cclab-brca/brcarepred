rm(list=ls())
## ER


S <- c(2, 5, 10, 15, 20)
TP <- list()
for (i in c("NEG", "POS")) {
    load(paste0("../Models/TP/ERTP_", i, ".RData"))
    times <- pt$S$Times
    ids.S <- sapply(S, function(x) which.min(abs(times-x)))
    times <- pt$LR$Times
    ids.LR <- sapply(S, function(x) which.min(abs(times-x)))
    times <- pt$DR$Times
    ids.DR <- sapply(S, function(x) which.min(abs(times-x)))
    TP[[i]] <- list(pt$S$p.s[ids.S], pt$S$p.s.l[ids.S], I(pt$S$p.s.l.d +
                        pt$S$p.s.d)[ids.S], I(pt$S$p.s.l.c +
                            pt$S$p.s.l.d.c + pt$S$p.s.d.c +
                                pt$S$p.s.c)[ids.S],
                    I(pt$S$p.s.l.o + pt$S$p.s.l.d.o +
                          pt$S$p.s.d.o + pt$S$p.s.o)[ids.S],
                    pt$LR$p.l[ids.LR], pt$LR$p.l.d[ids.LR],
                    I(pt$LR$p.l.c + pt$LR$p.l.d.c)[ids.LR],
                 I(pt$LR$p.l.o + pt$LR$p.l.d.o)[ids.LR],
                    pt$DR$p.d[ids.DR], pt$DR$p.d.c[ids.DR],
                 pt$DR$p.d.o[ids.DR])
}
pt <-  TP
rm(TP)

library(survival)
load(file="../Models/Bootstraps/BootstrapModel_ER_Results.RData")


boot.TP <- list()
for (i in 1:2) {
    boot.TP[[i]] <- list(ER.probs[[i]]$S$p.s[ids.S,], ER.probs[[i]]$S$p.s.l[ids.S,], I(ER.probs[[i]]$S$p.s.l.d +
                        ER.probs[[i]]$S$p.s.d)[ids.S,], I(ER.probs[[i]]$S$p.s.l.c +
                            ER.probs[[i]]$S$p.s.l.d.c + ER.probs[[i]]$S$p.s.d.c +
                                ER.probs[[i]]$S$p.s.c)[ids.S,],
                    I(ER.probs[[i]]$S$p.s.l.o + ER.probs[[i]]$S$p.s.l.d.o +
                          ER.probs[[i]]$S$p.s.d.o + ER.probs[[i]]$S$p.s.o)[ids.S,],
                    ER.probs[[i]]$LR$p.l[ids.LR,], ER.probs[[i]]$LR$p.l.d[ids.LR,],
                    I(ER.probs[[i]]$LR$p.l.c + ER.probs[[i]]$LR$p.l.d.c)[ids.LR,],
                 I(ER.probs[[i]]$LR$p.l.o + ER.probs[[i]]$LR$p.l.d.o)[ids.LR,],
                    ER.probs[[i]]$DR$p.d[ids.DR,], ER.probs[[i]]$DR$p.d.c[ids.DR,],
                 ER.probs[[i]]$DR$p.d.o[ids.DR,])
}

ER.mean <- lapply(boot.TP, function(x) lapply(x, function(y) apply(y, 1, mean)))
ER.sd <- lapply(boot.TP, function(x) lapply(x, function(y) apply(y, 1, sd)))

names(ER.sd) <- names(pt)
labs <- names(pt)
all.names <- c("Surgery-Surgery", "Surgery-LocoRegional Relapse", "Surgery-Distant Relapse", "Surgery-Cancer Death",
               "Surgery-Other Death", "LocoRegional Relapse-LocoRegional Relapse",
               "LocoRegional Relapse-Distant Relapse", "LocoRegional Relapse-Cancer Death",
               "LocoRegional Relapse-Other Death", "Distant Relapse-Distant Relapse",
               "Distant Relapse-Cancer Death", "Distant Relapse-Other Death")
for (j in labs) {
    sink(paste("TableS3_", sub("/", "_", j, fixed=T), ".txt", sep=""))
    cat("\t2y\t5y\t10y\t15y\t20y\t\t2y\t5y\t10y\t15y\t20y\n")
    for (counter in 1:length(all.names)) {
        cat(all.names[counter], "\t")
        for (i in 1:5) {
            cat(round(pt[[j]][[counter]][i], 3), "\t")
        }
        cat("\t")
        for (i in 1:5) {
            cat(round(ER.sd[[j]][[counter]][i], 3), "\t")
        }
        cat("\n")
    }
    sink()
}

rm(list=ls())
## Clinical groups

S <- c(2, 5, 10, 15, 20)
TP <- list()
for (i in 1:4) {
    load(paste0("../Models/TP/ERTP_4GROUPS_", i, ".RData"))
    times <- pt$S$Times
    ids.S <- sapply(S, function(x) which.min(abs(times-x)))
    times <- pt$LR$Times
    ids.LR <- sapply(S, function(x) which.min(abs(times-x)))
    times <- pt$DR$Times
    ids.DR <- sapply(S, function(x) which.min(abs(times-x)))
    TP[[i]] <- list(pt$S$p.s[ids.S], pt$S$p.s.l[ids.S], I(pt$S$p.s.l.d +
                        pt$S$p.s.d)[ids.S], I(pt$S$p.s.l.c +
                            pt$S$p.s.l.d.c + pt$S$p.s.d.c +
                                pt$S$p.s.c)[ids.S],
                    I(pt$S$p.s.l.o + pt$S$p.s.l.d.o +
                          pt$S$p.s.d.o + pt$S$p.s.o)[ids.S],
                    pt$LR$p.l[ids.LR], pt$LR$p.l.d[ids.LR],
                    I(pt$LR$p.l.c + pt$LR$p.l.d.c)[ids.LR],
                 I(pt$LR$p.l.o + pt$LR$p.l.d.o)[ids.LR],
                    pt$DR$p.d[ids.DR], pt$DR$p.d.c[ids.DR],
                 pt$DR$p.d.o[ids.DR])
}
pt <-  TP
rm(TP)

library(survival)
load(file="../Models/Bootstraps/BootstrapModel_4Groups_Results.RData")
boot.TP <- list()
for (i in 1:4) {
    boot.TP[[i]] <- list(FOURGROUPS.probs[[i]]$S$p.s[ids.S,], FOURGROUPS.probs[[i]]$S$p.s.l[ids.S,], I(FOURGROUPS.probs[[i]]$S$p.s.l.d +
                        FOURGROUPS.probs[[i]]$S$p.s.d)[ids.S,], I(FOURGROUPS.probs[[i]]$S$p.s.l.c +
                            FOURGROUPS.probs[[i]]$S$p.s.l.d.c + FOURGROUPS.probs[[i]]$S$p.s.d.c +
                                FOURGROUPS.probs[[i]]$S$p.s.c)[ids.S,],
                    I(FOURGROUPS.probs[[i]]$S$p.s.l.o + FOURGROUPS.probs[[i]]$S$p.s.l.d.o +
                          FOURGROUPS.probs[[i]]$S$p.s.d.o + FOURGROUPS.probs[[i]]$S$p.s.o)[ids.S,],
                    FOURGROUPS.probs[[i]]$LR$p.l[ids.LR,], FOURGROUPS.probs[[i]]$LR$p.l.d[ids.LR,],
                    I(FOURGROUPS.probs[[i]]$LR$p.l.c + FOURGROUPS.probs[[i]]$LR$p.l.d.c)[ids.LR,],
                 I(FOURGROUPS.probs[[i]]$LR$p.l.o + FOURGROUPS.probs[[i]]$LR$p.l.d.o)[ids.LR,],
                    FOURGROUPS.probs[[i]]$DR$p.d[ids.DR,], FOURGROUPS.probs[[i]]$DR$p.d.c[ids.DR,],
                 FOURGROUPS.probs[[i]]$DR$p.d.o[ids.DR,])
}

FOURGROUPS.mean <- lapply(boot.TP, function(x) lapply(x, function(y) apply(y, 1, mean)))
FOURGROUPS.sd <- lapply(boot.TP, function(x) lapply(x, function(y) apply(y, 1, sd)))

labs <- c("ER-/HER2-", "ER-/HER2+", "ER+/HER2-", "ER+/HER2+")
names(pt) <- labs


labs <- names(pt)
names(FOURGROUPS.sd) <- labs
all.names <- c("Surgery-Surgery", "Surgery-LocoRegional Relapse", "Surgery-Distant Relapse", "Surgery-Cancer Death",
               "Surgery-Other Death", "LocoRegional Relapse-LocoRegional Relapse",
               "LocoRegional Relapse-Distant Relapse", "LocoRegional Relapse-Cancer Death",
               "LocoRegional Relapse-Other Death", "Distant Relapse-Distant Relapse",
               "Distant Relapse-Cancer Death", "Distant Relapse-Other Death")
for (j in labs) {
    sink(paste("TableS3_", sub("/", "_", j, fixed=T), ".txt", sep=""))
    cat("\t2y\t5y\t10y\t15y\t20y\t\t2y\t5y\t10y\t15y\t20y\n")
    for (counter in 1:length(all.names)) {
        cat(all.names[counter], "\t")
        for (i in 1:5) {
            cat(round(pt[[j]][[counter]][i], 3), "\t")
        }
        cat("\t")
        for (i in 1:5) {
            cat(round(FOURGROUPS.sd[[j]][[counter]][i], 3), "\t")
        }
        cat("\n")
    }
    sink()
}

rm(list=ls())
## Pam50
S <- c(2, 5, 10, 15, 20)
TP <- list()
for (i in 1:5) {
    load(paste0("../Models/TP/TP_Pam50_", i, ".RData"))
    times <- pt$S$Times
    ids.S <- sapply(S, function(x) which.min(abs(times-x)))
    times <- pt$LR$Times
    ids.LR <- sapply(S, function(x) which.min(abs(times-x)))
    times <- pt$DR$Times
    ids.DR <- sapply(S, function(x) which.min(abs(times-x)))
    TP[[i]] <- list(pt$S$p.s[ids.S], pt$S$p.s.l[ids.S], I(pt$S$p.s.l.d +
                        pt$S$p.s.d)[ids.S], I(pt$S$p.s.l.c +
                            pt$S$p.s.l.d.c + pt$S$p.s.d.c +
                                pt$S$p.s.c)[ids.S],
                    I(pt$S$p.s.l.o + pt$S$p.s.l.d.o +
                          pt$S$p.s.d.o + pt$S$p.s.o)[ids.S],
                    pt$LR$p.l[ids.LR], pt$LR$p.l.d[ids.LR],
                    I(pt$LR$p.l.c + pt$LR$p.l.d.c)[ids.LR],
                 I(pt$LR$p.l.o + pt$LR$p.l.d.o)[ids.LR],
                    pt$DR$p.d[ids.DR], pt$DR$p.d.c[ids.DR],
                 pt$DR$p.d.o[ids.DR])
}
pt <-  TP
rm(TP)

library(survival)
load(file="../Models/Bootstraps/BootstrapModel_Pam50_Results.RData")

boot.TP <- list()
for (i in 1:5) {
    boot.TP[[i]] <- list(PAM50.probs[[i]]$S$p.s[ids.S,], PAM50.probs[[i]]$S$p.s.l[ids.S,], I(PAM50.probs[[i]]$S$p.s.l.d +
                        PAM50.probs[[i]]$S$p.s.d)[ids.S,], I(PAM50.probs[[i]]$S$p.s.l.c +
                            PAM50.probs[[i]]$S$p.s.l.d.c + PAM50.probs[[i]]$S$p.s.d.c +
                                PAM50.probs[[i]]$S$p.s.c)[ids.S,],
                    I(PAM50.probs[[i]]$S$p.s.l.o + PAM50.probs[[i]]$S$p.s.l.d.o +
                          PAM50.probs[[i]]$S$p.s.d.o + PAM50.probs[[i]]$S$p.s.o)[ids.S,],
                    PAM50.probs[[i]]$LR$p.l[ids.LR,], PAM50.probs[[i]]$LR$p.l.d[ids.LR,],
                    I(PAM50.probs[[i]]$LR$p.l.c + PAM50.probs[[i]]$LR$p.l.d.c)[ids.LR,],
                 I(PAM50.probs[[i]]$LR$p.l.o + PAM50.probs[[i]]$LR$p.l.d.o)[ids.LR,],
                    PAM50.probs[[i]]$DR$p.d[ids.DR,], PAM50.probs[[i]]$DR$p.d.c[ids.DR,],
                 PAM50.probs[[i]]$DR$p.d.o[ids.DR,])
}


PAM50.mean <- lapply(boot.TP, function(x) lapply(x, function(y) apply(y, 1, mean)))
PAM50.sd <- lapply(boot.TP, function(x) lapply(x, function(y) apply(y, 1, sd)))
labs <- c("Basal", "Her2", "LumA", "LumB", "Normal")
names(pt) <- labs


## Table

names(PAM50.sd) <- labs
all.names <- c("Surgery-Surgery", "Surgery-LocoRegional Relapse", "Surgery-Distant Relapse", "Surgery-Cancer Death",
               "Surgery-Other Death", "LocoRegional Relapse-LocoRegional Relapse",
               "LocoRegional Relapse-Distant Relapse", "LocoRegional Relapse-Cancer Death",
               "LocoRegional Relapse-Other Death", "Distant Relapse-Distant Relapse",
               "Distant Relapse-Cancer Death", "Distant Relapse-Other Death")
for (j in labs) {
    sink(paste("TableS3_", sub("/", "_", j, fixed=T), ".txt", sep=""))
    cat("\t2y\t5y\t10y\t15y\t20y\t\t2y\t5y\t10y\t15y\t20y\n")
    for (counter in 1:length(all.names)) {
        cat(all.names[counter], "\t")
        for (i in 1:5) {
            cat(round(pt[[j]][[counter]][i], 3), "\t")
        }
        cat("\t")
        for (i in 1:5) {
            cat(round(PAM50.sd[[j]][[counter]][i], 3), "\t")
        }
        cat("\n")
    }
    sink()
}


rm(list=ls())
## IntClust

TP <- list()
for (i in 1:11) {
    load(paste0("../Models/TP/TP_IntClust_", i, ".RData"))
    times <- pt$S$Times
    ids.S <- sapply(S, function(x) which.min(abs(times-x)))
    times <- pt$LR$Times
    ids.LR <- sapply(S, function(x) which.min(abs(times-x)))
    times <- pt$DR$Times
    ids.DR <- sapply(S, function(x) which.min(abs(times-x)))
    TP[[i]] <- list(pt$S$p.s[ids.S], pt$S$p.s.l[ids.S], I(pt$S$p.s.l.d +
                        pt$S$p.s.d)[ids.S], I(pt$S$p.s.l.c +
                            pt$S$p.s.l.d.c + pt$S$p.s.d.c +
                                pt$S$p.s.c)[ids.S],
                    I(pt$S$p.s.l.o + pt$S$p.s.l.d.o +
                          pt$S$p.s.d.o + pt$S$p.s.o)[ids.S],
                    pt$LR$p.l[ids.LR], pt$LR$p.l.d[ids.LR],
                    I(pt$LR$p.l.c + pt$LR$p.l.d.c)[ids.LR],
                 I(pt$LR$p.l.o + pt$LR$p.l.d.o)[ids.LR],
                    pt$DR$p.d[ids.DR], pt$DR$p.d.c[ids.DR],
                 pt$DR$p.d.o[ids.DR])
}
pt <-  TP
rm(TP)


load(file="../Models/Bootstraps/BootstrapModel_INTCLUST_Results.RData")


boot.TP <- list()
for (i in 1:11) {
    boot.TP[[i]] <- list(INTCLUST.probs[[i]]$S$p.s[ids.S,], INTCLUST.probs[[i]]$S$p.s.l[ids.S,], I(INTCLUST.probs[[i]]$S$p.s.l.d +
                        INTCLUST.probs[[i]]$S$p.s.d)[ids.S,], I(INTCLUST.probs[[i]]$S$p.s.l.c +
                            INTCLUST.probs[[i]]$S$p.s.l.d.c + INTCLUST.probs[[i]]$S$p.s.d.c +
                                INTCLUST.probs[[i]]$S$p.s.c)[ids.S,],
                    I(INTCLUST.probs[[i]]$S$p.s.l.o + INTCLUST.probs[[i]]$S$p.s.l.d.o +
                          INTCLUST.probs[[i]]$S$p.s.d.o + INTCLUST.probs[[i]]$S$p.s.o)[ids.S,],
                    INTCLUST.probs[[i]]$LR$p.l[ids.LR,], INTCLUST.probs[[i]]$LR$p.l.d[ids.LR,],
                    I(INTCLUST.probs[[i]]$LR$p.l.c + INTCLUST.probs[[i]]$LR$p.l.d.c)[ids.LR,],
                 I(INTCLUST.probs[[i]]$LR$p.l.o + INTCLUST.probs[[i]]$LR$p.l.d.o)[ids.LR,],
                    INTCLUST.probs[[i]]$DR$p.d[ids.DR,], INTCLUST.probs[[i]]$DR$p.d.c[ids.DR,],
                 INTCLUST.probs[[i]]$DR$p.d.o[ids.DR,])
}



BOOT.mean <- lapply(boot.TP, function(x) lapply(x, function(y) apply(y, 1, mean)))
BOOT.sd <- lapply(boot.TP, function(x) lapply(x, function(y) apply(y, 1, sd)))
labs <- c(1:3, "4ER+", "4ER-", 5:10)
names(pt) <- labs


## Table

names(BOOT.sd) <- labs
all.names <- c("Surgery-Surgery", "Surgery-LocoRegional Relapse", "Surgery-Distant Relapse", "Surgery-Cancer Death",
               "Surgery-Other Death", "LocoRegional Relapse-LocoRegional Relapse",
               "LocoRegional Relapse-Distant Relapse", "LocoRegional Relapse-Cancer Death",
               "LocoRegional Relapse-Other Death", "Distant Relapse-Distant Relapse",
               "Distant Relapse-Cancer Death", "Distant Relapse-Other Death")
for (j in labs) {
    sink(paste("TableS3_", sub("/", "_", j, fixed=T), ".txt", sep=""))
    cat("\t2y\t5y\t10y\t15y\t20y\t\t2y\t5y\t10y\t15y\t20y\n")
    for (counter in 1:length(all.names)) {
        cat(all.names[counter], "\t")
        for (i in 1:5) {
            cat(round(pt[[j]][[counter]][i], 3), "\t")
        }
        cat("\t")
        for (i in 1:5) {
            cat(round(BOOT.sd[[j]][[counter]][i], 3), "\t")
        }
        cat("\n")
    }
    sink()
}

