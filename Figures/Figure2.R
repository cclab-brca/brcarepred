library(mstate)
library(lattice)
load(file="../Models/ERM.RData")
fm <- fm2

model.cofs <- data.frame(Name=names(coef(fm)), HR=summary(fm)$conf.int[,1], LI=summary(fm)$conf.int[,3], UI=summary(fm)$conf.int[,4])
model.cofs$Effect <- sapply(strsplit(as.character(model.cofs$Name), ".", fixed=TRUE), function(x) x[1])
model.cofs$Time <- sapply(strsplit(as.character(model.cofs$Name), ".", fixed=TRUE), function(x) x[2])
model.cofs$ER <- sapply(strsplit(as.character(model.cofs$Name), ".", fixed=TRUE), function(x) x[3])
model.cofs$Time <- factor(model.cofs$Time, levels=c("PS", "LR", "DR"))
model.cofs$ER <- factor(model.cofs$ER, levels=c("NEG", "POS"), labels=c("ER-", "ER+"))
model.cofs <- model.cofs[which(!is.na(model.cofs$ER)),]
model.cofs$Effect[which(model.cofs$Effect=="TLastSurgery")] <- "Time from Surgery"
model.cofs$Effect[which(model.cofs$Effect=="TLastLocal")] <- "Time from LR"

## First plot: effect on time of each clinical parameter?
prepanel.ci <- function(x, y, lx, ux, subscripts, ...)
{
    x <- as.numeric(y)
    lx <- as.numeric(lx[subscripts])
    ux <- as.numeric(ux[subscripts])
    print(range(x, ux, lx, finite = TRUE))
    list(ylim = range(x, ux, lx, finite = TRUE))
}


panel.ci <- function(x, y, lx, ux, subscripts, pch = 16, ...)
{
    x <- as.numeric(x)
    y <- as.numeric(y)
    lx <- as.numeric(lx[subscripts])
    ux <- as.numeric(ux[subscripts])
    panel.dotplot(x, y, pch = pch, col.line="white", ...)
    panel.abline(h=0, col="grey")
    panel.arrows(x, lx, x, ux, col = 'black',
                 length = 0.25, unit = "native",
                 angle = 90, code = 3, lwd=0.5)

}
model.cofs$Effect <- factor(model.cofs$Effect, levels=c("GRADE", "LN", "SIZE", "Time from Surgery", "Time from LR"))
pdf("Figure2.pdf", width=9, height=7)
print(with(model.cofs,
     dotplot(log(HR) ~ Time | Effect*ER,
            lx = log(LI), ux = log(UI),
            prepanel = prepanel.ci,
             panel = panel.ci, ylab="Log-hazard Ratio")))
dev.off()

