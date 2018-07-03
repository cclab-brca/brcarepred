library(Rgraphviz)
library(RColorBrewer)
pdf("Figure1.pdf",
    width=5, height=8)
cols <- c("black", "lightblue", "darkblue", "red", "olivedrab")

Nodes <- c("Surgery","LocoRegionalRelapse", "DistantRelapse",
                            "Death/Cancer", "Death/Other")
Edges <- list("Surgery"=list(edges=c("Death/Other", "Death/Cancer", "LocoRegionalRelapse", "DistantRelapse")),
              "LocoRegionalRelapse"=list(edges=c("DistantRelapse", "Death/Other", "Death/Cancer")),
              "DistantRelapse"=list(edges=c("Death/Other", "Death/Cancer")),
              "Death/Other"=list(edges=NULL),
              "Death/Cancer"=list(edges=NULL))
treeobj <- new("graphNEL", nodes=Nodes, edgeL=Edges, edgemode="directed")
eAttrs<-list()
eAttrs$color <- c("Surgery~LocoRegionalRelapse" = cols[2],
                  "Surgery~DistantRelapse" = cols[3],
                  "Surgery~Death/Cancer" = cols[4],
                  "Surgery~Death/Other" = cols[5],
                  "LocoRegionalRelapse~DistantRelapse" = cols[3],
                  "LocoRegionalRelapse~Death/Cancer" = cols[4],
                  "LocoRegionalRelapse~Death/Other" = cols[5],
                  "DistantRelapse~Death/Cancer" = cols[4],
                  "DistantRelapse~Death/Other" = cols[5])
eAttrs$lwd <- c("Surgery~LocoRegionalRelapse" = 3,
                  "Surgery~DistantRelapse" = 3,
                  "Surgery~Death/Cancer" = 3,
                  "Surgery~Death/Other" = 3,
                  "LocoRegionalRelapse~DistantRelapse" = 3,
                  "LocoRegionalRelapse~Death/Cancer" = 3,
                  "LocoRegionalRelapse~Death/Other" = 3,
                  "DistantRelapse~Death/Cancer" = 3,
                  "DistantRelapse~Death/Other" = 3)
nAttrs<-list()
nAttrs$label <- c("Surgery"=expression("Surgery"),
                  "LocoRegionalRelapse"="LocoRegional\\\nRelapse",
                  "DistantRelapse"="Distant\\\nRelapse",
                  "Death/Cancer"="Death\\\nCancer",
                  "Death/Other"="Death\\\nOther")
nAttrs$fontsize <- c("Surgery"=15,
                  "LocoRegionalRelapse"=15,
                  "DistantRelapse"=15,
                  "Death/Cancer"=15,
                  "Death/Other"=15)
eAttrs$label <- c("Surgery~LocoRegionalRelapse" =
                      expression(lambda(t,G, S, LN)),
                  "Surgery~DistantRelapse" =
                     expression(lambda(t,G, S, LN)),
                  "Surgery~Death/Cancer" =
                     expression(lambda(t,G, S, LN)),
                  "Surgery~Death/Other" =
                     expression(lambda(t,A)),
                  "LocoRegionalRelapse~DistantRelapse" =
                     expression(lambda(t,G, S, LN, TD)),
                  "LocoRegionalRelapse~Death/Cancer" =
                     expression(lambda(t,G, S, LN, TD)),
                  "LocoRegionalRelapse~Death/Other" =
                     expression(lambda(t, A)),
                  "DistantRelapse~Death/Cancer" =
                     expression(lambda(t,G, S, LN, TD)),
                  "DistantRelapse~Death/Other" =
                      expression(lambda(t,A)))
plot(treeobj, cex=2, edgeAttrs = eAttrs, nodeAttrs=nAttrs)
dev.off()

