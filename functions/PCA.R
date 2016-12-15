library(FactoMineR)

PCAAnalysis <- function(data, groups, conf.interval = 0.95) {
    pca.exprs <- cbind(groups, as.data.frame(t(data)))
    pca.res <- PCA(pca.exprs, quali.sup = 1, ncp = 5, graph = F)

    concat <- cbind.data.frame(pca.exprs[, 1], pca.res$ind$coord)
    ellipse.coord <- coord.ellipse(concat, level.conf = conf.interval, 
        bary = T)

    plot.PCA(pca.res, axes=c(1, 2), choix = "ind", habillage = 1, 
        ellipse = ellipse.coord, label = "quali")
}
