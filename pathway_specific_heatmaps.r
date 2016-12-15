library(KEGGREST)

getMetabolites <- function(pathway) {
    names(keggGet(pathway)[[1]]$COMPOUND)
}

getGenes <- function(pathway) {
    temp <- keggLink('genes', pathway)
    gsub('(.+?\\:)', '', temp)
}

drawMetsHeatmapForPathway <- function(exprs, kegg.ids, groups, pathway) {
    temp <- merge(exprs, kegg.ids, by = 0)
    temp <- subset(temp, temp$KEGG %in% getMetabolites(pathway))

    rownames(temp) <- temp$BIOCHEMICAL
    temp <- temp[, 2:15]

    drawHeatMap(temp, groups)
}

drawGeneHeatmapForPathway <- function(exprs, groups, pathway = NULL, 
                                      genes = NULL) {
    genes <- if (is.null(genes)) getGenes(pathway) else genes
    temp <- subset(exprs,  exprs$EntrezGene %in% genes)

    rownames(temp) <- temp$Symbol
    temp <- temp[, -c(1, 2)]

    drawHeatMap(temp, groups)
}
