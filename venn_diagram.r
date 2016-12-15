setwd("~/Documents/brittany/myc_induction_regression_transcriptomics_062514")

tx.pathway.file <- "pathway_analysis_ORA_kegg.csv"
mx.pathway.file <- "mx_pathway_analysis_ORA_kegg.csv"

tx.pathways <- read.csv(tx.pathway.file, as.is = T, fill = T, check.names = F)
mx.pathways <- read.csv(mx.pathway.file, as.is = T, fill = T, check.names = F)

## filter TX pathways such that FDR < 0.05
tx.pathways <- subset(tx.pathways, tx.pathways$p.value < 0.05)

## filter MX pathways such that p.value < 0.05
mx.pathways <- subset(mx.pathways, mx.pathways$p.value < 0.05)

pathways <- unique(c(tx.pathways$Pathway, mx.pathways$Pathway))
venn <- sapply(pathways, function(p) {
    result = ""

    if (p %in% tx.pathways$Pathway)
        result <- paste(result, "T", sep = "")

    if (p %in% mx.pathways$Pathway)
        result <- paste(result, "M", sep = "")

    result
})

out <- data.frame(Pathway = names(venn), Occurrence = venn, row.names = NULL)
write.csv(out, file = "mx_tx_pathway_analysis_common_hits.csv", row.names = F)
