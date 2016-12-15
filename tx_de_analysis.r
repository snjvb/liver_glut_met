setwd("~/Documents/brittany/myc_induction_regression_transcriptomics_062514")

library(cluster)
library(combinat)
library(FactoMineR)
library(ggplot2)
library(gplots)
library(grid)
library(KEGGREST)
library(limma)
library(plyr)
library(RColorBrewer)

source("~/Documents/brittany/myc_induction_regression_transcriptomics_062514/functions/KEGGLRpath.R")
source("~/Documents/brittany/myc_induction_regression_transcriptomics_062514/functions/ReactomeLRpath.R")
source("~/Documents/brittany/myc_induction_regression_transcriptomics_062514/functions/FishersExactTest.R")
source("~/Documents/brittany/myc_induction_regression_transcriptomics_062514/functions/PCA.R")
source("~/Documents/brittany/myc_induction_regression_transcriptomics_062514/gene_set.r")
source("~/Documents/brittany/myc_induction_regression_transcriptomics_062514/fig_3e_gene_set.r")

exprs.file <- "~/Documents/brittany/myc_induction_regression_transcriptomics_062514/raw_data/tx_exprs.txt"

## functions
drawHeatMap <- function(data, groups) {
    d <- data.matrix(data)

    ## cluster rows and samples by pearson distance and ward linkage
    hc <- as.dendrogram(agnes(as.dist(1 - cor(d, method = "pearson")),
        method = "ward"))
    hr <- as.dendrogram(agnes(as.dist(1 - cor(t(d), method = "pearson")),
        method = "ward"))

    unique.groups <- unique(groups)
    unique.cols <- brewer.pal(8, "Accent")[1:length(unique.groups)]
    col.cols <- as.character(mapvalues(groups, from = unique.groups,
        to = unique.cols))

    heatmap.2(data.matrix(d), margins = c(5, 25),
        Colv = hc, Rowv = hr, scale = "row",
        trace = "none", col = bluered(50), keysize = 1,
        ColSideColors = col.cols, labCol = groups,
        density.info = 'none',
        cexRow = 0.8 + 1/log10(nrow(d)))
}

drawVolcanoPlot <- function(data) {
    temp <- data[, c("logFC", "rawP")]
    temp$legend <- factor(apply(temp, 1, function(row) {
        log.fc <- row[["logFC"]]
        p <- row[["rawP"]]

        if (p < 0.05) {
            if (abs(log.fc) >= 1) {
                return("A")
            } else {
                return("B")
            }
        } else {
            return("C")
        }
    }))

    ## legend labels
    legend.counts <- table(temp$legend)
    legend.labels <- c(paste("Log2 Fold Change >= 1\nRaw p-value < 0.05\n(N=",
                             legend.counts["A"], ")", sep = ""),
                       paste("Log2 Fold Change < 1\nRaw p-value < 0.05\n(N=",
                             legend.counts["B"], ")", sep = ""),
                       paste("Raw p-value >= 0.05\n(N=", legend.counts["C"],
                             ")", sep = ""))

    ggplot(data = temp, aes(x = logFC, y = -log10(rawP),
        colour = legend)) +
        geom_point(alpha = 0.7, size = 4) +
        scale_colour_brewer(name = "Legend", breaks = c("A", "B", "C"),
            labels = legend.labels, palette = "Set1") +
        geom_vline(xintercept = c(-1, 1), linetype = "longdash") +
        geom_hline(yintercept = -log10(0.05), linetype = "longdash") +
        xlab("Log2 Fold Change") + ylab("-log10(Raw p-value)") +
        theme_bw() +
        theme(axis.title.x = element_text(size = 14, face = "bold"),
            axis.text.x = element_text(size = 12),
            axis.title.y = element_text(size = 14, face = "bold"),
            axis.text.y = element_text(size = 12),
            legend.title = element_text(size = 16),
            legend.text = element_text(size = 12),
            legend.key.size = unit(1.5, "cm"))

}

## read in normalized expression values
norm <- read.delim(exprs.file, as.is = T, fill = T, check.names = F,
    row.names = 1)
exprs <- data.matrix(norm[, 11:25])
cond <- gsub("(.*)\\.[0-9]*", "\\1", colnames(exprs))

## collapse probes onto entrez id
norm <- norm[order(rowMeans(exprs), decreasing = T), ]
exprs2 <- norm[!duplicated(norm$EntrezGene), ]

## differenfial expression analysis using limma
design <- model.matrix(~0 + cond)
colnames(design) <- gsub("cond", "", colnames(design))

fit <- lmFit(exprs2[, 11:25], design = design)
contrasts <- makeContrasts(ULSTvsCTRL = ULST - CTRL, TLSTvsULST = TLST - ULST,
    levels = design)
fit2 <- contrasts.fit(fit, contrasts = contrasts)
fit2 <- eBayes(fit2)

res <- data.frame(
  ULSTvsCTRL.logFC=fit2$coefficients[, "ULSTvsCTRL"],
  ULSTvsCTRL.raw.p=fit2$p.value[, "ULSTvsCTRL"],
  ULSTvsCTRL.fdr.p=p.adjust(fit2$p.value[, "ULSTvsCTRL"], "fdr"),
  ULSTvsCTRL.holm.p=p.adjust(fit2$p.value[, "ULSTvsCTRL"], "holm"),
  TLSTvsULST.logFC=fit2$coefficients[, "TLSTvsULST"],
  TLSTvsULST.raw.p=fit2$p.value[, "TLSTvsULST"],
  TLSTvsULST.fdr.p=p.adjust(fit2$p.value[, "TLSTvsULST"], "fdr"),
  TLSTvsULST.holm.p=p.adjust(fit2$p.value[, "TLSTvsULST"], "holm")
)

## generate heatmap for Ctrl, ULST and TLST
temp.exprs2 <- subset(exprs2,
    rownames(exprs2) %in% rownames(subset(res, res$ULSTvsCTRL.fdr.p < 0.05)))
temp.exprs2 <- temp.exprs2[, c(11:13, 18:25)]

pdf("tx_de_analysis_heatmap_CTRL_ULST_TLST.pdf", 8, 10)
drawHeatMap(temp.exprs2, cond[-c(4:7)])
dev.off()

## volcano plot of differential expression analysis of Ctrl vs ULST
temp.res <- res[, 1:2]
names(temp.res) <- c("logFC", "rawP")

pdf("tx_de_analysis_volcano_plots_CTRL_vs_ULST.pdf", 8, 10)
drawVolcanoPlot(temp.res)
dev.off()

## generate heatmap for Ctrl, ULST and TLST (GLUTATHIONE PATHWAY ALL)
temp.exprs2 <- exprs2[exprs2$EntrezGene %in% glut.tx, ]
rownames(temp.exprs2) <- temp.exprs2$Symbol
temp.exprs2 <- temp.exprs2[, c(11:13, 18:25)]

pdf("tx_de_analysis_heatmap_CTRL_ULST_TLST_glut_all.pdf", 8, 10)
drawHeatMap(temp.exprs2, cond[-c(4:7)])
dev.off()

## generate heatmap for Ctrl and ULST (GLUTATHIONE PATHWAY ALL)
temp.exprs3 <- exprs2[exprs2$EntrezGene %in% glut.tx, ]
rownames(temp.exprs3) <- temp.exprs3$Symbol
temp.exprs3 <- temp.exprs3[, c(11:13, 18:21)]

pdf("tx_de_analysis_heatmap_CTRL_ULST_glut_all.pdf", 8, 10)
drawHeatMap(temp.exprs3, cond[-c(4:7, 12:15)])
dev.off()

## generate heatmap for Ctrl and ULST (FIG 3E Gene Set)
temp.exprs4 <- exprs2[match(fig_3e_gene_set, exprs2$EntrezGene), ]
rownames(temp.exprs4) <- temp.exprs4$Symbol
temp.exprs4 <- temp.exprs4[, c(11:13, 18:21)]

pdf('tx_de_analysis_heatmap_CTRL_ULST_fig_3e_gene_set.pdf', 8, 6)
drawHeatMap(temp.exprs4, cond[-c(4:7, 12:15)])
dev.off()

## PCA
pdf("tx_glut_PCA.pdf", 8, 10)
PCAAnalysis(temp.exprs2, cond[-c(4:7)], conf.interval = 0.99)
dev.off()

## generate heatmap for Ctrl, ULST and TLST (GLUTATHIONE PATHWAY SIG ONLY)
temp.exprs2 <- exprs2[exprs2$EntrezGene %in% glut.tx, ]
temp.exprs2 <- subset(
    temp.exprs2,
    rownames(temp.exprs2) %in% rownames(
        subset(res, res$ULSTvsCTRL.fdr.p < 0.05)))
rownames(temp.exprs2) <- temp.exprs2$Symbol
temp.exprs2 <- temp.exprs2[, c(11:13, 18:25)]

pdf("tx_de_analysis_heatmap_CTRL_ULST_TLST_glut_sig.pdf", 8, 10)
drawHeatMap(temp.exprs2, cond[-c(4:7)])
dev.off()

## generate heatmap for Ctrl and ULST (GLUTATHIONE PATHWAY SIG ONLY)
temp.exprs3 <- exprs2[exprs2$EntrezGene %in% glut.tx, ]
temp.exprs3 <- subset(
    temp.exprs3,
    rownames(temp.exprs3) %in% rownames(
        subset(res, res$ULSTvsCTRL.fdr.p < 0.05)))
rownames(temp.exprs3) <- temp.exprs3$Symbol
temp.exprs3 <- temp.exprs3[, c(11:13, 18:21)]

pdf("tx_de_analysis_heatmap_CTRL_ULST_glut_sig.pdf", 8, 10)
drawHeatMap(temp.exprs3, cond[-c(4:7, 12:15)])
dev.off()

###############################################################################
source('pathway_specific_heatmaps.r')

heatmap.exprs <- exprs2[, c(1, 3, 11:13, 18:21)]
heatmap.groups <- cond[-c(4:7, 12:15)]

pdf("tx_de_analysis_heatmap_CTRL_ULST_glyserthr.pdf", 8, 10)
drawGeneHeatmapForPathway(heatmap.exprs, heatmap.groups, 'mmu00260')
dev.off()

pdf("tx_de_analysis_heatmap_CTRL_ULST_a-trna.pdf", 8, 10)
drawGeneHeatmapForPathway(heatmap.exprs, heatmap.groups, 'mmu00970')
dev.off()

pdf("tx_de_analysis_heatmap_CTRL_ULST_cysmet.pdf", 8, 10)
drawGeneHeatmapForPathway(heatmap.exprs, heatmap.groups, 'mmu00270')
dev.off()

pdf("tx_de_analysis_heatmap_CTRL_ULST_abc.pdf", 8, 10)
drawGeneHeatmapForPathway(heatmap.exprs, heatmap.groups, 'mmu02010')
dev.off()

pdf("tx_de_analysis_heatmap_CTRL_ULST_mineralabs.pdf", 8, 10)
drawGeneHeatmapForPathway(heatmap.exprs, heatmap.groups, 'mmu04978')
dev.off()

pdf("tx_de_analysis_heatmap_CTRL_ULST_tca.pdf", 8, 10)
drawGeneHeatmapForPathway(heatmap.exprs, heatmap.groups, 'mmu00020')
dev.off()

pdf("tx_de_analysis_heatmap_CTRL_ULST_centralcarbonmet.pdf", 8, 11)
drawGeneHeatmapForPathway(heatmap.exprs, heatmap.groups, 'mmu05230')
dev.off()

pdf("tx_de_analysis_heatmap_CTRL_ULST_nrf_response.pdf", 8, 10)
source('nrf2_response_genes.r')
drawGeneHeatmapForPathway(heatmap.exprs, heatmap.groups, genes = nrf2)
dev.off()
###############################################################################
###############################################################################
source('pathway_specific_barplots.r')

barplot.exprs <- exprs2[rownames(subset(res, res$ULSTvsCTRL.fdr.p < 0.05)),
                        c(1, 3, 11:13, 18:21)]
barplot.exprs <- cbind(barplot.exprs[, 1:2], 2**barplot.exprs[, 3:9])
barplot.groups <- cond[-c(4:7, 12:15)]

pdf("tx_de_analysis_barplot_CTRL_ULST_glyserthr.pdf", 8, 5)
print(drawGeneBarplotsForPathway(barplot.exprs, barplot.groups, 'mmu00260'))
dev.off()

pdf("tx_de_analysis_barplot_CTRL_ULST_a-trna.pdf", 8, 5)
print(drawGeneBarplotsForPathway(barplot.exprs, barplot.groups, 'mmu00970'))
dev.off()

pdf("tx_de_analysis_barplot_CTRL_ULST_cysmet.pdf", 8, 5)
print(drawGeneBarplotsForPathway(barplot.exprs, barplot.groups, 'mmu00270'))
dev.off()

pdf("tx_de_analysis_barplot_CTRL_ULST_abc.pdf", 8, 5)
print(drawGeneBarplotsForPathway(barplot.exprs, barplot.groups, 'mmu02010'))
dev.off()

pdf("tx_de_analysis_barplot_CTRL_ULST_mineralabs.pdf", 8, 5)
print(drawGeneBarplotsForPathway(barplot.exprs, barplot.groups, 'mmu04978'))
dev.off()

pdf("tx_de_analysis_barplot_CTRL_ULST_glut_all.pdf", 14, 5)
print(drawGeneBarplotsForPathway(barplot.exprs, barplot.groups, 'mmu00480'))
dev.off()
###############################################################################

## integrated human/mouse analysis - mouse heatmap
source("human_mouse_integrated_analysis.r")

temp.exprs2 <- exprs2[exprs2$EntrezGene %in% glut.tx, ]
temp.exprs2 <- subset(
    temp.exprs2,
    rownames(temp.exprs2) %in% rownames(
        subset(res, res$ULSTvsCTRL.fdr.p < 0.05)) &
    temp.exprs2$EntrezGene %in% orthologs$mouse.entrez.id)
rownames(temp.exprs2) <- temp.exprs2$Symbol
temp.exprs2 <- temp.exprs2[, c(11:13, 18:21)]

pdf("integrative_analysis_mouse_heatmap.pdf", 8, 10)
drawHeatMap(temp.exprs2, cond[-c(4:7, 12:15)])
dev.off()

stop('done.')

## compile and output DE results
res2 <- merge(exprs2[1], res, by = 0)
res2 <- res2[, -1]
write.csv(res2, file = "tx_de_analysis_summary.csv", row.names = F)

## pathway analysis
res3 <- merge(exprs2[3], res, by = 0)
res3 <- subset(res3, !is.na(res3$EntrezGene))
rownames(res3) <- res3[, 2]
res3 <- res3[, -c(1:2)]

## logistic regression with KEGG pathways
kegg <- KEGGLRpath(res3$ULSTvsCTRL.fdr.p, rownames(res3),
    member.type = "gene", species = "mmu")
kegg$pathway <- gsub(" - Mus musculus \\(mouse\\)", "", kegg$pathway)
kegg <- kegg[order(kegg$p.value), ]
write.csv(kegg[, -8], file = "pathway_analysis_kegg.csv", row.names = F)

## logistic regression with Reactome pathways
reactome <- ReactomeLRpath(res3$ULSTvsCTRL.fdr.p, rownames(res3),
    species = "Mus musculus")
reactome$pathway <- gsub("Mus musculus: ", "", reactome$pathway)
reactome <- reactome[order(reactome$p.value), ]
write.csv(reactome[, -7], file = "pathway_analysis_reactome.csv",
    row.names = F)

## hypergeometric test for pathway analysis using kegg
kegg2 <- ORAnalysis(rownames(res3), res3$ULSTvsCTRL.fdr.p, species = "mmu")
write.csv(kegg2[, -4], file = "pathway_analysis_ORA_kegg.csv",
    row.names = F)

## hypergeometric test for pathway analysis using reactome
reactome2 <- ORAnalysis(rownames(res3), res3$ULSTvsCTRL.fdr.p,
    db = "reactome", species = "Mus musculus")
write.csv(reactome2[, -4], file = "pathway_analysis_ORA_reactome.csv",
    row.names = F)
