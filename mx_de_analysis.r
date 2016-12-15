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

source("~/Documents/brittany/myc_induction_regression_transcriptomics_062514/functions/FishersExactTest.R")
source("~/Documents/brittany/myc_induction_regression_transcriptomics_062514/functions/PCA.R")
source("~/Documents/brittany/myc_induction_regression_transcriptomics_062514/gene_set.r")

exprs.file <- "~/Documents/brittany/myc_induction_regression_transcriptomics_062514/raw_data/mx_exprs.txt"

## read in expression values
raw <- read.delim(exprs.file, as.is = T, fill = T, check.names = F,
    row.names = 5)[, c(1, 9, 11:40)]

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

## extract and clean KEGG IDs map
kegg.ids <- raw[, 1:2]
kegg.ids$KEGG[137] <- "C00268"
kegg.ids$KEGG[197] <- "C00263"
kegg.ids$KEGG[333] <- "C05411"
kegg.ids <- kegg.ids[!duplicated(kegg.ids), ]

raw <- raw[, -c(1:2)]

## extract group information
groups <- gsub("(.*)\\_.*", "\\1", colnames(raw))

## convert all expression values to numerics
raw2 <- apply(raw, 2, function(col) {
    col <- as.numeric(gsub(",", "", col))
})
rownames(raw2) <- rownames(raw)

## replace missing values with the minimum expression value for the
## corresponding metabolite
raw3 <- t(apply(raw2, 1, function(row) {
    row[is.na(row)] <- min(row, na.rm = T)

    row
}))

## log transform data matrix to prepare for limma analysis
exprs <- log2(raw3)

## limma DE analysis
design <- model.matrix(~0 + groups)
colnames(design) <- gsub("groups", "", colnames(design))

fit <- lmFit(exprs, design = design)
contrasts <- makeContrasts(MycEffect = MYC - LT2, levels = design)
fit2 <- contrasts.fit(fit, contrasts = contrasts)
fit2 <- eBayes(fit2)

res <- data.frame(
  logFC=fit2$coefficients[, "MycEffect"],
  rawP=fit2$p.value[, "MycEffect"],
  fdr=p.adjust(fit2$p.value[, "MycEffect"], "fdr"),
  holm=p.adjust(fit2$p.value[, "MycEffect"], "holm")
)

## generate heatmap for LT2 vs MYC
temp.exprs <- subset(exprs,
    rownames(exprs) %in% rownames(subset(res, res$fdr < 0.05)))
temp.exprs <- temp.exprs[, c(1:14)]

pdf("mx_de_analysis_heatmap_LT2_MYC.pdf", 8, 10)
drawHeatMap(temp.exprs, groups[1:14])
dev.off()

# volcano plot of differential expression analysis of LT2 vs MYC
pdf("mx_de_analysis_volcano_plots_LT2_MYC.pdf", 8, 10)
drawVolcanoPlot(res)
dev.off()

## compile and output DE results
res2 <- merge(kegg.ids[1], res, by = 0)
res2 <- res2[, -1]
res2 <- res2[order(res2$rawP), ]
write.csv(res2, file = "mx_de_analysis_summary.csv", row.names = F)

## merge KEGG IDs and clean up duplicates
res3 <- merge(kegg.ids[2], res, by = 0)
res3 <- subset(res3, res3$KEGG != "")
res3 <- res3[order(res3$rawP), ]
res3 <- res3[!duplicated(res3$KEGG), ]
rownames(res3) <- res3$KEGG
res3 <- res3[, -c(1, 2)]

## generate dataset for subsequent generation of boxplots
boxplot.exprs <- merge(exprs, kegg.ids, by = 0)
boxplot.exprs <- subset(boxplot.exprs, boxplot.exprs$KEGG != '')
boxplot.exprs <- merge(boxplot.exprs, res3, by.x = 'KEGG', by.y = 0)

## generate heatmap for LT2, MYC (GLUTATHIONE PATHWAY ALL)
temp.exprs <- merge(exprs, kegg.ids, by = 0)
temp.exprs <- subset(temp.exprs, temp.exprs$KEGG %in% glut.mx)
rownames(temp.exprs) <- temp.exprs$BIOCHEMICAL
temp.exprs <- temp.exprs[, 2:15]

pdf("mx_de_analysis_heatmap_LT2_MYC_glut_all.pdf", 8, 10)
drawHeatMap(temp.exprs, groups[1:14])
dev.off()

## PCA
pdf("mx_glut_PCA.pdf", 8, 10)
PCAAnalysis(temp.exprs, groups[1:14], conf.interval = 0.99)
dev.off()

## generate heatmap for LT2, MYC (GLUTATHIONE PATHWAY SIG ONLY)
temp.exprs <- merge(exprs, kegg.ids, by = 0)
temp.exprs <- subset(temp.exprs,
    temp.exprs$KEGG %in% glut.mx &
    temp.exprs$KEGG %in% rownames(res3[res3$fdr < 0.05, ]))
rownames(temp.exprs) <- temp.exprs$BIOCHEMICAL
temp.exprs <- temp.exprs[, 2:15]

pdf("mx_de_analysis_heatmap_LT2_MYC_glut_sig.pdf", 8, 10)
drawHeatMap(temp.exprs, groups[1:14])
dev.off()

###############################################################################
source('pathway_specific_heatmaps.r')

heatmap.exprs <- exprs
heatmap.groups <- groups[1:14]

pdf("mx_de_analysis_heatmap_LT2_MYC_glyserthr.pdf", 8, 10)
drawMetsHeatmapForPathway(heatmap.exprs, kegg.ids, heatmap.groups, 'mmu00260')
dev.off()

pdf("mx_de_analysis_heatmap_LT2_MYC_a-trna.pdf", 8, 10)
drawMetsHeatmapForPathway(heatmap.exprs, kegg.ids, heatmap.groups, 'mmu00970')
dev.off()

pdf("mx_de_analysis_heatmap_LT2_MYC_cysmet.pdf", 8, 10)
drawMetsHeatmapForPathway(heatmap.exprs, kegg.ids, heatmap.groups, 'mmu00270')
dev.off()

pdf("mx_de_analysis_heatmap_LT2_MYC_abc.pdf", 8, 10)
drawMetsHeatmapForPathway(heatmap.exprs, kegg.ids, heatmap.groups, 'mmu02010')
dev.off()

pdf("mx_de_analysis_heatmap_LT2_MYC_mineralabs.pdf", 8, 10)
drawMetsHeatmapForPathway(heatmap.exprs, kegg.ids, heatmap.groups, 'mmu04978')
dev.off()
###############################################################################
###############################################################################
source('pathway_specific_barplots.r')

barplot.exprs <- raw3[rownames(subset(res, res$fdr < 0.05)), ]
barplot.groups <- groups[1:14]

pdf("mx_de_analysis_barplot_LT2_MYC_glyserthr.pdf", 8, 5)
print(drawMetsBarplotsForPathway(barplot.exprs, kegg.ids, barplot.groups,
                                 'mmu00260'))
dev.off()

pdf("mx_de_analysis_barplot_LT2_MYC_a-trna.pdf", 8, 5)
print(drawMetsBarplotsForPathway(barplot.exprs, kegg.ids, barplot.groups,
                                 'mmu00970'))
dev.off()

pdf("mx_de_analysis_barplot_LT2_MYC_cysmet.pdf", 8, 5)
print(drawMetsBarplotsForPathway(barplot.exprs, kegg.ids, barplot.groups,
                                 'mmu00270'))
dev.off()

pdf("mx_de_analysis_barplot_LT2_MYC_abc.pdf", 8, 5)
print(drawMetsBarplotsForPathway(barplot.exprs, kegg.ids, barplot.groups,
                                 'mmu02010'))
dev.off()

pdf("mx_de_analysis_barplot_LT2_MYC_mineralabs.pdf", 8, 5)
print(drawMetsBarplotsForPathway(barplot.exprs, kegg.ids, barplot.groups,
                                 'mmu04978'))
dev.off()

pdf("mx_de_analysis_barplot_LT2_MYC_glutmet.pdf", 8, 5)
print(drawMetsBarplotsForPathway(barplot.exprs, kegg.ids, barplot.groups,
                                 'mmu00480'))
dev.off()
###############################################################################

## hypergeometric test for pathway analysis using kegg
kegg <- ORAnalysis(rownames(res3), res3$fdr, type = "COMPOUND",
    species = "mmu")
write.csv(kegg[, -4], file = "mx_pathway_analysis_ORA_kegg.csv",
    row.names = F)
