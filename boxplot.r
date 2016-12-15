library(ggplot2)
library(reshape2)

moi <- c('C00051', 'C00127')
# moi <- c('C00072', 'C05422')
boxplot.exprs2 <- boxplot.exprs[boxplot.exprs$KEGG %in% moi, ]
data <- melt(data.frame(Metabolite = boxplot.exprs2$BIOCHEMICAL, 
                        boxplot.exprs2[, 3:16],
                        row.names = NULL,  
                        stringsAsFactors = F), 
             id.vars = c('Metabolite'), value.name = "Expression")
data$variable <- gsub("(.+?)_.*", "\\1", data$variable)

sig <- data.frame(Metabolite = boxplot.exprs2$BIOCHEMICAL, 
                  FDR = boxplot.exprs2$fdr, 
                  row.names = NULL, 
                  stringsAsFactors = F)
sig$Y <- apply(boxplot.exprs2[, 3:16], 1, function(x) max(x) + 1.1)
sig$Symbol <- as.character(symnum(sig$FDR, corr = FALSE, na = FALSE, 
                     cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                     symbols = c("***", "**", "*", ".", " ")))

path <- data.frame(Metabolite = rep(boxplot.exprs2$BIOCHEMICAL, each = 4), 
                   X = c(1, 1, 2, 2), 
                   Y = rep(apply(boxplot.exprs2[, 3:16], 1, max), each = 4) + 
                   c(0.5, 0.9, 0.9, 0.5))

png("gsh_gssg_boxplot.png", 640, 320, 'px')
print(ggplot(data = data, aes(x = factor(variable), y = Expression)) + 
    geom_boxplot(outlier.size = 0, color = c('red', 'blue'), 
        fill = 'lightgrey', size = 1.1) + 
    geom_jitter(position = position_jitter(), size = 4) + 
    geom_text(data = sig, aes(x = 1.5, y = Y, label = Symbol), size = 14) + 
    geom_path(data = path, aes(x = X, y = Y), size = 1.5) + 
    xlab("") + 
    ylab("Log2 Abundance") + 
    ylim(min(boxplot.exprs2[, 3:16]) - 0.2, max(boxplot.exprs2[, 3:16]) + 2) + 
    facet_grid(. ~ Metabolite) + 
    theme_bw() + 
    theme(
        axis.text.y = element_text(size = 20), 
        axis.title.y = element_text(size = 24, face = 'bold'), 
        axis.text.x = element_text(size = 30), 
        strip.text.x = element_text(size = 18, face = 'bold')))
dev.off()
