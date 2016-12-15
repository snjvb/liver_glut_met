library(bear)
library(KEGGREST)
library(reshape2)

getMetabolites <- function(pathway) {
    names(keggGet(pathway)[[1]]$COMPOUND)
}

getGenes <- function(pathway) {
    temp <- keggLink('genes', pathway)
    gsub('(.+?\\:)', '', temp)
}

generateDataForPlot <- function(gene.exprs, groups, control) {
    .gene.exprs <- apply(gene.exprs, 1, function(met) {
        scale <- tapply(met, groups, mean)[control]
        met / scale
    })

    .gene.exprs <- data.frame(Group = groups, .gene.exprs)
    .gene.exprs <- melt(.gene.exprs, id.vars = 'Group')

    df <- summarySE(.gene.exprs, measurevar = 'value', 
                     groupvars = c('Group', 'variable'))
    df$variable <- as.factor(df$variable)
    df$Group <- as.factor(df$Group)

    df
}

drawBarplots <- function(exprs, groups, control)
{
    df <- generateDataForPlot(exprs, groups, control)

    ggplot(data = df, aes(x = variable, y = value, fill = Group)) + 
        geom_bar(position = position_dodge(), 
                 stat = 'identity', color = 'black', size = 1.3, 
                 width = 0.6) +
        geom_errorbar(aes(ymin = ifelse(value < 0, value - sd, value), 
                          ymax = ifelse(value < 0, value, value + sd)), 
                      width = 0.2, size = 1, position = position_dodge(0.65)) +
        scale_y_sqrt() + 
        scale_fill_manual(values = c('black', 'white')) + 
        theme_bw() + 
        theme(
            legend.position = 'none', 
            axis.line = element_line(size = 1.3), 
            axis.ticks = element_line(size = 1.5), 
            axis.ticks.length = unit(3, 'mm'), 
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
            # axis.text.x = element_blank(), 
            axis.text.y = element_text(size = 20), 
            axis.title.x = element_blank(), 
            axis.title.y = element_blank(), 
            panel.border = element_blank(), 
            panel.grid = element_blank())
}

drawMetsBarplotsForPathway <- function(exprs, kegg.ids, groups, pathway) {
    temp <- merge(exprs, kegg.ids, by = 0)
    temp <- subset(temp, temp$KEGG %in% getMetabolites(pathway))

    rownames(temp) <- temp$BIOCHEMICAL
    temp <- temp[, 2:15]

    drawBarplots(temp, groups, 'LT2')
}

drawGeneBarplotsForPathway <- function(exprs, groups, pathway) {
    temp <- subset(exprs, exprs$EntrezGene %in% getGenes(pathway))

    rownames(temp) <- temp$Symbol
    temp <- temp[, -c(1, 2)]

    drawBarplots(temp, groups, 'CTRL')
}
