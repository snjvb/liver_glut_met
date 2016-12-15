orthologs.file <- "~/Documents/brittany/myc_induction_regression_transcriptomics_062514/HOM_MouseHumanSequence.rpt"

## get mouse glutathione pathway sig DE genes

genes <- c('14775', '14870', '14859', '56615', '14857', '14862', '14864', '110208', '269951', '18263', '20810', '68312', '76263', '14873', '14629', '14866', '66447', '14867', '14863', '14858', '15926', '14598', '20603', '14872', '75475', '14778', '14380', '14860', '14381', '20133', '211666', '20135')

## prepare orthologs dataframe (sourced from MGI Informatics)
orthologs <- read.delim(orthologs.file, as.is = T, 
    fill = T)[, c(1, 2, 4, 5)]
orthologs <- split(orthologs[, -2], orthologs$Common.Organism.Name)
names(orthologs) <- c("human", "mouse")
orthologs <- merge(orthologs$mouse, orthologs$human, by = "HomoloGene.ID")
names(orthologs) <- c("homologene.id", "mouse.symbol", "mouse.entrez.id", 
    "human.symbol", "human.entrez.id")

## subset orthologs to only glutathione pathway genes
orthologs <- subset(orthologs, orthologs$mouse.entrez.id %in% genes)

## manually filter out duplicate orthologs
orthologs <- orthologs[-which(rownames(orthologs) %in% c("11109", '14913')), ]
