library(KEGGREST)

glut.tx <- keggLink("genes", "mmu00480")
glut.tx <- gsub("mmu\\:", "", glut.tx)

glut.mx <- names(keggGet("mmu00480")[[1]]$COMPOUND)
