library(combinat)
library(KEGGREST)
library(reactome.db)

calcFishersP <- function(M, N, k, r) {
    pval <- 0
    for (i in r:k) {
        pval <- pval + dhyper(i, M, N, k)
    }

    pval
}

reactomePathways <- function(member.ids, species = "Homo sapiens") {
    species.prefix <- paste(species, ":", sep = "")

    p <- keys(reactome.db, "PATHNAME")
    p <- p[grepl(species.prefix, p)]
    p <- select(reactome.db, p, c("PATHNAME", "ENTREZID"), "PATHNAME")
    p$PATHNAME <- gsub(species.prefix, "", p$PATHNAME)

    p <- aggregate(p$ENTREZID ~ p$PATHNAME, FUN = "as.vector")
    p2 <- p[, 2]
    names(p2) <- p[, 1]

    p2
}

keggPathways <- function(member.ids, type = "GENE", species = "hsa") {
    pathways <- keggList("pathway", species)
    pathway.ids <- gsub("path:", "", names(pathways))
    pathway.names <- gsub("(.*) - .+$", "\\1", pathways)

    pathway.members <- lapply(pathway.ids, function(p) {
        tryCatch({
            temp <- keggGet(p)[[1]][[type]]
            stopifnot(!is.null(temp))
            if (type == "GENE") {
                intersect(temp[seq(1, length(temp), 2)], member.ids)
            } else if (type == "COMPOUND") {
                intersect(names(temp), member.ids)
            }
        }, error = function(err) {
            NULL
        })
    })
    names(pathway.members) <- paste(pathway.names, " (", pathway.ids, ")", 
        sep = "")

    pathway.members
}

ORAnalysis <- function(member.ids, sig.values, sig.cutoff = 0.05, 
    db = "kegg", ...) {
    all.members <- member.ids
    sig.members <- member.ids[sig.values < sig.cutoff]

    if (db == "kegg") {
        pathways <- keggPathways(member.ids, ...)
    } else if (db == "reactome") {
        pathways <- reactomePathways(member.ids, ...)
    }

    stats <- lapply(names(pathways), function(p) {
        p.genes <- pathways[[p]]
        p.sig <- intersect(p.genes, sig.members)

        M <- length(p.genes)
        N <- length(all.members) - M
        k <- length(sig.members)
        r <- length(p.sig)

        list(
            Pathway = p, 
            Total.Members = M, 
            Sig.Members = r, 
            Sig.IDs = paste(p.sig, collapse = ", "), 
            p.value = calcFishersP(M, N, k, r))
    })

    res <- as.data.frame(do.call("rbind", lapply(stats, unlist)), 
        stringsAsFactors = F)
    res$p.value <- as.numeric(res$p.value)
    res$FDR <- p.adjust(res$p.value, "BH")

    res[order(res$p.value), ]
}
