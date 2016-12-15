library(reactome.db)
# library(CLEAN)

ReactomeLRpath <- function (sig.values, member.ids, min.g = 10, max.g = 99999,
                            sig.cutoff = 0.05, odds.min.max = c(0.001, 0.5),
                            species = "Homo sapiens")
{
    temp <- keys(reactome.db, "PATHNAME")
    temp <- temp[grepl(species, temp)]
    temp <- select(reactome.db, temp, c("PATHNAME", "ENTREZID"),
                       "PATHNAME")
    temp <- aggregate(ENTREZID ~ PATHNAME, data=temp, FUN="c")

    pathways <- temp$ENTREZID
    names(pathways) <- temp$PATHNAME

    cat("Running analysis")

    ENTREZ.in.DB <- unique(unlist(pathways))
    if (length(intersect(ENTREZ.in.DB, unique(member.ids))) < min.g)
      warning(paste("min.g=", min.g, " is larger than number of gene IDs overlapping with Reactome. Please consider lower min.g.",
        sep = ""))

    sig.values[sig.values == 0] <- 10^(-15)
    member.ids2 <- member.ids[!is.na(member.ids) & !is.na(sig.values) & member.ids %in% ENTREZ.in.DB]
    sig.values2 <- sig.values[!is.na(member.ids) & !is.na(sig.values) & member.ids %in% ENTREZ.in.DB]
    uniqids <- unique(member.ids2)
    numuniq <- length(uniqids)
    LOR.mult <- log(odds.min.max[2]) - log(odds.min.max[1])
    nlp <- (-1) * log(sig.values2)
    newp <- NA
    for (i in (1:numuniq)) {
      current <- nlp[member.ids2 == uniqids[i]]
      newp[i] <- mean(current)
    }
    catsizes <- sapply(pathways, length)
    yy <- pathways[catsizes >= min.g]
    ncats <- length(yy)
    siggenes <- uniqids[exp(-newp) < sig.cutoff]
    LRcoeff <- NA
    LRpval <- NA
    cattots <- NA
    yyind <- NA
    GOnums <- NA
    catsigIDs <- NA
    ind <- 0

    for (i in (1:ncats)) {
      catgenes <- as.character(yy[[i]])
      catgenes <- unique(catgenes)
      catpoprows <- match(catgenes, uniqids)
      catpoprows <- catpoprows[!is.na(catpoprows)]
      cattot <- length(catpoprows)
      if (cattot >= min.g & cattot <= max.g) {
        ind <- ind + 1
        cat <- rep(0, numuniq)
        cat[catpoprows] <- 1
        forLR <- as.data.frame(cbind(cat, newp))
        names(forLR) <- c("cat", "nlogp")
        withCallingHandlers(glm.lrGO <- glm(cat ~
          nlogp, family = binomial(link = "logit"),
          forLR), warning = function(w) {
          if (length(grep("fitted probabilities numerically 0 or 1 occurred",
            as.character(w))) == 1)
            warning("fitted probabilities numerically 0 or 1 occurred: functional category appears to induce linear separation of p-values.")
        })
        lrGO <- summary(glm.lrGO)
        LRcoeff[ind] <- lrGO$coefficients["nlogp",
          "Estimate"]
        LRpval[ind] <- lrGO$coefficients["nlogp",
          "Pr(>|z|)"]
        cattots[ind] <- cattot
        GOnums[ind] <- names(yy[i])
        catsig.IDlist <- intersect(siggenes, catgenes)
        catsigIDs[ind] <- paste(catsig.IDlist[order(catsig.IDlist)],
          collapse = ", ")
      }

      if (i%%100 == 0)
        cat(".")
    }

    cat("\n")

    odds.ratio <- exp(LOR.mult * LRcoeff)
    BenjFDR <- p.adjust(LRpval, "BH")

    if (all(is.na(GOnums))) {
      warning("no results found.")
      LRres <- NA
    }
    else {
      LRres <- cbind(GOnums, cattots,
        LRcoeff, odds.ratio, as.data.frame(LRpval),
        BenjFDR, catsigIDs)
      names(LRres) <- c("pathway",
        "n.members", "coeff", "odds.ratio", "p.value",
        "FDR", "sig.members")
      rownames(LRres) <- NULL
    }

    invisible(LRres)
}
