library(KEGGREST)
# library(CLEAN)

GetPathwaysForSpecies <- function(species) {
  pathways <- keggList("pathway", species)
  pathway.ids <- gsub("path:", "", names(pathways))
  pathways <- data.frame(
    ID = pathway.ids,
    Name = pathways,
    row.names = pathway.ids,
    stringsAsFactors = FALSE
  )
  as.matrix(pathways)
}

GetCompoundsForPathway <- function(pathway, members) {
  compounds <- tryCatch({
    temp <- keggGet(pathway)[[1]]$COMPOUND
    stopifnot(!is.null(temp))
    intersect(names(temp), members)
  }, error = function(err) {
    NULL
  })
  compounds
}

GetGenesForPathway <- function(pathway, members) {
  genes <- tryCatch({
    temp <- keggGet(pathway)[[1]]$GENE
    stopifnot(!is.null(temp))
    intersect(temp[seq(1, length(temp), 2)], members)
  }, error = function(err) {
    NULL
  })
  genes
}

GetMembersForPathway <- function(pathway, members) {
  members <- tryCatch({
    kegg.record <- keggGet(pathway)[[1]]
    genes <- kegg.record$GENE[seq(1, length(kegg.record$GENE), 2)]
    compounds <- names(kegg.record$COMPOUND)
    intersect(append(genes, compounds), members)
  }, error = function(err) {
    NULL
  })
  members
}

KEGGLRpath <- function (sig.values, member.ids, member.type = "both",
                        min.g = 10, max.g = 99999, sig.cutoff = 0.05, 
                        odds.min.max = c(0.001, 0.5), species = "hsa")
{
  if (is.factor(member.ids))
    member.ids <- as.character(member.ids)

  if (!member.type %in% c("both", "compound", "gene"))
    stop("Invalid member type or member type not specified.")

  cat(paste("Fetching pathways for species: ", species, "...\n", sep = ""))
  pathway.info <- GetPathwaysForSpecies(species)


  if (member.type == "compound") {
    f <- GetCompoundsForPathway
  } else if (member.type == "gene") {
      f <- GetGenesForPathway
  } else {
      f <- GetMembersForPathway
  }

  cat(paste("Building pathway sets...\n", sep = ""))
  pathway.members <- sapply(pathway.info[, "ID"], f, member.ids)

  cat("Running analysis")

  ENTREZ.in.DB <- unique(unlist(pathway.members))
  if (length(intersect(ENTREZ.in.DB, unique(member.ids))) < min.g)
    warning(paste("min.g=", min.g, " is larger than number of gene IDs
            overlapping with KEGG. Please consider lower min.g.",
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
  catsizes <- sapply(pathway.members, length)
  yy <- pathway.members[catsizes >= min.g]
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
          warning("fitted probabilities numerically 0 or 1 occurred:
                  functional category appears to induce linear separation of
                  p-values.")
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
  kterms <- pathway.info[, 2]
  k.ids <- pathway.info[, 1]
  keggrows <- match(GOnums, k.ids)
  kegg.annot <- kterms[keggrows]

  if (all(is.na(GOnums))) {
    warning("no results found.")
    LRres <- NA
  }
  else {
    LRres <- cbind(GOnums, kegg.annot, cattots,
      LRcoeff, odds.ratio, as.data.frame(LRpval),
      BenjFDR, catsigIDs)
    names(LRres) <- c("pathway.id", "pathway",
      "n.members", "coeff", "odds.ratio", "p.value",
      "FDR", "sig.members")
    rownames(LRres) <- NULL
  }

  invisible(LRres)
}
