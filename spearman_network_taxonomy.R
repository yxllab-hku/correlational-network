print("reading count table")
dat <- read.csv("df_count.csv")
taxLevels <- c(
# can be modified to include more taxonomy levels
"genus"
)
for (taxLevel in taxLevels) {
    print(paste0("Analyzing taxonomy level: ", taxLevel))
    precursor <- NULL
    peptidase <- NULL
    genomeInTaxonmy <- NULL
    taxonomy <- NULL
    spearman <- NULL
    spearmanP <- NULL
    p <- NULL
    odds <- NULL
    precursorCount <- NULL
    peptidaseCount <- NULL
    precursorContainingGenomes <- NULL
    peptidaseContainingGenomes <- NULL
    intersectionContainingGenomes <- NULL
    outputIndex <- 1

    is_valid <- function(x) ifelse(is.integer(x), length(unique(x)) != 1, TRUE)

    for (tax in as.list(t(unique(dat[taxLevel])))) {
        subDat <- dat[dat[taxLevel] == tax,]
        # remove column with all values > 0 or all values = 0
        subDat <- subDat[, sapply(subDat, is_valid)]
        # find columns of precursors and peptidases
        precursor_indexes <- which(grepl("pre_", colnames(subDat)))
        peptidase_indexes <- which(grepl("peptidase_", colnames(subDat)))
        if (length(precursor_indexes) == 0 | length(peptidase_indexes) == 0) { next }
        print(paste0("analyzing ", taxLevel, ": ", tax))
        for (i in precursor_indexes) {
            for (j in peptidase_indexes) {
                precursorContainingGenome <- as.integer(as.logical(subDat[, i]))
                peptidaseContainingGenome <- as.integer(as.logical(subDat[, j]))
                if (length(unique(precursorContainingGenome)) != 1 & length(unique(peptidaseContainingGenome)) != 1){
                    f <- fisher.test(precursorContainingGenome, peptidaseContainingGenome, alternative = "greater")
                    p[outputIndex] <- f$p
                    odds[outputIndex] <- f$estimate
                } else {
                    p[outputIndex] <- NA
                    odds[outputIndex] <- NA
                }
                s <- cor.test(subDat[, i], subDat[, j], method = "spearman", alternative = "greater", exact = FALSE)
                spearman[outputIndex] <- s$estimate
                spearmanP[outputIndex] <- s$p.value
                precursor[outputIndex] <- names(subDat[i])
                peptidase[outputIndex] <- names(subDat[j])
                precursorCount[outputIndex] <- sum(subDat[, i])
                peptidaseCount[outputIndex] <- sum(subDat[, j])
                intersectionContainingGenomes[outputIndex] <- sum(precursorContainingGenome & peptidaseContainingGenome)
                precursorContainingGenomes[outputIndex] <- sum(precursorContainingGenome)
                peptidaseContainingGenomes[outputIndex] <- sum(peptidaseContainingGenome)
                genomeInTaxonmy[outputIndex] <- nrow(subDat)
                taxonomy[outputIndex] <- tax
                outputIndex <- outputIndex + 1
            }
        }
    }

    pAdjusted <- p.adjust(p, "BH") #Benjamini–Hochberg procedure, false discovery rate
    spearmanPAdjusted <- p.adjust(spearmanP, "BH")
    resultTable <- data.frame(taxonomy, precursor, peptidase, spearman, spearmanP, spearmanPAdjusted, intersectionContainingGenomes, precursorCount, peptidaseCount, precursorContainingGenomes, peptidaseContainingGenomes, genomeInTaxonmy)
    write.csv(resultTable, paste0("correlation_spearman_", taxLevel, ".csv"), row.names = FALSE, quote = FALSE)
    resultTable <- read.csv(paste0("correlation_spearman_", taxLevel, ".csv"))

    Fig1_1 <- resultTable[(resultTable$spearman > 0) & (resultTable$taxonomy != ""),]
    write.csv(Fig1_1, paste0("Fig1_1_", taxLevel, ".csv"), row.names = FALSE, quote = FALSE)

    # generate correlation netowrk, node, edge tables
    networkTable <- resultTable[(resultTable$spearman > 0.3) & (resultTable$spearmanPAdjusted < 0.00001) & (resultTable$intersectionContainingGenomes >=10) & (resultTable$taxonomy != ""),]
    networkTable <- networkTable[rowSums(is.na(networkTable)) != ncol(networkTable),] # remove NA rows
    networkTable <- networkTable[networkTable$taxonomy != "", ] # remove empty taxanomy
    edgeTable <- networkTable
    if (nrow(networkTable) > 0) {
        print("creating network tables")
        networkTable$interaction <- "|"
        nodeTablePre <- unique(data.frame(nodes = networkTable$precursor, number = networkTable$precursorCount))
        nodeTablepeptidase <- unique(data.frame(nodes = networkTable$peptidase, number = networkTable$peptidaseCount))
        nodeTablePre$color <- "#ff6666"
        nodeTablepeptidase$color <- "#3399ff"
        nodeTablePre$nodeSize <- apply(nodeTablePre['number'], 1, function(x) log(x, 10) * 20)
        nodeTablepeptidase$nodeSize <- apply(nodeTablepeptidase['number'], 1, function(x) log(x, 10) * 15)
        nodeTable <- rbind(nodeTablePre, nodeTablepeptidase)
        edgeTable["shared name"] <- paste0(edgeTable$precursor, " (|) ", edgeTable$peptidase)
        write.csv(networkTable, paste0("network_spearman_", taxLevel, ".csv"), row.names = FALSE, quote = FALSE)
        write.csv(nodeTable, paste0("node_spearman_", taxLevel, ".csv"), row.names = FALSE, quote = FALSE)
        write.csv(edgeTable, paste0("edge_spearman_", taxLevel, ".csv"), row.names = FALSE, quote = FALSE)

    } else {
        print("no significant correlation found, no network will be created")
    }
}
print("finished")
