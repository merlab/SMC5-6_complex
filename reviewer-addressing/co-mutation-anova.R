# purpose: makes the dataframe for the alteration analysis of the complexes in the genes
source("./R/routine_tasks.R")
library(ggplot2)
library(viridis)
library(ggpubr)
# library(cooccur)
complexGenes <- c("NSMCE2", "SMC6", "SMC5", "NSMCE1", "NSMCE3", "NSMCE4A", "EID3")
instabilityGenes <- c(
    "TP53", "BRCA1", "BRCA2", "NBN", "TTK", "AURKA", "PLK1",
    "CHEK2", "CCNE1", "RB1", "RECQL4", "BLM"
)



generatePlotDf <- function(subcbpd) {
    d <- data.frame()
    for (j in instabilityGenes) {
        complexAlt <- as.numeric(subcbpd[, "isalt"])
        instabilityAlt <- as.numeric(subcbpd[, j])
        # shared <- sum(complexAlt == 1 & instabilityAlt == 1) / sum(complexAlt == 1)
        tp <- sum(complexAlt == 1 & instabilityAlt == 1)
        tn <- sum(complexAlt == 0 & instabilityAlt == 0)
        fn <- sum(complexAlt == 1 & instabilityAlt == 0)
        fp <- sum(complexAlt == 0 & instabilityAlt == 1)
        # print(j)
        # print(tp)
        # print(tn)
        # print(fn)
        # print(fp)
        # plotdf <- rbind(plotdf, c(gene = j, shared = shared, notshared = 1 - shared))
        d <- rbind(d, c(gene = j, tp = tp, tn = tn, fp = fp, fn = fn))
    }
    colnames(d) <- c("gene", "tp", "tn", "fp", "fn")
    d$tp <- as.numeric(d$tp)
    d$tn <- as.numeric(d$tn)
    d$fp <- as.numeric(d$fp)
    d$fn <- as.numeric(d$fn)

    # colnames(plotdf) <- c("gene", "shared", "notshared")
    # plotdf$shared <- as.numeric(plotdf$shared)
    # plotdf$notshared <- as.numeric(plotdf$notshared)


    return(d)
}

# # cbpd <- readRDS("./data/cbioportal/format_exOther.rds")
cbpd <- readRDS("./reviewer-addressing/cbioportal/format_exOther.rds")
instabilityGenes <- names(sort(colSums(cbpd[, instabilityGenes])))

d <- generatePlotDf(cbpd)
d[, 2:5] <- d[, 2:5] / rowSums(d[, 2:5])
d$ppv <- d$tp / (d$tp + d$fp)
d$FDR <- 1 - d$ppv
print(d)
print("done")
# cancerTypes <- c("All", unique(levels(cbpd$major)))
# cancerTypes <- cancerTypes[cancerTypes != "Other"]

# cancerTypes <- c(
#     "Breast",
#     "Ovarian",
#     "Prostate"
# )
# print(cancerTypes)
# for (cancerType in cancerTypes) {
#     print(cancerType)
# subcbpd <- cbpd[cbpd$major == cancerType, ]
# subcbpd <- subcbpd[, c("isalt", instabilityGenes)]
