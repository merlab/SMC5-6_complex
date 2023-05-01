row <- "isalt"
library(limma)
library(ggplot2)
library(ggpubr)
library(viridis)
instabilityT <- 3
source("./R/routine_tasks.R")
doLimmma <- function(cbpd, design) {
    fit <- lmFit(expmat, design)
    fitted.ebayes <- eBayes(fit)
    out <- topTable(fitted.ebayes, number = Inf)
    out$logFDR <- -log10(out$adj.P.Val)
    out$gene <- rownames(out)
    out$Sig <- "Not Significant"
    out$Sig[out$adj.P.Val < 0.01] <- "Significant"
    return(out)
}
expmat <- readRDS("./data/metabric-brca/microarray-metagx.rds")
samples <- colnames(expmat)
cbpd <- readRDS("./reviewer-addressing/cbioportal/format_exOther.rds")
rownames(cbpd) <- cbpd$name
cbpd <- cbpd[cbpd$name %in% samples, ]
expmat <- expmat[, colnames(expmat) %in% cbpd$name]
cbpd <- cbpd[colnames(expmat), ]
expmat <- expmat[, cbpd$name]

cbpd$isalt <- as.numeric(cbpd$isalt)
cbpd$isalt <- ifelse(cbpd$isalt == 1, "Mutated", "Wild")
cbpd$isalt <- factor(cbpd$isalt, levels = c("Wild", "Mutated"))

cbpd$instabilityScore <- as.numeric(cbpd$instabilityScore)
# cbpd$instabilityScore <- ifelse(cbpd$instabilityScore >= 3, "Unstable", "Stable")
cbpd$instabilityScore <- ifelse(cbpd$instabilityScore >= instabilityT, "Unstable", "Stable")
cbpd$instabilityScore <- factor(cbpd$instabilityScore, levels = c("Stable", "Unstable"))
# complex
er <- as.numeric(cbpd$isalt)
design <- model.matrix(~er)
design <- as.matrix(design)
row.names(design) <- cbpd$name
complexRes <- doLimmma(expmat, design)
# instabiltiy
er <- as.numeric(cbpd$instabilityScore)
design <- model.matrix(~er)
design <- as.matrix(design)
row.names(design) <- cbpd$name
instabilityRes <- doLimmma(expmat, design)
mergeRes <- merge(complexRes, instabilityRes, by = "gene")
mergeRes$bothSig <- "Not Significant"
mergeRes$bothSig[mergeRes$Sig.x == "Not Significant" & mergeRes$Sig.y == "Significant"] <- "Instability Only"
mergeRes$bothSig[mergeRes$Sig.x == "Significant" & mergeRes$Sig.y == "Not Significant"] <- "Complex Only"
mergeRes$bothSig[mergeRes$Sig.x == "Significant" & mergeRes$Sig.y == "Significant"] <- "Both"
# x is complex, y is instability
pdf(sprintf("./reviewer-addressing/lfc-%s.pdf", instabilityT), height = 6, width = 6)
plot(
    ggplot(mergeRes, aes(x = logFC.x, y = logFC.y, color = bothSig)) +
        geom_point(alpha = 0.5) + # , size = .25) +
        scale_color_viridis(discrete = TRUE) +
        xlab("Complex LFC") +
        ggtitle(paste("Instability >=", instabilityT)) +
        ylab("Instability LFC") +
        theme_classic()
)
dev.off()
print("done")

# saveRDS(out, "./data/metabric-brca/dgea-limma.rds")
