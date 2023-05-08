row <- "isalt"
library(limma)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(viridis)
library(writexl)
library(DEGreport)
instabilityT <- 3
source("./R/routine_tasks.R")
doLimmma <- function(cbpd, design) {
    fit <- lmFit(expmat, design)
    fitted.ebayes <- eBayes(fit)
    out <- topTable(fitted.ebayes, number = Inf)
    out$logFDR <- -log10(out$adj.P.Val)
    out$gene <- rownames(out)
    out$Sig <- "Not Significant"
    out$Sig[out$adj.P.Val < 0.01 & abs(out$logFC) > 0.15] <- "Significant"
    return(out)
}
expmat <- readRDS("./data/metabric-brca/microarray-metagx.rds")
samples <- colnames(expmat)
cbpd <- readRDS("./reviewer-addressing/cbioportal/format_exOther.rds")
# t <- cbpd[cbpd$isalt == 1, ]
# t <- t[t$study == "Breast Invasive Carcinoma (TCGA, PanCancer Atlas)", ]
# stop()
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

show_genes <- c(
    "RAD54B", "CCNE2", "RECQL4", "MCM4", "AURKA", #' NSMCE2',
    "KIF13B", "PPP2R2A", "PARP3", "TP53BP1"
)

df2 <- mergeRes[mergeRes$gene %in% show_genes, ]

# x is complex, y is instability
set.seed(123)
pdf(sprintf("./reviewer-addressing/plot/lfc-%s.pdf", instabilityT), height = 6, width = 6)
plot(
    ggplot(mergeRes, aes(x = logFC.x, y = logFC.y, color = bothSig, alpha = bothSig)) +
        geom_point() + # , size = .25) + alpha = 0.5
        # scale_color_viridis(discrete = TRUE) +
        # geom_label_repel(
        #     data = df2, mapping = aes(label = gene),
        #     color = "red",
        #     size = 4, min.segment.length = 0,
        #     nudge_x = 1, nudge_y = 1,
        #     label.size = NA, fill = NA
        # ) +
        geom_text_repel(
            data = df2[df2$logFC.x < 0, ], mapping = aes(label = gene),
            size = 4, # min.segment.length = 0,
            label.size = NA, fill = NA, color = "red",
            nudge_x = 0,
            xlim = c(-1.12, -1.11),
            alpha = 1
        ) +
        geom_text_repel(
            data = df2[df2$logFC.x > 0, ], mapping = aes(label = gene),
            size = 4, # min.segment.length = 0 # 1.75
            label.size = NA, fill = NA, , color = "red", # nudge_x = -0.1,
            # xlim = c(0.95, 1)
            # nudge_x = 0,
            xlim = c(1.05, 1.06),
            alpha = 1
        ) +
        scale_color_manual(
            values = c("#d95f02", "#1b9e77", "#7570b3", "grey50"),
            labels = c("Both", "SMC5-6 Complex", "Instability Score", "Not Significant")
        ) + # "grey70"
        # scale_alpha_manual(values = c(0.85, 0.5, 0.5, 0.1)) + # 0.1
        scale_alpha_manual(values = c(0.85, 0.4, 0.4, 0.08)) + # 0.1
        labs(
            color = "",
            y = "Instability LFC",
            x = "Complex LFC"
        ) +
        scale_x_continuous(
            limits = c(-1.65, 1.4),
            breaks = c(-1.5, -1, -.5, 0, .5, 1, 1.5)
        ) +
        scale_y_continuous(
            limits = c(-1.3, .8),
            breaks = c(-1.25, -1, -.75, -.5, -.25, 0, .25, .5, .75)
        ) +
        geom_segment(
            x = -.15, xend = .15, y = .15, yend = .15,
            color = "black", linetype = "dashed"
        ) +
        geom_segment(
            x = -.15, xend = .15, y = -.15, yend = -.15,
            color = "black", linetype = "dashed"
        ) +
        geom_segment(
            x = -.15, xend = -.15, y = .15, yend = -.15,
            color = "black", linetype = "dashed"
        ) +
        geom_segment(
            x = .15, xend = .15, y = -.15, yend = 0.15,
            color = "black", linetype = "dashed"
        ) +
        guides(color = guide_legend(override.aes = list(size = 4)), alpha = "none") +
        theme_classic() +
        theme(
            legend.position = "top",
            legend.text = element_text(size = 10),
            legend.direction = "horizontal"
        )
)
dev.off()

colnames(mergeRes) <- gsub("\\.x$", "-SMC5-6", colnames(mergeRes))
colnames(mergeRes) <- gsub("\\.y$", "-instability", colnames(mergeRes))
print(dim(mergeRes))
mergeRes <- mergeRes[, -grep("^t-", colnames(mergeRes))]
print(dim(mergeRes))
mergeRes <- mergeRes[, -grep("^logFDR", colnames(mergeRes))]
print(dim(mergeRes))
mergeRes <- mergeRes[, -grep("^B-", colnames(mergeRes))]
print(dim(mergeRes))
mergeRes <- mergeRes[, -grep("^Sig-", colnames(mergeRes))]
print(dim(mergeRes))
dim(mergeRes[mergeRes$bothSig == "Both", ])
write_xlsx(mergeRes, "./reviewer-addressing/METABRIC-lfc-analysis.xlsx")
print("done")
