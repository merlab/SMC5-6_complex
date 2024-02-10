row <- "isalt"
library(limma)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(viridis)
library(writexl)
library(DEGreport)
# instabilityT <- 2
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
# cbpd$instabilityScore <- ifelse(cbpd$instabilityScore > instabilityT, "Unstable", "Stable")
cbpd$instabilityScore <- factor(cbpd$instabilityScore, levels = c("Stable", "Unstable"))
# print(table(cbpd$instabilityScore))
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
p <- (
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
    geom_point(data = df2[df2$logFC.x > 0, ], size = 0.75, color = "red", fill = "red", shape = 15) +
    geom_point(data = df2[df2$logFC.x < 0, ], size = 0.75, color = "#1966ffff", fill = "#1966ffff", shape = 15) +
    geom_text_repel(
      data = df2[df2$logFC.x < 0, ], mapping = aes(label = gene),
      size = 4, # min.segment.length = 0,
      label.size = NA, fill = NA, color = "#1966ffff",
      nudge_x = 0,
      xlim = c(-1.12, -1.11),
      alpha = 1
    ) +
    geom_text_repel(
      data = df2[df2$logFC.x > 0, ], mapping = aes(label = gene),
      size = 4, # min.segment.length = 0 # 1.75
      label.size = NA, fill = NA, , color = "red",
      # xlim = c(0.95, 1)
      # nudge_x = 0,
      xlim = c(1.05, 1.06),
      alpha = 1
    ) +
    scale_color_manual(
      values = c("#d95f02", "#1b9e77", "#7570b3", "grey50"),
      # labels = c("Both", "SMC5-6 Complex", "Instability Score", "Not Significant")
      labels = c(
        "Significant in both",
        "Significant in SMC5-6 Complex",
        "Significant in instability",
        "Not Significant"
      ),
      name = "Significance"
    ) + # "grey70"
    # scale_alpha_manual(values = c(0.85, 0.5, 0.5, 0.1)) + # 0.1
    scale_alpha_manual(values = c(0.85, 0.4, 0.4, 0.08)) + # 0.1
    labs(
      color = "",
      y = "Instability log Fold-Change",
      x = "SMC5/6 complex log Fold-Change"
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
      # legend.position = "top",
      # legend.direction = "horizontal"
      legend.position = "right",
      legend.direction = "vertical",
      legend.text = element_text(size = 10)
    )
)
p <- rmbg(p)
pdf("./figures/figS5a.pdf", height = 6, width = 8)
plot(p)
dev.off()
# # if you want to save the raw results:
# colnames(mergeRes) <- gsub("\\.x$", "-SMC5-6", colnames(mergeRes))
# colnames(mergeRes) <- gsub("\\.y$", "-instability", colnames(mergeRes))
# mergeRes <- mergeRes[, -grep("^t-", colnames(mergeRes))]
# mergeRes <- mergeRes[, -grep("^logFDR", colnames(mergeRes))]
# mergeRes <- mergeRes[, -grep("^B-", colnames(mergeRes))]
# mergeRes <- mergeRes[, -grep("^Sig-", colnames(mergeRes))]
# dim(mergeRes[mergeRes$bothSig == "Both", ])
# table(mergeRes$bothSig)
# write_xlsx(mergeRes, "./reviewer-addressing/METABRIC-lfc-analysis.xlsx")
print("done")
