###################
## Configuration ##
###################
source("./R/R_rainclouds.R")
source("./R/routine_tasks.R")
library(readxl)
library(cowplot)
library(ggpubr)
library(ggplot2)
rainCloudPlot <- function(ref_df) {
    colors <- c("#DE3B1C", "#707176")
    p <- ggplot(ref_df, aes(x = isalt, y = SignatureVal, color = isalt, fill = isalt)) +
        geom_flat_violin(
            position = position_nudge(x = .25, y = 0),
            adjust = 2, trim = TRUE
        ) +
        geom_point(position = position_jitter(width = .15), size = .25) +
        geom_boxplot(aes(x = as.numeric(isalt) + 0.25),
            outlier.shape = NA, alpha = 0.3, width = .1, colour = "BLACK"
        ) +
        ylab("Expression Signature") +
        xlab("") +
        scale_colour_manual(values = colors) +
        scale_fill_manual(values = colors) +
        theme_cowplot() +
        # wilcox t test raw output. comment to add custom p value
        stat_compare_means(method = "t.test") +
        theme(plot.title = element_text(hjust = 0.5, size = 12, face = "plain"))
    p <- rmbg(p)
}
current_study <- c("Breast Cancer (METABRIC, Nature 2012 & Nat Commun 2016)")
complexGenes <- c("SMC5", "SMC6", "NSMCE1", "NSMCE2", "NSMCE3", "NSMCE4A", "EID3")
# columns used for analysis0
cols <- c(Signature = "SMC5/6 derived Expression Signature") # , NSMCE2 = "NSMCE2")
#
expmat <- readRDS("./data/metabric-brca/microarray-metagx.rds")
# normalize the experssion matrix
expmat <- apply(expmat, 1, function(x) {
    return((x - mean(x)) / sd(x))
})

signatureGenes <- readLines("./results/metabricSignature.txt")
signatureScore <- rowMeans(expmat[, colnames(expmat) %in% signatureGenes])
# signatureScore <- rowMeans(expmat[, colnames(expmat) %in% complexGenes])
ref_df <- readRDS("./data/cbioportal/formatted.rds")
ref_df <- ref_df[rownames(expmat), ]
ref_df$SignatureVal <- signatureScore
ref_df <- ref_df[!is.na(ref_df$isalt), ]
ref_df$isalt <- as.factor(ifelse(ref_df$isalt == 1, "Altered", "Wild-type"))


pdf("./reviewer-addressing/plot/lfcExprSig-rainCloud.pdf",
    width = 5, height = 5
)
plot(rainCloudPlot(ref_df))
dev.off()


signatureGenes <- readLines("./results/metabricSignature.txt")
# signatureScore <- rowMeans(expmat[, colnames(expmat) %in% signatureGenes])
signatureScore <- rowMeans(expmat[, colnames(expmat) %in% complexGenes])
ref_df <- readRDS("./data/cbioportal/formatted.rds")
ref_df <- ref_df[rownames(expmat), ]
ref_df$SignatureVal <- signatureScore
ref_df <- ref_df[!is.na(ref_df$isalt), ]
ref_df$isalt <- as.factor(ifelse(ref_df$isalt == 1, "Altered", "Wild-type"))

pdf("./reviewer-addressing/plot/smcExprSig-rainCloud.pdf",
    width = 5, height = 5
)
plot(p)
dev.off()
print("done")
