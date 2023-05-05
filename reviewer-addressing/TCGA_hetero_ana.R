# purpose: download RNASeq of TCGA PRAD from TCGA biolinker
source("./R/routine_tasks.R")
source("./R/R_rainclouds.R")
complexGenes <- c("NSMCE2", "SMC6", "SMC5", "NSMCE1", "NSMCE3", "NSMCE4A", "EID3")
folder_check("./reviewer-addressing/tcga/")
library(TCGAbiolinks)
library(ggplot2)
library(dplyr)
library(viridis)
# saveRDS(dataPrep, "./data/TCGA_PRAD/STAR-counts.rds")
allTypes <- c(
    "BRCA", "PRAD", "LUAD", "LUSC", "SKCM", "OV", "ESCA",
    "LIHC", "UCEC", "PAAD", "BLCA"
)
for (type in allTypes) {
    f <- sprintf("./reviewer-addressing/tcga/%s-mut.rds", type)
    if (file.exists(f)) next()
    query <- GDCquery(
        project = paste0("TCGA-", type),
        data.category = "Simple Nucleotide Variation",
        access = "open",
        legacy = FALSE,
        data.type = "Masked Somatic Mutation",
        workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
    )
    GDCdownload(
        query = query
    )
    dataPrep <- GDCprepare(
        query = query,
        save = TRUE
    )
    colnames(dataPrep)
    saveRDS(dataPrep, f)
}
pdf("./reviewer-addressing/plot/hetero.pdf")
plotDf <- data.frame()
for (type in allTypes) {
    print(type)
    f <- sprintf("./reviewer-addressing/tcga/%s-mut.rds", type)
    mat <- readRDS(f)
    mat <- mat[mat$Hugo_Symbol %in% complexGenes, ]
    plot(hist(mat$t_alt_count / (mat$t_ref_count + mat$t_alt_count)))
    # plot(hist(mat$t_alt_count / mat$t_ref_count), main = type)
    t <- data.frame(type = type, variationFreq = mat$t_alt_count / (mat$t_ref_count + mat$t_alt_count), gene = mat$Hugo_Symbol)
    plotDf <- rbind(plotDf, t)
}
dev.off()
print("done")
# t_ref_count t_alt_count n_depth n_ref_count n

# range(brca$t_alt_count / (brca$t_ref_count + brca$t_alt_count))
# brca <- brca[brca$Hugo_Symbol %in% complexGenes, ]
#
# # hist(maf$t_alt_count / maf$t_ref_count)
# table(brca$Hugo_Symbol)
# colnames(brca)
# brca$Tumor_Sample_Barcode
plotDf$gene <- factor(plotDf$gene, levels = complexGenes)
# pdf("./reviewer-addressing/plot/hetero-sum.pdf", height = 6, width = 12)
# pdf("./reviewer-addressing/plot/hetero-sum.pdf", height = 6, width = 8)
pdf("./reviewer-addressing/plot/hetero-sum.pdf", height = 5, width = 8)
plot(
    ggplot(plotDf, aes(x = gene, y = variationFreq, fill = gene)) + # x = type
        # geom_boxplot(alpha = .75, outline.color = NA) + # , width = .7) +
        geom_flat_violin(
            # position = position_nudge(x = .25, y = 0),
            position = position_nudge(x = .35, y = 0),
            adjust = 2, trim = TRUE
        ) +
        geom_point(
            position = position_jitter(width = .15),
            size = .25
        ) +
        geom_boxplot(aes(x = as.numeric(gene) + 0.25),
            outlier.shape = NA, alpha = 0.3, width = .1, colour = "black",
        ) +
        geom_hline(yintercept = 0.5, color = "red", linetype = "dashed") +

        # https://datavizpyr.com/how-to-make-grouped-boxplot-with-jittered-data-points-in-ggplot2/
        # geom_point(position = position_jitterdodge(dodge.width = .05)) +
        # scale_color_viridis(discrete = TRUE) +
        scale_fill_viridis(discrete = TRUE, name = "Gene") +
        # xlab("Cancer Types") +
        xlab("Genes") +
        ylab("Allele Variation Frequency") +
        guides(fill = "none") +
        # scale_color_viridis(discrete = TRUE) +
        # scale_y_continuous(expand = c(0, 0), limits = c(0, 1), breaks = seq(0, .9, .2)) +
        # scale_y_continuous(expand = c(0, 0), limits = c(0, .5), breaks = seq(0, .5, .1)) +
        scale_y_continuous(expand = c(0, 0), limits = c(0, 1), breaks = seq(0, 1, .1)) +
        theme_classic() +
        theme(
            legend.position = "top",
            legend.direction = "horizontal"
        )
)
dev.off()
print(length(plotDf[plotDf$variationFreq >= 0.5, ]) / nrow(plotDf))
