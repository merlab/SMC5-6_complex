###################
## Configuration ##
###################
# code taken from:
# https://stackoverflow.com/questions/3483203/create-a-boxplot-in-r-that-labels-a-box-with-the-sample-size-n
# rain cloud plots from
# https://github.com/RainCloudPlots/RainCloudPlots
# libs
library(readxl)
library(ComplexHeatmap)
library(ggplot2)
library(survminer)
library(survival)
library(ggpubr)
library(dplyr)

KM_survival_plot <- function(sur_df, title, xlab = TRUE, ylab = TRUE, strata = FALSE) {
    colors <- c("#DE3B1C", "#707176")
    gc()
    diff <- survdiff(Surv(OVT, OVS) ~ group, data = sur_df)
    p.val <- 1 - pchisq(diff$chisq, length(diff$n) - 1)
    print(p.val)
    print(max(sur_df$OVT))

    title <- title
    fit <- NA
    fit <- survfit(Surv(OVT, OVS) ~ group, data = sur_df)
    # https://stat.ethz.ch/pipermail/r-help/2007-April/130676.html
    # survival plot
    ggsurv <-
        ggsurvplot(fit,
            conf.int = TRUE,
            pval = FALSE,
            palette = colors,
            xlab = "",
            ylab = ifelse(ylab, "Probability of overall survival", ""),
            title = title,
            legend.title = "Status",
            legend = c(.85, .85),
            legend.labs = c("High expr", "Low expr"),
            risk.table = TRUE,
            axes.offset = FALSE,
            risk.table.height = 0.22
        )
    ggsurv$table <- ggrisktable(fit,
        data = sur_df,
        ylab = "",
        xlab = ifelse(xlab, "Time (Months)", ""),
        risk.table.title = "",
        # add color
        palette = rev(colors),
        color = "strata",
        legend = "none",
        axes.offset = TRUE,
        tables.theme = theme_classic(),
        fontsize = 3.25,
        risk.table.col = "strata",
        break.time.by = 50,
    ) + scale_y_discrete(labels = c("High expr", "Low expr")) +
        theme(
            axis.text.x = element_text(size = 12),
            axis.title.x = element_text(size = 12),
            legend.text = element_blank(),
            legend.title = element_blank()
        )
    ggsurv$plot <- ggsurv$plot + theme(plot.title = element_text(hjust = 0.5, size = 12))
    ggsurv$plot <- ggsurv$plot + theme(axis.text.x = element_blank())
    ggsurv$plot <- ggsurv$plot + theme(plot.title = element_text(hjust = 0.5, size = 12))
    ggsurv$plot <- ggsurv$plot + theme(plot.margin = unit(c(5, 5, 0, 5), "points"))
    ggsurv$table <- ggsurv$table + theme(plot.margin = unit(c(5, 5, 0, 5), "points"))

    ggsurv$plot <- ggsurv$plot + annotate(
        geom = "text", x = 80, y = 0.1,
        color = "black", size = 5,
        label = paste0("P = ", signif(p.val, 2))
    )
    # this is to bring risk table up
    rm(fit, sur_df)
    fit <- NA
    # return(ggsurv)
    print(ggsurv)
}
current_study <- c("Breast Cancer (METABRIC, Nature 2012 & Nat Commun 2016)")
complexGenes <- c("SMC5", "SMC6", "NSMCE1", "NSMCE2", "NSMCE3", "NSMCE4A", "EID3")
# columns used for analysis0
cols <- c(Signature = "SMC5/6 derived Expression Signature") # , NSMCE2 = "NSMCE2")

dgea <- readRDS("./data/metabric-brca/dgea-limma.rds")
dgea <- dgea[order(abs(dgea$logFC), decreasing = TRUE), ]
geneList <- dgea$gene
#
expmat <- readRDS("./data/metabric-brca/microarray-metagx.rds")
# normalize the experssion matrix
expmat <- apply(expmat, 1, function(x) {
    return((x - mean(x)) / sd(x))
})

ref_df <- readRDS("./data/cbioportal/formatted.rds")
ref_df <- ref_df[rownames(expmat), ]
ref_df$OVT <- as.numeric(ref_df$OVT)
ref_df$OVS <- as.numeric(ref_df$OVS)
ref_df <- ref_df[ref_df$study == current_study, ]




sur_dfs <- list()
titles <- list()
# signatureGeneList <- list(
#     "top10" = geneList[seq_len(10)],
#     "top20" = geneList[seq_len(20)],
#     "top25" = geneList[seq_len(25)],
#     "top50" = geneList[seq_len(50)],
#     "top100" = geneList[seq_len(100)],
#     "top200" = geneList[seq_len(200)],
#     "top250" = geneList[seq_len(250)],
#     "top500" = geneList[seq_len(250)],
#     "top1000" = geneList[seq_len(1000)],
#     "top2000" = geneList[seq_len(2000)],
#     "fdr0.01" = geneList[dgea$adj.P.Val <= 0.01],
#     "fdr0.05" = geneList[dgea$adj.P.Val <= 0.05]
# )
pdf("./reviewer-addressing/plot/lfcExprSig-Fig3cd.pdf", width = 5, height = 5)
# for (sigGeneName in names(signatureGeneList)) {
#     print(sigGeneName)
sur_df <- data.frame(
    # study = ref_df$study,
    # tissue = ref_df$tissue,
    OVT = ref_df$OVT,
    OVS = ref_df$OVS
)
signatureScore <- rowMeans(expmat[
    ,
    # colnames(expmat) %in% signatureGeneList[[sigGeneName]]
    colnames(expmat) %in% geneList[seq_len(10)]
])

sur_df$group <- signatureScore
sur_df <- na.omit(sur_df)
# hist(sur_df$group)
print(median(sur_df$group))
sur_df$group <- as.numeric(sur_df$group >= median(sur_df$group))
sur_df$group <- ifelse(sur_df$group == 1, "High expr", "Low expr") # namecol
KM_survival_plot(sur_df, title = "")
# KM_survival_plot(sur_df, title = sigGeneName)
# sur_dfs[[paste(current_study, namecol)]] <- na.omit(sur_df)
# titles[[paste(current_study, namecol)]] <- namecol
# }


dev.off()
print("done")
