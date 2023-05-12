# code taken from:
# https://stackoverflow.com/questions/3483203/create-a-boxplot-in-r-that-labels-a-box-with-the-sample-size-n
# rain cloud plots from
# https://github.com/RainCloudPlots/RainCloudPlots
library(ComplexHeatmap)
library(ggplot2)
library(ggpubr)
library(survminer)
library(survival)
library(dplyr)
library(cowplot)
source("./R/R_rainclouds.R")
dir <- getwd()
source("./R/routine_tasks.R")
genes <- c("NSMCE2", "Complex")
tissue <- "Prostate"
colors <- c("#DE3B1C", "#707176")
complexGenes <- c("SMC5", "SMC6", "NSMCE1", "NSMCE2", "NSMCE3", "NSMCE4A", "EID3", "SLF1", "SLF2")
studies <- c("Prostate Adenocarcinoma (TCGA, PanCancer Atlas)")
# o
# columns used for analysis
cols <- c(Signature = "SMC5/6 derived Expression Signature") # , NSMCE2 = "NSMCE2")
#
expmat <- readRDS("./data/tcga-prad/rnaseq.rds")
expmat <- apply(expmat, 1, function(x) {
    return((x - mean(x)) / sd(x))
})

signatureGenes <- readLines("./results/metabricSignature.txt")
signatureScore <- rowMeans(expmat[, colnames(expmat) %in% signatureGenes])
#
ref_df <- readRDS("./data/cbioportal/format_exOther.rds")
ref_df <- ref_df[rownames(expmat), ]
ref_df$Signature <- signatureScore
ref_df$Signature <- ref_df$Signature > median(ref_df$Signature)

KM_survival_plot <- function(sur_df, colors, title, xlab = TRUE, ylab = TRUE, strata = FALSE) {
    # censoring function
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
            legend = c(.35, .25),
            legend.labs = c("Altered", "Wild"),
            risk.table = TRUE,
            axes.offset = FALSE,
            risk.table.height = 0.22
        )
    ggsurv$table <- ggrisktable(fit,
        data = sur_df,
        ylab = "",
        xlab = ifelse(xlab, "Time (Months)", ""),
        risk.table.title = ""
        # add color
        , palette = colors,
        color = "strata",
        legend = "none",
        axes.offset = TRUE,
        tables.theme = theme_classic(),
        fontsize = 3.25,
        risk.table.col = "strata"
        # NEEDED TO GET THE RISK TABLE CORRECTLY
        , break.time.by = 40,
        xlim = c(0, 160)
    ) + scale_y_discrete(labels = c("Wild", "Altered")) +
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
    ggsurv$plot <- ggsurv$plot + scale_x_continuous(limits = c(0, 170), breaks = seq(0, 160, 40)) #
    ggsurv$table <- ggsurv$table + scale_x_continuous(limits = c(0, 170), breaks = seq(0, 160, 40)) #

    ggsurv$plot <- ggsurv$plot + annotate(
        geom = "text", x = 30, y = 0.1,
        color = "black", size = 5,
        label = paste0("P = ", signif(p.val, 2))
    )
    # this is to bring risk table up
    rm(fit, sur_df)
    fit <- NA
    return(ggsurv)
}

sur_dfs <- list()
titles <- list()
for (current_study in studies) {
    print(current_study)

    df <- ref_df[ref_df$study == current_study, ]
    for (l in 1:length(cols)) {
        col <- names(cols[l])
        namecol <- as.character(cols[l])
        print(namecol)
        sur_df <- data.frame(
            study = df$study, tissue = df$tissue, OVT = df$OVT,
            OVS = df$OVS, group = df$Signature
        )

        sur_df$OVT <- as.numeric(sur_df$OVT)
        sur_df$OVS <- as.numeric(sur_df$OVS)
        sur_df$group <- as.numeric(sur_df$group)
        sur_df$group <- ifelse(sur_df$group == 1, "Altered", "Wild") # namecol
        # The top row are for breast cancer, the bottom row is for prostate cancer. #
        sur_dfs[[paste(current_study, namecol)]] <- na.omit(sur_df)
        titles[[paste(current_study, namecol)]] <- namecol
    }
}

pdf("./reviewer-addressing/plot/lfcExprSig-figS4cd.pdf", width = 5, height = 5)
print(KM_survival_plot(sur_dfs[[1]], color = colors, title = titles[[1]]))
dev.off()
print("done")
