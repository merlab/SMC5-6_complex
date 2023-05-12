###################
## Configuration ##
###################
library(readxl)
current_study <- c("Breast Cancer (METABRIC, Nature 2012 & Nat Commun 2016)")
complexGenes <- c("SMC5", "SMC6", "NSMCE1", "NSMCE2", "NSMCE3", "NSMCE4A", "EID3", "SLF1", "SLF2")
# columns used for analysis
cols <- c(Signature = "Expression Signature") # , NSMCE2 = "NSMCE2")
#
expmat <- readRDS("./data/metabric-brca/microarray-metagx.rds")
# normalize the experssion matrix
expmat <- apply(expmat, 1, function(x) {
    return((x - mean(x)) / sd(x))
})

signatureGenes <- readLines("./results/metabricSignature.txt")
signatureScore <- rowMeans(expmat[, signatureGenes])

#

#
# ref_df <- readRDS("./data/cbioportal/format_exOther.rds")
# ref_df <- readRDS("./data/cbioportal/format_exOther.rds")
ref_df <- readRDS("./reviewer-addressing/cbioportal/formatted.rds")
ref_df <- ref_df[rownames(expmat), ]
ref_df$Signature <- signatureScore
ref_df$SignatureVal <- signatureScore
ref_df$Signature <- ref_df$Signature > median(ref_df$Signature)

# code taken from:
# https://stackoverflow.com/questions/3483203/create-a-boxplot-in-r-that-labels-a-box-with-the-sample-size-n
# rain cloud plots from
# https://github.com/RainCloudPlots/RainCloudPlots
# libs
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
    return(ggsurv)
}


sur_dfs <- list()
titles <- list()
print(current_study)

df <- ref_df[ref_df$study == current_study, ]
df$OVT <- as.numeric(df$OVT)
df$OVS <- as.numeric(df$OVS)
for (l in seq_along(cols)) {
    col <- names(cols[l])
    namecol <- as.character(cols[l])
    print(namecol)
    sur_df <- data.frame(
        study = df$study, tissue = df$tissue, OVT = df$OVT,
        OVS = df$OVS, group = df[, col]
    )
    sur_df$group <- as.numeric(sur_df$group)
    sur_df$group <- ifelse(sur_df$group == 1, "High expr", "Low expr") # namecol
    sur_dfs[[paste(current_study, namecol)]] <- na.omit(sur_df)
    titles[[paste(current_study, namecol)]] <- namecol
}


# pdf("./reviewer-addressing/plot/fig3cd.pdf", width = 10, height = 5)
pdf("./reviewer-addressing/plot/expmat-Fig3cd.pdf", width = 5, height = 5)
print(KM_survival_plot(sur_dfs[[1]], tit = "Expression Signature"))
dev.off()
print("done")
cox <- coxph(Surv(OVT, OVS) ~ SignatureVal + isalt + NSMCE2 + grade, data = df)
print(cox)
