###################
## Configuration ##
###################
studies <- c("Breast Cancer (METABRIC, Nature 2012 & Nat Commun 2016)")
complexGenes <- c("SMC5", "SMC6", "NSMCE1", "NSMCE2", "NSMCE3", "NSMCE4A", "EID3")
# columns used for analysis
# cols <- c(Complex = "Complex", NSMCE2 = "NSMCE2")
## NOTE: all NSMCE2 is alteration in metabric dataset
cols <- c(Complex = "Complex") # , NSMCE2 = "NSMCE2")
# color of the plot
ref_df <- readRDS("./data/cbioportal/formatted.rds")
ref_df$NSMCE2 <- ref_df$NSMCE2_det
ref_df$Complex <- apply(ref_df, 1, function(x) {
    return(paste(x[paste0(complexGenes, "_det")], collapse = ";"))
})
for (i in 1:10) {
    ref_df$Complex <- gsub(";;", ";", ref_df$Complex)
    ref_df$Complex <- gsub("NA", "", ref_df$Complex)
    ref_df$Complex <- gsub("^;", "", ref_df$Complex)
    ref_df$Complex[ref_df$Complex %in% c(";;", ";")] <- ""
}
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

KM_survival_plot <- function(sur_df, tit = "", xlab = TRUE, ylab = TRUE, strata = FALSE, censor = TRUE) {
    colors <- c("#DE3B1C", "#707176")
    gc()
    # censoring function
    if (censor == TRUE) {
        sur_df$OVS[sur_df$OVT > 60] <- 0
        sur_df$OVT[sur_df$OVT > 60] <- 60
    }
    diff <- survdiff(Surv(OVT, OVS) ~ group, data = sur_df)
    p.val <- 1 - pchisq(diff$chisq, length(diff$n) - 1)
    print(p.val)

    fit <- NA
    fit <- survfit(Surv(OVT, OVS) ~ group, data = sur_df)
    # https://stat.ethz.ch/pipermail/r-help/2007-April/130676.html
    # survival plot
    ggsurv <- ggsurvplot(fit,
        conf.int = TRUE,
        pval = FALSE,
        palette = colors,
        xlab = "",
        ylab = ifelse(ylab, "Probability of overall survival", ""),
        title = tit,
        legend.title = "Status",
        legend = c(.85, .4),
        legend.labs = c("Amplified", "Wild"),
        risk.table = TRUE,
        axes.offset = FALSE,
        risk.table.height = 0.22
    )
    ggsurv$plot <- ggsurv$plot + annotate(
        geom = "text", x = 60, y = 0.1,
        color = "black", size = 5,
        label = paste0("P = ", signif(p.val, 2))
    )
    ggsurv$table <- ggrisktable(fit,
        data = sur_df,
        ylab = "",
        xlab = ifelse(xlab, "Time (Months)", ""),
        risk.table.title = "",
        palette = colors,
        color = "strata",
        legend = "none",
        axes.offset = TRUE,
        tables.theme = theme_classic(),
        fontsize = 3.25,
        risk.table.col = "strata"
        # NEEDED TO GET THE RISK TABLE CORRECTLY
        , break.time.by = ifelse(censor, 10, 50)
    ) + scale_y_discrete(labels = c("Wild", "Altered")) +
        theme(
            axis.text.x = element_text(size = 12),
            axis.title.x = element_text(size = 12, face = "bold"),
            axis.title.y = element_text(size = 12, face = "bold"),
            legend.text = element_blank(),
            legend.title = element_blank()
        )
    ggsurv$plot <- ggsurv$plot + theme(plot.title = element_text(hjust = 0.5, size = 12))
    ggsurv$plot <- ggsurv$plot + theme(axis.text.x = element_blank())
    ggsurv$plot <- ggsurv$plot + theme(plot.title = element_text(hjust = 0.5, size = 12))
    ggsurv$plot <- ggsurv$plot + theme(plot.margin = unit(c(5, 5, 0, 5), "points"))
    ggsurv$table <- ggsurv$table + theme(plot.margin = unit(c(5, 5, 0, 5), "points"))
    if (censor == TRUE) {
        ggsurv$plot <- ggsurv$plot + scale_x_continuous(limits = c(0, 65), breaks = seq(0, 60, 10))
        ggsurv$table <- ggsurv$table + scale_x_continuous(limits = c(0, 65), breaks = seq(0, 60, 10))
    } else {
        ggsurv$table <- ggsurv$table + scale_x_continuous(limits = c(0, 360), breaks = seq(0, 350, 50))
        ggsurv$plot <- ggsurv$plot + scale_x_continuous(limits = c(0, 360), breaks = seq(0, 350, 50))
    }
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
            OVS = df$OVS
        )
        sur_df$group <- rep(1, nrow(df))
        print(table(sur_df$group))
        sur_df$group[df[, col] == ""] <- 0
        print(table(sur_df$group))
        sur_df$group[grep("Amplification", df[, col])] <- NA
        print(table(sur_df$group))
        sur_df <- na.omit(sur_df)
        sur_df$OVT <- as.numeric(sur_df$OVT)

        sur_df$OVS <- as.numeric(sur_df$OVS)
        sur_df$group <- as.numeric(sur_df$group)
        sur_df$group <- ifelse(sur_df$group == 1, "Altered", "Wild") # namecol
        sur_dfs[[paste(current_study, namecol)]] <- na.omit(sur_df)
        titles[[paste(current_study, namecol)]] <- namecol
    }
}


# pdf("./reviewer-addressing/plot/fig3cd.pdf", width = 10, height = 5)
pdf("./reviewer-addressing/plot/noAmp-Fig3cd.pdf", width = 5, height = 5)
print(KM_survival_plot(sur_dfs[[1]], tit = "Complex", censor = FALSE))
dev.off()
print("done")
