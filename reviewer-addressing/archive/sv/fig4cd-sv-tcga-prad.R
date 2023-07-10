###################
## Configuration ##
###################
complexGenes <- c("SMC5", "SMC6", "NSMCE1", "NSMCE2", "NSMCE3", "NSMCE4A", "EID3", "SLF1", "SLF2")
complexGenesDet <- paste0(complexGenes, "_det")
studies <- c("TCGA")
# columns used for analysis
cols <- c(isalt_det = "Complex_det") # , NSMCE2_det == "NSMCE2_det")
# color of the plot
ref_df <- readRDS("./data/cbioportal/cbpdDataWInst.rds")
print(table(ref_df$major))
# ref_df <- ref_df[ref_df$major == "Breast", ]
ref_df <- ref_df[ref_df$major == "Prostate", ]

ref_df$isalt_det <- apply(ref_df, 1, function(x) {
    paste(x[complexGenesDet], collapse = ";")
})
table(grep("sv", ref_df$isalt_det, value = TRUE))
# ref_df <- ref_df[ref_df$MYC == 0, ]
# code taken from:
# https://stackoverflow.com/questions/3483203/create-a-boxplot-in-r-that-labels-a-box-with-the-sample-size-n
# rain cloud plots from
# https://github.com/RainCloudPlots/RainCloudPlots
# libs
library(ComplexHeatmap)
library(ggplot2)
library(survminer)
library(survival)

KM_survival_plot <- function(sur_df, title, xlab = TRUE, ylab = TRUE, strata = FALSE, censor = TRUE) {
    # colors <- c("#DE3B1C", "#707176")
    gc()
    # censoring function
    if (censor == TRUE) {
        sur_df$OVS[sur_df$OVT > 60] <- 0
        sur_df$OVT[sur_df$OVT > 60] <- 60
    }

    diff <- survdiff(Surv(OVT, OVS) ~ group, data = sur_df)
    p.val <- 1 - pchisq(diff$chisq, length(diff$n) - 1)
    p.val_text <- paste("P = ", signif(p.val, 2))

    title <- title
    fit <- NA
    fit <- survfit(Surv(OVT, OVS) ~ group, data = sur_df)
    print(fit)
    # https://stat.ethz.ch/pipermail/r-help/2007-April/130676.html
    # survival plot
    return(ggsurvplot(fit))
}



sur_dfs <- list()
titles <- list()
for (current_study in studies) {
    print(current_study)

    df <- ref_df[grep(current_study, ref_df$study), ]

    for (l in 1:length(cols)) {
        col <- names(cols[l])
        namecol <- as.character(cols[l])
        print(namecol)
        sur_df <- data.frame(
            study = df$study, tissue = df$tissue, OVT = df$OVT,
            OVS = df$OVS, col = df[, col]
        )
        sur_df$group <- "Wild-type"
        sur_df$group[grepl("Amplification", sur_df$col)] <- "Amplification"
        sur_df$group[grepl("sv", sur_df$col)] <- "SV"
        print(table(sur_df$group))
        # sur_df <- sur_df[sur_df$]
        sur_df$group <- as.factor(sur_df$group)
        sur_df$OVT <- as.numeric(sur_df$OVT)

        sur_df$OVS <- as.numeric(sur_df$OVS)
        #         sur_df$group <- ifelse(sur_df$group == 1, "Altered", "Wild") # namecol
        sur_dfs[[paste(current_study, namecol)]] <- na.omit(sur_df)
        titles[[paste(current_study, namecol)]] <- namecol
    }
}


pdf("./reviewer-addressing/fig4cd-sv-tcga-prad.pdf", width = 6, height = 5)
p <- list(
    (KM_survival_plot(sur_dfs[[1]], title = "SMC5/6 complex", censor = FALSE))
    # (KM_survival_plot(sur_dfs[[2]], title = "NSMCE2", censor = FALSE))
)
arrange_ggsurvplots(p,
    print = TRUE,
    ncol = 1,
    nrow = 1
)
dev.off()
print("done")
table(df$NSMCE2)
