###################
## Configuration ##
###################
library(ggplot2)
library(survminer)
library(survival)
# columns used for analysis
complexGenes <- c("SMC5", "SMC6", "NSMCE1", "NSMCE2", "NSMCE3", "NSMCE4A", "EID3", "SLF1", "SLF2")
complexGenesDet <- paste0(complexGenes, "_det")
cols <- c("isalt", "NSMCE2")
# color of the plot
raw_df <- readRDS("./data/cbioportal/cbpdDataWInst.rds")
table(raw_df$NSMCE2, raw_df$MYC)
raw_df$OVT <- as.numeric(raw_df$OVT)
raw_df$OVS <- as.numeric(raw_df$OVS)
# ref_df <- ref_df[ref_df$major == "Breast", ]
# ref_df <- ref_df[ref_df$major != "Other", ]
# table(ref_df$NSMCE2, ref_df$MYC)
# ref_df[ref_df$NSMCE2 == 1 & ref_df$MYC == 0, ]

# for (i in cols) {
#     raw_df[, i] <- ifelse(raw_df[, i] == 1, "Altered", "Wild-type")
#     raw_df[, i] <- factor(raw_df[, i], levels = c("Altered", "Wild-type"))
# }


isalt_det <- apply(raw_df, 1, function(x) {
    v <- x[complexGenesDet]
    v <- gsub(" \\(putative passenger\\)", "", v)
    v <- tolower(v)
    v <- v[!duplicated(v)]
    v <- paste(v, collapse = ";")
    v <- gsub(";;", ";", v)
    v <- gsub(";;", ";", v)
    v <- gsub(";;", ";", v)
    v <- gsub(";;", ";", v)
    v <- gsub("^;", "", v)
    v <- gsub(";$", "", v)
    v <- gsub("amplification;amplification", "amplification", v)
    v <- gsub(" \\(putative passenger\\)", "", v)
})
# print(table(isalt_det))

raw_df$isalt_det <- isalt_det
# ref_df$isalt_det <- gsub("(putative passenger)", "", ref_df$isalt_det)


raw_df$isaltNew <- NA
raw_df$isaltNew[-grep("amplification", raw_df$isalt_det, ignore.case = TRUE)] <- "Mutation"
raw_df$isaltNew[raw_df$isalt_det == ""] <- "Wild-type"
raw_df$isaltNew[raw_df$isalt_det == "amplification"] <- "Amplification"
raw_df$isaltNew <- factor(raw_df$isaltNew, levels = c("Wild-type", "Amplification", "Mutation"))
# x <- ref_df[!is.na(ref_df$OVS), ]
# print(table(x$isaltNew, x$MYC))

# code taken from:
# https://stackoverflow.com/questions/3483203/create-a-boxplot-in-r-that-labels-a-box-with-the-sample-size-n
# rain cloud plots from
# https://github.com/RainCloudPlots/RainCloudPlots
# libs

KM_survival_plot <- function(sur_df, title, censor = TRUE) {
    colors <- c("#DE3B1C", "#707176")
    gc()
    # censoring function
    if (censor == TRUE) {
        sur_df <- sur_df[!is.na(sur_df$OVS) & !is.na(sur_df$OVT), ]
        # censorT <- 60
        # censorInt <- 10
        # censorT <- 240
        # censorInt <- 40
        # censorT <- 200
        # censorInt <- 40
        # censorT <- 180
        # censorInt <- 30
        #
        # censorT <- 250
        censorT <- 250
        # censorT <- 240
        censorInt <- 50
        sur_df$OVS[sur_df$OVT >= censorT] <- 0
        sur_df$OVT[sur_df$OVT >= censorT] <- censorT
    }

    diff <- survdiff(Surv(OVT, OVS) ~ group, data = sur_df)
    p.val <- 1 - pchisq(diff$chisq, length(diff$n) - 1)
    p.val_text <- paste("P = ", signif(p.val, 2))
    print(p.val)

    fit <- NA
    fit <- survfit(Surv(OVT, OVS) ~ group, data = sur_df)
    # https://stat.ethz.ch/pipermail/r-help/2007-April/130676.html
    # survival plot

    x <- unique(as.character(sur_df$group))
    x <- na.omit(x)
    x <- length(x)
    if (x == 2) {
        # legend.labs <- c("Wild-type", "Amplification") # , "Mutation")
        legend.labs <- c("Not amplified", "Amplified")
        colors <- c("#707176", "#DE3B1C") # , "#377eb8")
    }
    if (x == 3) {
        legend.labs <- c("Wild-type", "Amplification", "Mutation")

        colors <- c("#707176", "#DE3B1C", "#377eb8")
    }
    ggsurv <-
        ggsurvplot(fit,
            conf.int = TRUE,
            censor = FALSE,
            pval = FALSE,
            palette = colors,
            xlab = "",
            ylab = "Probability of overall survival",
            title = title,
            legend.title = "Status",
            legend = c(.85, .85),
            legend.labs = legend.labs,
            risk.table = TRUE,
            axes.offset = FALSE,
            risk.table.height = 0.22
        )
    ggsurv$table <- ggrisktable(fit,
        data = sur_df,
        ylab = "",
        xlab = "Time (Months)",
        risk.table.title = "",
        palette = colors,
        color = "strata",
        legend = "none",
        axes.offset = TRUE,
        tables.theme = theme_classic(),
        fontsize = 3.25,
        risk.table.col = "strata",
        # NEEDED TO GET THE RISK TABLE CORRECTLY
        # break.time.by = ifelse(censor, 10, 50)
        break.time.by = ifelse(censor, censorInt, 50)
    ) +
        scale_y_discrete(labels = rev(legend.labs)) +
        theme(
            axis.text.x = element_text(size = 12),
            axis.title.x = element_text(size = 12),
            axis.title.y = element_text(size = 12),
            legend.text = element_blank(),
            legend.title = element_blank()
        )
    ggsurv$plot <- ggsurv$plot + theme(plot.title = element_text(hjust = 0.5, size = 12))
    ggsurv$plot <- ggsurv$plot + theme(axis.text.x = element_blank())
    ggsurv$plot <- ggsurv$plot + theme(plot.margin = unit(c(5, 5, 0, 5), "points"))
    ggsurv$table <- ggsurv$table + theme(plot.margin = unit(c(5, 5, 0, 5), "points"))
    if (censor == TRUE) {
        ggsurv$plot <- ggsurv$plot + scale_x_continuous(limits = c(0, censorT), breaks = c(seq(0, censorT, censorInt), censorT))
        ggsurv$table <- ggsurv$table + scale_x_continuous(limits = c(0, censorT), breaks = c(seq(0, censorT, censorInt), censorT))
    }
    # else {
    #     ggsurv$table <- ggsurv$table + scale_x_continuous(limits = c(0, 360), breaks = seq(0, 350, 50))
    #     ggsurv$plot <- ggsurv$plot + scale_x_continuous(limits = c(0, 360), breaks = seq(0, 350, 50))
    # }
    if (signif(p.val, 2) == 1.7e-5) {
        print(p.val)
        p.val_text <- bquote("P = " ~ "1.7" ~ "x" ~ 10^"-5")
    }

    ggsurv$plot <- ggsurv$plot + annotate(
        geom = "text",
        x = ifelse(censor, censorInt, 50),
        y = 0.1,
        color = "black",
        size = 3,
        label = p.val_text
    )
    # this is to bring risk table up
    rm(fit, sur_df)
    fit <- NA
    return(ggsurv)
}





raw_df <- raw_df
# raw_df$MYC <- as.numeric(grepl("amp_rec", raw_df$MYC_det))

pdf("./reviewer-addressing/custom.pdf", width = 5, height = 5)

df <- raw_df
df <- df[df$major == "Breast", ]
# df <- df[grep("amp", df$MYC_det), ]
sur_df <- df
table(sur_df$NSMCE2_det)
sur_df$group <- ifelse(grepl("amp", sur_df$NSMCE2_det, ignore.case = TRUE), "Amplified", "Not amplified")
sur_df$group <- factor(sur_df$group, levels = c("Not amplified", "Amplified"))

# sur_df$group <- NA
# sur_df$group[-grep("amp", sur_df$NSMCE2_det, ignore.case = TRUE)] <- "Mutation"
# sur_df$group[sur_df$NSMCE2_det == ""] <- "Wild-type"
# sur_df$group[grep("amp", sur_df$NSMCE2_det, ignore.case = TRUE)] <- "Amplification"
# sur_df$group <- factor(sur_df$group, levels = c("Wild-type", "Amplification", "Mutation"))
# print(table(grepl("amp", sur_df$MYC_det), sur_df$group))
# sur_df <- sur_df[grep("amp", sur_df$MYC_det), ]
# print(KM_survival_plot(sur_df = sur_df, title = "NSMCE2 - MYC amplified - Breast", censor = TRUE))
# sur_df <- sur_df[sur_df$group != "Mutation", ]
print(KM_survival_plot(sur_df = sur_df, title = "NSMCE2 - MYC amplified - Breast", censor = TRUE))

dev.off()
print("done")
# df <- df[!is.na(df$OVT) & !is.na(df$OVS), ]

# print(table(grepl("amp", df$MYC_det), grepl("amp", df$NSMCE2_det, ignore.case = TRUE)))
# print(table(df$MYC, df$NSMCE2))
