###################
## Configuration ##
###################
library(ggplot2)
library(survminer)
library(survival)
# columns used for analysis
cols <- c("isalt", "NSMCE2")
# color of the plot
ref_df <- readRDS("./data/cbioportal/cbpdDataWInst.rds")
ref_df$OVT <- as.numeric(ref_df$OVT)
ref_df$OVS <- as.numeric(ref_df$OVS)
# ref_df <- ref_df[ref_df$major == "Breast", ]
# ref_df <- ref_df[ref_df$major != "Other", ]
table(ref_df$NSMCE2, ref_df$MYC)
ref_df[ref_df$NSMCE2 == 1 & ref_df$MYC == 0, ]

for (i in cols) {
    ref_df[, i] <- ifelse(ref_df[, i] == 1, "Altered", "Wild-type")
    ref_df[, i] <- factor(ref_df[, i], levels = c("Altered", "Wild-type"))
}

# code taken from:
# https://stackoverflow.com/questions/3483203/create-a-boxplot-in-r-that-labels-a-box-with-the-sample-size-n
# rain cloud plots from
# https://github.com/RainCloudPlots/RainCloudPlots
# libs

KM_survival_plot <- function(sur_df, title, censor = TRUE) {
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

    fit <- NA
    fit <- survfit(Surv(OVT, OVS) ~ group, data = sur_df)
    # https://stat.ethz.ch/pipermail/r-help/2007-April/130676.html
    # survival plot
    ggsurv <-
        ggsurvplot(fit,
            conf.int = TRUE,
            pval = FALSE,
            # palette = colors,
            xlab = "",
            ylab = "Probability of overall survival",
            title = title,
            legend.title = "Status",
            legend = c(.85, .4),
            # legend.labs = c("Altered", "Wild"),
            risk.table = TRUE,
            axes.offset = FALSE,
            risk.table.height = 0.22
        )
    ggsurv$table <- ggrisktable(fit,
        data = sur_df,
        ylab = "",
        xlab = "Time (Months)",
        risk.table.title = "",
        # palette = colors,
        color = "strata",
        legend = "none",
        axes.offset = TRUE,
        tables.theme = theme_classic(),
        fontsize = 3.25,
        risk.table.col = "strata"
        # NEEDED TO GET THE RISK TABLE CORRECTLY
        , break.time.by = ifelse(censor, 10, 50)
    ) +
        theme(
            axis.text.x = element_text(size = 12),
            axis.title.x = element_text(size = 12),
            axis.title.y = element_text(size = 12),
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
    if (signif(p.val, 1) == 1e-6) {
        print(p.val)
        p.val_text <- bquote("p = " ~ "1" ~ "x" ~ 10^-6)
    }

    ggsurv$plot <- ggsurv$plot + annotate(
        geom = "text", x = ifelse(censor, 10, 50), y = 0.1,
        color = "black", size = 3,
        label = p.val_text
    )
    # this is to bring risk table up
    rm(fit, sur_df)
    fit <- NA
    return(ggsurv)
}



sur_dfs <- list()
titles <- list()






raw_df <- ref_df
# raw_df$MYC <- as.numeric(grepl("amp_rec", raw_df$MYC_det))




pdf("./reviewer-addressing/fig4cd-mycMut-metabric.pdf", width = 5, height = 5)

tryCatch(expr = {
    df <- ref_df
    df <- df[df$MYC == 1, ]
    sur_df <- df
    sur_df$group <- sur_df$isalt
    print(KM_survival_plot(sur_df = sur_df, title = "SMC5/6 complex - MYC altered - all", censor = FALSE))
    print(KM_survival_plot(sur_df = sur_df, title = "SMC5/6 complex - MYC altered - all", censor = TRUE))
    sur_df$group <- sur_df$NSMCE2
    # print(KM_survival_plot(sur_df = sur_df, title = "NSMCE2 - MYC altered - all", censor = FALSE))
    # print(KM_survival_plot(sur_df = sur_df, title = "NSMCE2 - MYC altered - all", censor = TRUE))
}, error = function(cond) message(cond))


tryCatch(expr = {
    df <- ref_df
    df <- df[df$MYC == 0, ]
    sur_df <- df
    sur_df$group <- sur_df$isalt
    print(KM_survival_plot(sur_df = sur_df, title = "SMC5/6 complex - MYC WT - all", censor = FALSE))
    print(KM_survival_plot(sur_df = sur_df, title = "SMC5/6 complex - MYC WT - all", censor = TRUE))
    sur_df$group <- sur_df$NSMCE2
    # print(KM_survival_plot(sur_df = sur_df, title = "NSMCE2 - MYC WT - all", censor = FALSE))
    # print(KM_survival_plot(sur_df = sur_df, title = "NSMCE2 - MYC WT - all", censor = TRUE))
}, error = function(cond) message(cond))

for (i in unique(raw_df$major)) {
    tryCatch(expr = {
        df <- raw_df[raw_df$major == i, ]
        df <- df[df$MYC == 1, ]
        sur_df <- df
        sur_df$group <- sur_df$isalt
        print(KM_survival_plot(sur_df = sur_df, title = paste("SMC5/6 complex - MYC altered", i), censor = FALSE))
        # print(KM_survival_plot(sur_df = sur_df, title = paste("SMC5/6 complex - MYC altered", i), censor = TRUE))
    }, error = function(cond) message(cond))

    # tryCatch(expr = {
    #     df <- raw_df[raw_df$major == i, ]
    #     df <- df[df$MYC == 1, ]
    #     sur_df <- df
    #     sur_df$group <- sur_df$NSMCE2
    #     # print(KM_survival_plot(sur_df = sur_df, title = paste("NSMCE2 - MYC altered", i), censor = FALSE))
    #     # print(KM_survival_plot(sur_df = sur_df, title = paste("NSMCE2 - MYC altered", i), censor = TRUE))
    # }, error = function(cond) message(cond))

    tryCatch(expr = {
        df <- raw_df[raw_df$major == i, ]
        df <- df[df$MYC == 0, ]
        sur_df <- df
        sur_df$group <- sur_df$isalt
        print(KM_survival_plot(sur_df = sur_df, title = paste("SMC5/6 complex - MYC WT", i), censor = FALSE))
        # print(KM_survival_plot(sur_df = sur_df, title = paste("SMC5/6 complex - MYC WT", i), censor = TRUE))
    }, error = function(cond) message(cond))
    # tryCatch(expr = {
    #     df <- raw_df[raw_df$major == i, ]
    #     df <- df[df$MYC == 0, ]
    #     sur_df <- df
    #     sur_df$group <- sur_df$NSMCE2
    #     # print(ggsurvplot(survfit(Surv(OVT, OVS) ~ group, data = sur_df)))
    #     # sur_df$OVS[sur_df$OVT > 60] <- 0
    #     # sur_df$OVT[sur_df$OVT > 60] <- 60
    #     # print(ggsurvplot(survfit(Surv(OVT, OVS) ~ group, data = sur_df)))
    #     print(KM_survival_plot(sur_df = sur_df, title = paste("NSMCE2 - MYC WT", i), censor = FALSE))
    #     # print(KM_survival_plot(sur_df = sur_df, title = paste("NSMCE2 - MYC WT", i), censor = TRUE))
    # }, error = function(cond) message(cond))
}

dev.off()
print("done")
