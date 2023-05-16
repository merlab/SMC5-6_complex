# purpose: makes the dataframe for the alteration analysis of the complexes in the genes
source("./R/routine_tasks.R")
library(ggplot2)
library(viridis)
library(ggpubr)
# library(cooccur)
complexGenes <- c("NSMCE2", "SMC6", "SMC5", "NSMCE1", "NSMCE3", "NSMCE4A", "EID3")
instabilityGenes <- c(
    "TP53", "BRCA1", "BRCA2", "NBN", "TTK", "AURKA", "PLK1",
    "CHEK2", "CCNE1", "RB1", "RECQL4", "BLM"
)

makeDotHeatmap <- function(df, title = NA, instabilityGenes, complexGenes) {
    # df$pval <- as.numeric(df$pval)
    # df$pval <- p.adjust(df$pval, method = "fdr")
    # df$pval <- p.adjust(df$pval, method = "bonferroni")
    # print(df[df$complex == "NSMCE2" & df$instability == "TP53", ])
    # df$shared[df$pval > 0.05] <- NA
    # df$pval[df$pval > 0.05] <- NA
    # df$pval[df$pval < 1e-100] <- 1e-100
    # df$shared <- as.numeric(df$shared) * 100
    # df$log10pval <- -log10(as.numeric(df$pval))
    # df$fdr[df$fdr > 0.05] <- NA
    df$ppv <- df$ppv * 100
    df$fn[df$fn > 0.05] <- NA
    df$ppv[df$fn > 0.05] <- NA
    # df$ppv[df$fdr > 0.05] <- NA
    df$log10fdr <- -log10(df$fdr)
    # print(df[df$complex == "NSMCE2" & df$instability == "TP53", ])
    complexOrder <- sort(table(df$complex), decreasing = FALSE)
    # print(complexOrder)
    df$complex <- factor(df$complex, levels = complexGenes)
    instabilityOrder <- sort(table(df$instability), decreasing = TRUE)
    # print(instabilityOrder)
    df$instability <- factor(df$instability, levels = instabilityGenes)
    df$x <- as.numeric(df$complex)
    df$y <- as.numeric(df$instability)
    # p <- ggplot(df, aes(x = y, y = x, size = shared, color = log10pval)) +
    # p <- ggplot(df, aes(x = y, y = x, size = ppv, color = log10fdr)) +
    p <- ggplot(df, aes(x = y, y = x, size = ppv, color = fn)) +
        geom_point() +
        scale_size(
            # name = "% co-occurance",
            name = "Positive Predictive Value",
            breaks = seq(0, 100, 20),
            limits = c(0, 105),
            # range = c(3, 8.5)
            range = c(2, 7)
        ) +
        scale_color_viridis(
            # name = expression("-" ~ "log"[10] ~ (FDR)),
            name = "False Negative",
            limits = c(0, .05)
            # limits = c(-log10(0.05), -log10(1e-100)),
            # # breaks = c(-log10(0.05), -4, -8, -12, -14),
            # breaks = c(-log10(0.05), 10, 25, 50, 75, 100),
            # labels = c(0.05, 1e-10, 1e-25, 1e-50, 1e-75, 1e-100)
        ) +
        theme_bw() +
        coord_flip() +
        ggtitle(title) +
        scale_y_continuous(
            position = "right",
            expand = c(0, 0),
            breaks = seq(1, max(df$x)),
            limits = c(0.5, max(df$x) + 0.5),
            labels = levels(df$complex),
            minor_breaks = seq(0.5, max(df$x))
        ) +
        scale_x_continuous(
            position = "bottom",
            expand = c(0, 0),
            breaks = seq(1, max(df$y)),
            limits = c(0.5, max(df$y) + 0.5),
            labels = levels(df$instability),
            minor_breaks = seq(0.5, max(df$y))
        ) +
        xlab("") +
        ylab("")
    p <- rmbg(p)
    p <- p + theme(
        axis.text.x = element_text(
            angle = 90,
            hjust = -0.05,
            color = "black",
            size = 10
        ),
        axis.text.y = element_text(color = "black", size = 10),
        panel.grid.major = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.minor = element_line(color = "white", linewidth = 1),
        panel.background = element_rect(fill = "#E2E2E2", color = "#E2E2E2"),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
        # plot.title = element_text(size = 12)
    )
    return(p)
}


generatePlotDf <- function(subcbpd) {
    d <- data.frame()
    for (i in complexGenes) {
        for (j in instabilityGenes) {
            complexAlt <- as.numeric(subcbpd[, i])
            instabilityAlt <- as.numeric(subcbpd[, j])
            tp <- sum(complexAlt == 1 & instabilityAlt == 1)
            tn <- sum(complexAlt == 0 & instabilityAlt == 0)
            fn <- sum(complexAlt == 1 & instabilityAlt == 0)
            fp <- sum(complexAlt == 0 & instabilityAlt == 1)
            # plotdf <- rbind(plotdf, c(complex = i, instability = j, pval = pval, shared = shared))
            d <- rbind(d, c(complex = i, instability = j, tp = tp, tn = tn, fp = fp, fn = fn))
        }
    }
    colnames(d) <- c("complex", "instability", "tp", "tn", "fp", "fn")
    d$tp <- as.numeric(d$tp)
    d$tn <- as.numeric(d$tn)
    d$fp <- as.numeric(d$fp)
    d$fn <- as.numeric(d$fn)
    d[, 3:6] <- d[, 3:6] / rowSums(d[, 3:6])
    d$ppv <- d$tp / (d$tp + d$fp)
    d$fdr <- 1 - d$ppv

    return(d)
}

# cbpd <- readRDS("./data/cbioportal/format_exOther.rds")
cbpd <- readRDS("./reviewer-addressing/cbioportal/format_exOther.rds")
instabilityGenes <- names(sort(colSums(cbpd[, instabilityGenes])))

cancerTypes <- c("All", unique(levels(cbpd$major)))
cancerTypes <- cancerTypes[cancerTypes != "Other"]

cancerTypes <- c(
    "Breast",
    "Ovarian",
    "Prostate"
)
print(cancerTypes)
l <- list()
p <- list()
plotdfs <- data.frame()
for (cancerType in cancerTypes) {
    print(cancerType)
    subcbpd <- cbpd[cbpd$major == cancerType, ]
    subcbpd <- subcbpd[, c(complexGenes, instabilityGenes)]
    plotdf <- generatePlotDf(subcbpd)
    p[[cancerType]] <- makeDotHeatmap(plotdf, cancerType, instabilityGenes, complexGenes)
    plotdf$cancerType <- cancerType
    plotdfs <- rbind(plotdf, plotdfs)
}
pdf("./reviewer-addressing/plot/co-mutation.pdf", height = 5, width = 10.2, onefile = TRUE)
plot(
    ggarrange(
        plotlist = p, nrow = 1, ncol = 3,
        common.legend = TRUE,
        legend = "right",
        # labels = LETTERS[1:3]
        align = "hv"
    )
)

dev.off()

pdf("./reviewer-addressing/plot/co-mutation-all.pdf", height = 5, width = 5, onefile = TRUE)
plotdf <- generatePlotDf(cbpd)
plot(makeDotHeatmap(plotdf, "All data", instabilityGenes, complexGenes))
plotdf$cancerType <- "all"
plotdfs <- rbind(plotdf, plotdfs)
dev.off()

print("done")
# saveRDS(l, "./cooccur.rds")
# scale_x_continuous(
#     position = "bottom",
#     expand = c(0, 0),
#     breaks = seq(1, max(df$x)),
#     limits = c(0.5, max(df$x) + 0.5),
#     labels = levels(df$complex),
#     minor_breaks = seq(0.5, max(df$x))
# ) +
# scale_y_continuous(
#     position = "right",
#     expand = c(0, 0),
#     breaks = seq(1, max(df$y)),
#     limits = c(0.5, max(df$y) + 0.5),
#     labels = levels(df$instability),
#     minor_breaks = seq(0.5, max(df$y))
# ) +
#
# install.packages("Rediscover")
library(Rediscover)


range(plotdfs$fdr)
plotdfs[which(plotdfs$fdr == min(plotdfs$fdr)), ]
