# purpose: makes the dataframe for the alteration analysis of the complexes in the genes
source("./R/routine_tasks.R")
library(ggplot2)
library(gridExtra)
library(viridis)
library(ggpubr)
library(tidyr)
# library(cooccur)
complexGenes <- c("NSMCE2", "SMC6", "SMC5", "NSMCE1", "NSMCE3", "NSMCE4A", "EID3")
complexGenes <- c("SMC5", "SMC6", "NSMCE1", "NSMCE2", "NSMCE3", "NSMCE4A", "EID3")
instabilityGenes <- c(
    "TP53", "BRCA1", "BRCA2", "NBN", "TTK", "AURKA", "PLK1",
    "CHEK2", "CCNE1", "RB1", "RECQL4", "BLM"
)


labels <- c("False negative", "True negative", "False positive", "True positive")
colors <- c("#8DC73F", "#006DA040", "#FF9027", "#FF4454")


generatePlotDf <- function(subcbpd, cancerType) {
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
            v <- c(complex = i, instability = j, tp = tp, tn = tn, fp = fp, fn = fn)
            d <- rbind(d, v)
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

    x <- pivot_longer(d, c("tp", "tn", "fp", "fn"))
    p <- list()
    n <- 0
    n <- n + 1
    p[[n]] <- NULL
    for (i in instabilityGenes) {
        n <- n + 1
        # p[[n]] <- ggplot() +
        #     annotate(geom = "text", lable = i, size = 20, x = 1, y = 1) +
        #     theme_void()
        p[[n]] <- ggplot(data.frame(x = 1, y = 1), aes(x, y)) +
            geom_text(label = i, size = 5) +
            theme_void()
    }
    for (i in complexGenes) {
        n <- n + 1
        p[[n]] <- ggplot(data.frame(x = 1, y = 1), aes(x, y)) +
            geom_text(label = i, size = 5) +
            theme_void()
        for (j in instabilityGenes) {
            n <- n + 1
            subx <- x[x$complex == i & x$instability == j, ]
            subx$name <- factor(subx$name, levels = c("fn", "tn", "fp", "tp"))
            # print(head(x))
            p[[n]] <- ggplot(subx, aes(x = "", y = value, fill = name)) +
                geom_bar(stat = "identity", width = 1, color = NA) +
                # scale_fill_viridis(
                scale_fill_manual(
                    values = colors,
                    name = " ",
                    labels = c("False negative", "True negative", "False positive", "True positive")
                ) +
                coord_polar("y", start = 0) +
                # guides(fill = "none") +
                theme_void() +
                theme(
                    plot.title = element_text(size = 1),
                    plot.background = element_rect(fill = "#D3D3D380", color = "white"),
                    plot.margin = unit(c(1, 1, 1, 1), "pt")
                )
        }
    }

    PLOT <- ggarrange(
        plotlist = p,
        nrow = length(complexGenes) + 1,
        ncol = length(instabilityGenes) + 1,
        common.legend = TRUE,
        align = "hv"
    )
    annotate_figure(PLOT,
        top = text_grob(cancerType,
            face = "bold",
            size = 12
        )
    )
    return(PLOT)
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
pdf("./reviewer-addressing/plot/co-mutation.pdf", height = 9, width = 13, onefile = TRUE)
for (cancerType in cancerTypes[1]) {
    print(cancerType)
    subcbpd <- cbpd[cbpd$major == cancerType, ]
    subcbpd <- subcbpd[, c(complexGenes, instabilityGenes)]
    plot(generatePlotDf(subcbpd, cancerType))
}
dev.off()

# pdf("./reviewer-addressing/plot/co-mutation-all.pdf", height = 9, width = 13, onefile = TRUE)
pdf("./reviewer-addressing/plot/co-mutation-all.pdf", height = 9, width = 15, onefile = TRUE)
plot(generatePlotDf(cbpd, "All"))
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
# library(Rediscover)
# makeDotHeatmap <- function(df, title = NA, instabilityGenes, complexGenes) {
#     # df$pval <- as.numeric(df$pval)
#     # df$pval <- p.adjust(df$pval, method = "fdr")
#     # df$pval <- p.adjust(df$pval, method = "bonferroni")
#     # print(df[df$complex == "NSMCE2" & df$instability == "TP53", ])
#     # df$shared[df$pval > 0.05] <- NA
#     # df$pval[df$pval > 0.05] <- NA
#     # df$pval[df$pval < 1e-100] <- 1e-100
#     # df$shared <- as.numeric(df$shared) * 100
#     # df$log10pval <- -log10(as.numeric(df$pval))
#     # df$fdr[df$fdr > 0.05] <- NA
#     df$ppv <- df$ppv * 100
#     # df$fn[df$fn > 0.05] <- NA
#     # df$ppv[df$fn > 0.05] <- NA
#     # df$ppv[df$fdr > 0.05] <- NA
#     df$log10fdr <- -log10(df$fdr)
#     # print(df[df$complex == "NSMCE2" & df$instability == "TP53", ])
#     complexOrder <- sort(table(df$complex), decreasing = FALSE)
#     # print(complexOrder)
#     df$complex <- factor(df$complex, levels = complexGenes)
#     instabilityOrder <- sort(table(df$instability), decreasing = TRUE)
#     # print(instabilityOrder)
#     df$instability <- factor(df$instability, levels = instabilityGenes)
#     df$x <- as.numeric(df$complex)
#     df$y <- as.numeric(df$instability)
#     # p <- ggplot(df, aes(x = y, y = x, size = shared, color = log10pval)) +
#     # p <- ggplot(df, aes(x = y, y = x, size = ppv, color = log10fdr)) +
#     p <- ggplot(df, aes(x = y, y = x, size = ppv, color = fn)) +
#         geom_point() +
#         scale_size(
#             # name = "% co-occurance",
#             name = "Positive Predictive Value",
#             breaks = seq(0, 100, 20),
#             limits = c(0, 105),
#             # range = c(3, 8.5)
#             range = c(2, 7)
#         ) +
#         scale_color_viridis(
#             # name = expression("-" ~ "log"[10] ~ (FDR)),
#             name = "False Negative",
#             limits = c(0, .05)
#             # limits = c(-log10(0.05), -log10(1e-100)),
#             # # breaks = c(-log10(0.05), -4, -8, -12, -14),
#             # breaks = c(-log10(0.05), 10, 25, 50, 75, 100),
#             # labels = c(0.05, 1e-10, 1e-25, 1e-50, 1e-75, 1e-100)
#         ) +
#         theme_bw() +
#         coord_flip() +
#         ggtitle(title) +
#         scale_y_continuous(
#             position = "right",
#             expand = c(0, 0),
#             breaks = seq(1, max(df$x)),
#             limits = c(0.5, max(df$x) + 0.5),
#             labels = levels(df$complex),
#             minor_breaks = seq(0.5, max(df$x))
#         ) +
#         scale_x_continuous(
#             position = "bottom",
#             expand = c(0, 0),
#             breaks = seq(1, max(df$y)),
#             limits = c(0.5, max(df$y) + 0.5),
#             labels = levels(df$instability),
#             minor_breaks = seq(0.5, max(df$y))
#         ) +
#         xlab("") +
#         ylab("")
#     p <- rmbg(p)
#     p <- p + theme(
#         axis.text.x = element_text(
#             angle = 90,
#             hjust = -0.05,
#             color = "black",
#             size = 10
#         ),
#         axis.text.y = element_text(color = "black", size = 10),
#         panel.grid.major = element_blank(),
#         axis.ticks = element_blank(),
#         panel.grid.minor = element_line(color = "white", linewidth = 1),
#         panel.background = element_rect(fill = "#E2E2E2", color = "#E2E2E2"),
#         plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
#         # plot.title = element_text(size = 12)
#     )
#     return(p)
# }
#
