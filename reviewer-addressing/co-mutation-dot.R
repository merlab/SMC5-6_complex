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
    df$shared[df$pval > 0.05] <- NA
    df$pval[df$pval > 0.05] <- NA
    df$shared <- as.numeric(df$shared) * 100
    df$log10pval <- -log10(as.numeric(df$pval))
    complexOrder <- sort(table(df$complex), decreasing = FALSE)
    # print(complexOrder)
    df$complex <- factor(df$complex, levels = complexGenes)
    instabilityOrder <- sort(table(df$instability), decreasing = TRUE)
    # print(instabilityOrder)
    df$instability <- factor(df$instability, levels = instabilityGenes)
    df$x <- as.numeric(df$complex)
    df$y <- as.numeric(df$instability)
    # p <- ggplot(df, aes(x = x, y = y, size = shared, color = log10pval)) +
    p <- ggplot(df, aes(x = y, y = x, size = shared, color = log10pval)) +
        geom_point() +
        scale_size(
            name = "% co-occurance",
            breaks = seq(0, 100, 20),
            limits = c(0, 105),
            # range = c(3, 8.5)
            range = c(2, 7)
        ) +
        scale_color_viridis(
            name = expression("-" ~ "log"[10] ~ (P)),
            # limits = c(-log10(0.05), -log10(1e-6)),
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
    plotdf <- data.frame()
    for (i in complexGenes) {
        for (j in instabilityGenes) {
            complexAlt <- as.numeric(subcbpd[, i])
            instabilityAlt <- as.numeric(subcbpd[, j])
            if (length(unique(instabilityAlt)) == 1 || length(unique(complexAlt)) == 1) {
                pval <- 1
                shared <- 0
                next()
            } else {
                pval <- fisher.test(complexAlt, instabilityAlt)$p.value
                shared <- sum(complexAlt == 1 & instabilityAlt == 1) / sum(complexAlt == 1)
            }
            plotdf <- rbind(plotdf, c(complex = i, instability = j, pval = pval, shared = shared))
        }
    }
    colnames(plotdf) <- c("complex", "instability", "pval", "shared")

    return(plotdf)
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
# pdf("temp.pdf", height = 4.5, width = 5.5, onefile = TRUE)
l <- list()
p <- list()
for (cancerType in cancerTypes) {
    print(cancerType)
    subcbpd <- cbpd[cbpd$major == cancerType, ]
    subcbpd <- subcbpd[, c(complexGenes, instabilityGenes)]
    plotdf <- generatePlotDf(subcbpd)
    p[[cancerType]] <- makeDotHeatmap(plotdf, cancerType, instabilityGenes, complexGenes)
}
# p[[2]] <- p[[2]] + theme(axis.text.y = element_blank())
# p[[3]] <- p[[3]] + theme(axis.text.y = element_blank())
# pdf("temp.pdf", height = 4, width = 6, onefile = TRUE)
# pdf("./reviewer-addressing/co-mutation.pdf", height = 8, width = 12, onefile = TRUE)
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
plot(makeDotHeatmap(plotdf, cancerType, instabilityGenes, complexGenes))
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
