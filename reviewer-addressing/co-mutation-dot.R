# purpose: makes the dataframe for the alteration analysis of the complexes in the genes
# library(Rediscover)
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


generatePlotDf <- function(subcbpd, cancerType, color, bgcolor) {
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

    geneList <- c(complexGenes, instabilityGenes)
    x <- data.frame(subcbpd[, geneList])
    x <- data.matrix(t(x))
    # PMA <- getPM(x)
    # mymutex <- getMutex(A=x, PM=PMA)
    # colnames(mymutex) <- geneList
    # rownames(mymutex) <- geneList
    # mymutex <- mymutex[instabilityGenes, complexGenes]
    # nr <- nrow(mymutex)
    # nc <- ncol(mymutex)
    # # mymutex <- matrix(p.adjust(mymutex, method = "fdr"), nrow = nr, ncol = nc)
    # mymutex <- matrix(p.adjust(mymutex, method = "bonferroni"), nrow = nr, ncol = nc)
    # colnames(mymutex) <- complexGenes
    # rownames(mymutex) <- instabilityGenes
    # print(mymutex)
    # mymutex[which(mymutex <= 0.01, arr.ind = TRUE)] <- 0
    # mymutex[which(mymutex > 0.01, arr.ind = TRUE)] <- 1
    # print(table(mymutex))
    # print(mymutex)

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
                geom_bar(stat = "identity", width = 1, color = "white", size= 0.05) +
                # scale_fill_viridis(
                scale_fill_manual(
                    # values = colors,
                    values = color,
                    name = " ",
                    labels = c("False negative", "True negative", "False positive", "True positive")
                ) +
                coord_polar("y", start = 0) +
                # guides(fill = "none") +
                theme_void() +
                theme(
                    plot.title = element_text(size = 1),
                    #plot.background = element_rect(fill = "#D3D3D380", color = "white"),
                    plot.background = element_blank(),
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
    ) +
                      theme(plot.background = element_rect(color = bgcolor, fill = bgcolor))
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
# colors1 <- c("#8DC73F80", "#7f7f7f80", "#ffa60080", "#FF4454")
# colors2 <- c("#8DC73F", "#006DA040", "#FF9027", "#FF4454")
# colors3 <- c("#8DC73F", "#00456640", "#ffa600", "#FF4454")
# colors4 <- c("#8DC73F80", "#00456633", "#ffa60080", "#FF4454")
# colors5 <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3")
# colors6 <- c("#fbb4ae", "#b3cde3", "#ccebc5", "#decbe4")
# colors7 <- c("#1b9e77", "#7570b3", "#d95f02", "#e7298a")
# colors8 <- c("#4daf4a", "#377eb880", "#984ea3", "#e41a1c")
bgcolor1 <- "#D5E4EB"
bgcolor2 <- "#F8F2E4"
bgcolor3 <- "white"
colors1 <- c("#1CD0BB","#DFDFDF","#00A9D7","#e41a1c")
colors2 <- c("#1CD0BB","#6E7C7C","#00A9D7","#e41a1c")
colors3 <- c("#00A9D7","#DFDFDF","#1CD0BB","#e41a1c")
colors4 <- c("#00A9D7","#6E7C7C","#1CD0BB","#e41a1c")

labels = c("False negative", "True negative", "False positive", "True positive")


colorsf <- c("#00A9D7", "#6E7C7C","#e41a1c","#1CD0BB")
# labels = c("False negative", "True negative", "False positive", "True positive")
for (cancerType in cancerTypes) {
    print(cancerType)
    subcbpd <- cbpd[cbpd$major == cancerType, ]
    subcbpd <- subcbpd[, c(complexGenes, instabilityGenes)]
    # plot(generatePlotDf(subcbpd, cancerType, color = colors1, bgcolor = bgcolor1))
    # plot(generatePlotDf(subcbpd, cancerType, color = colors2, bgcolor = bgcolor1))
    # plot(generatePlotDf(subcbpd, cancerType, color = colors3, bgcolor = bgcolor1))
    # plot(generatePlotDf(subcbpd, cancerType, color = colors4, bgcolor = bgcolor1))
    #
    # plot(generatePlotDf(subcbpd, cancerType, color = colors1, bgcolor = bgcolor2))
    # plot(generatePlotDf(subcbpd, cancerType, color = colors2, bgcolor = bgcolor2))
    # plot(generatePlotDf(subcbpd, cancerType, color = colors3, bgcolor = bgcolor2))
    # plot(generatePlotDf(subcbpd, cancerType, color = colors4, bgcolor = bgcolor2))
    #
    # plot(generatePlotDf(subcbpd, cancerType, color = colors1, bgcolor = bgcolor3))
    # plot(generatePlotDf(subcbpd, cancerType, color = colors2, bgcolor = bgcolor3))
    # plot(generatePlotDf(subcbpd, cancerType, color = colors3, bgcolor = bgcolor3))
    # plot(generatePlotDf(subcbpd, cancerType, color = colors4, bgcolor = bgcolor3))
    # plot(generatePlotDf(subcbpd, cancerType, color = colors4))
    # plot(generatePlotDf(subcbpd, cancerType, color = colors5))
    # plot(generatePlotDf(subcbpd, cancerType, color = colors6))
    # plot(generatePlotDf(subcbpd, cancerType, color = colors7))
    # plot(generatePlotDf(subcbpd, cancerType, color = colors8))
    plot(generatePlotDf(subcbpd, cancerType, color = colorsf, bgcolor = bgcolor2))
}
dev.off()

pdf("./reviewer-addressing/plot/co-mutation-all.pdf", height = 9, width = 15, onefile = TRUE)
plot(generatePlotDf(cbpd, "All", color = colorsf, bgcolor = bgcolor2))
dev.off()

print("done")
