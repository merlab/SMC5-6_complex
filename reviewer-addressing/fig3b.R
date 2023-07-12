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
# instabilityGenes <- c(
#     "TP53", "BRCA1", "BRCA2", "NBN", "TTK", "AURKA", "PLK1",
#     "CHEK2", "CCNE1", "RB1", "RECQL4", "BLM", "MYC"
# )

instabilityGenes <- c(
    "PLK1", "TTK", "CHEK2", "BLM", "AURKA", "BRCA1", "CCNE1", "NBN",
    "RECQL4", "BRCA2", "RB1", "TP53", "MYC"
)
bgcolor <- "#F8F2E4"
colors <- c("#00A9D7", "#6E7C7C", "#e41a1c", "#1CD0BB")
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
            v <- c(complex = i, instability = j, tp = tp, tn = tn, fp = fp, fn = fn)
            d <- rbind(d, v)
        }
    }

    geneList <- c(complexGenes, instabilityGenes)
    x <- data.frame(subcbpd[, geneList])
    x <- data.matrix(t(x))

    colnames(d) <- c("complex", "instability", "tp", "tn", "fp", "fn")
    d$tp <- as.numeric(d$tp)
    d$tn <- as.numeric(d$tn)
    d$fp <- as.numeric(d$fp)
    d$fn <- as.numeric(d$fn)
    d[, 3:6] <- d[, 3:6] / rowSums(d[, 3:6])
    d$ppv <- d$tp / (d$tp + d$fp)
    d$fdr <- 1 - d$ppv
    print(range(d$fdr))
    print(d[which(d$fdr == min(d$fdr)), ])
    print("")
    # print(d[which(d$instability == "TP53"), ])
    print(d[which(d$instability == "MYC"), ])
    # print(range(d$fp[which(d$instability == "TP53")]))
    print("")
    # print(d[which(d$complex == "NSMCE2" & d$instability == "TP53"), ])

    x <- pivot_longer(d, c("tp", "tn", "fp", "fn"))
    p <- list()
    n <- 0
    n <- n + 1
    p[[n]] <- NULL
    for (i in instabilityGenes) {
        n <- n + 1
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
            p[[n]] <- ggplot(subx, aes(x = "", y = value, fill = name)) +
                geom_bar(stat = "identity", width = 1, color = "white", size = 0.05) +
                scale_fill_manual(
                    values = color,
                    name = " ",
                    labels = c("False negative", "True negative", "False positive", "True positive")
                ) +
                coord_polar("y", start = 0) +
                theme_void() +
                theme(
                    plot.title = element_text(size = 1),
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

cbpd <- readRDS("./data/cbioportal/cbpdDataWInst.rds")



pdf("./reviewer-addressing/fig3b.pdf", height = 9, width = 15, onefile = TRUE)
plot(generatePlotDf(cbpd, "All", color = colors, bgcolor = bgcolor))
dev.off()

print("done")
