# load libs
source("./R/R_rainclouds.R")
library(ggplot2)
library(ggpubr)
library(cowplot)
dir <- getwd()
source("./R/routine_tasks.R")
genes <- c("NSMCE2", "Complex")
tissue <- "Breast"
colors <- c("#DE3B1C", "#707176")
# read data
complexGenes <- c("SMC5", "SMC6", "NSMCE1", "NSMCE2", "NSMCE3", "NSMCE4A", "EID3", "SLF1", "SLF2")
# table(unlist(df[, paste0(complexGenes, "_det")]))
df <- readRDS("./data/cbioportal/format_exOther.rds")
df$NSMCE2 <- df$NSMCE2_det
# df$Complex <- df$isalt
df$Complex <- apply(df, 1, function(x) {
    return(paste(x[paste0(complexGenes, "_det")], collapse = ";"))
})
for (i in 1:10) {
    df$Complex <- gsub(";;", ";", df$Complex)
    df$Complex <- gsub("^;", "", df$Complex)
    df$Complex <- gsub("^;", "", df$Complex)
    df$Complex <- gsub("NA", "", df$Complex)
    df$Complex[df$Complex %in% c(";;", ";")] <- ""
}


pa <- list()
pp <- list()
plot_df <- data.frame(
    gene = NA, aneuploidyScore = NA,
    ploidy = NA, altStat = NA
)
for (gene in genes) {
    sDf <- df[
        grep(tissue, df$tissue, ignore.case = TRUE),
        c(gene, "aneuploidyScore", "ploidy")
    ]
    # NOTE: have to change this
    # altStat <- sDf[, gene]
    altStat <- rep(1, nrow(sDf))
    altStat[sDf[, gene] == ""] <- 0
    altStat[grep("Amplification", sDf[, gene])] <- NA
    print(gene)

    plot_df <- rbind(
        plot_df,
        data.frame(
            gene = gene,
            aneuploidyScore = as.numeric(sDf$aneuploidyScore),
            ploidy = as.numeric(sDf$ploidy),
            altStat = altStat
        )
    )
}
plot_df$altStat <- as.factor(ifelse(plot_df$altStat == "1", "Altered", "Wild"))
a_plot_df <- na.omit(data.frame(
    aneuploidyScore = plot_df$aneuploidyScore,
    altStat = plot_df$altStat, gene = plot_df$gene
))
print(table(a_plot_df$altStat, a_plot_df$gene))
p_plot_df <- na.omit(data.frame(
    ploidy = plot_df$ploidy,
    altStat = plot_df$altStat, gene = plot_df$gene
))
print(table(p_plot_df$altStat, p_plot_df$gene))

for (i in c("Complex", "NSMCE2")) {
    pa[[i]] <- ggplot(a_plot_df[a_plot_df$gene == i, ], aes(
        x = altStat, y = aneuploidyScore,
        fill = altStat, colour = altStat
    )) +
        geom_flat_violin(
            position = position_nudge(x = .25, y = 0),
            adjust = 2, trim = TRUE
        ) +
        geom_point(position = position_jitter(width = .15), size = .25) +
        geom_boxplot(aes(x = as.numeric(altStat) + 0.25, y = aneuploidyScore),
            outlier.shape = NA, alpha = 0.3, width = .1, colour = "BLACK"
        ) +
        ylab("Aneuploidy score") +
        xlab("") +
        ylim(0, 40) +
        ggtitle(i) +
        scale_colour_manual(values = colors) +
        scale_fill_manual(values = colors) +
        theme_cowplot() +
        # wilcox t test raw output. comment to add custom p value
        stat_compare_means(method = "wilcox.test") +
        theme(plot.title = element_text(hjust = 0.5, size = 12, face = "plain"))
    # ploidy
    pp[[i]] <- ggplot(p_plot_df[p_plot_df$gene == i, ], aes(
        x = altStat, y = ploidy,
        fill = altStat, colour = altStat
    )) +
        geom_flat_violin(
            position = position_nudge(x = .25, y = 0),
            adjust = 2, trim = TRUE
        ) +
        geom_point(position = position_jitter(width = .15), size = .25) +
        geom_boxplot(aes(x = as.numeric(altStat) + 0.25, y = ploidy),
            outlier.shape = NA, alpha = 0.3, width = .1, colour = "BLACK"
        ) +
        ylab("Ploidy score") +
        xlab("") +
        ylim(1, 7) +
        ggtitle(i) +
        scale_colour_manual(values = colors) +
        scale_fill_manual(values = colors) +
        theme_cowplot() +
        # wilcox t test raw output. comment to add custom p value
        stat_compare_means(method = "wilcox.test") +
        theme(plot.title = element_text(hjust = 0.5, size = 12, face = "plain"))

    if (i == "Complex") {
        pa[[i]] <- pa[[i]] + guides(fill = "none", color = "none")
        pp[[i]] <- pp[[i]] + guides(fill = "none", color = "none")
    } else {
        pa[[i]] <- pa[[i]] + guides(color = "none", fill = guide_legend(title = "Status"))
        pp[[i]] <- pp[[i]] + guides(color = "none", fill = guide_legend(title = "Status"))
    }
}

pa <- ggarrange(pa[[1]], pa[[2]], ncol = 2, nrow = 1, widths = c(0.4, 0.6))
pp <- ggarrange(pp[[1]], pp[[2]], ncol = 2, nrow = 1, widths = c(0.4, 0.6))

# pdf("./reviewer-addressing/plot/newFig3ab.pdf", width = 13, height = 5)
pdf("./reviewer-addressing/plot/noAmp-Fig3ab.pdf", width = 14, height = 5)
plot(ggarrange(pa, pp, nrow = 1, ncol = 2, align = "hv"))
dev.off()
