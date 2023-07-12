# load libs
source("./R/R_rainclouds.R")
library(ggplot2)
library(ggpubr)
library(cowplot)
dir <- getwd()
source("./R/routine_tasks.R")
gene <- "NSMCE2_det"
tissue <- "Breast"
colors <- c("#DE3B1C", "#707176")
# read data
df <- readRDS("./data/cbioportal/cbpdDataWInst.rds")
df <- df[grepl("amp", df$MYC_det), ]
df$Complex <- df$isalt

pa <- list()
pp <- list()
sDf <- df[
    grep(tissue, df$tissue, ignore.case = TRUE),
    c(gene, "aneuploidyScore", "ploidy")
]
plot_df <- data.frame(
    gene = gene,
    aneuploidyScore = as.numeric(sDf$aneuploidyScore),
    ploidy = as.numeric(sDf$ploidy),
    altStat = sDf[, gene]
)
# plot_df$altStat <- as.factor(ifelse(plot_df$altStat == "1", "Altered", "Wild"))
plot_df$altStat <- ifelse(grepl("amp", plot_df$altStat, ignore.case = TRUE), "Amplified", "Not Amplified")
a_plot_df <- na.omit(data.frame(
    aneuploidyScore = plot_df$aneuploidyScore,
    altStat = plot_df$altStat, gene = plot_df$gene
))
p_plot_df <- na.omit(data.frame(
    ploidy = plot_df$ploidy,
    altStat = plot_df$altStat, gene = plot_df$gene
))

# custom p value text. uncomment to add custom p value
pa <- ggplot(a_plot_df, aes(
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
pp <- ggplot(p_plot_df, aes(
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
# custom p value test. uncomment to add custom p value
# pa[[i]] <- pa[[i]] + annotate(geom="text", x=1.5, y=38, color="black", size = 4,label=p2)
# pp[[i]] <- pp[[i]] + annotate(geom="text", x=1.5, y=6.7, color="black", size = 4,label=p4)
pa <- pa + guides(color = "none", fill = "none")
pp <- pp + guides(color = "none", fill = guide_legend(title = "Status"))



pdf("./reviewer-addressing/ploidy.pdf", width = 6.5, height = 5)
plot(
    ggarrange(pa, pp, nrow = 1, ncol = 2, align = "hv", common.legend = TRUE) +
        theme(
            legend.position = "left",
            legend.direction = "vertical"
        )
)
dev.off()
print("done")
