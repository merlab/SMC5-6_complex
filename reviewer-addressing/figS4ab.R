# code taken from:
# https://stackoverflow.com/questions/3483203/create-a-boxplot-in-r-that-labels-a-box-with-the-sample-size-n
# rain cloud plots from
# https://github.com/RainCloudPlots/RainCloudPlots
library(ComplexHeatmap)
library(ggplot2)
library(ggpubr)
library(survminer)
library(survival)
library(dplyr)
library(cowplot)
source("./R/R_rainclouds.R")
dir <- getwd()
source("./R/routine_tasks.R")
genes <- c("NSMCE2", "Complex")
tissue <- "Prostate"
colors <- c("#DE3B1C", "#707176")


df <- readRDS("./reviewer-addressing/cbioportal/formatted.rds")
df$Complex <- df$isalt

pa <- list()
pp <- list()
plot_df <- data.frame(gene = NA, aneuploidyScore = NA
                , ploidy = NA, altStat = NA)
for (gene in genes) {
    sDf <- df[
        grep(tissue, df$tissue, ignore.case = TRUE)
        , c(gene, "aneuploidyScore", "ploidy")]
    plot_df <- rbind(plot_df,
                data.frame(gene = gene,
                    aneuploidyScore = as.numeric(sDf$aneuploidyScore),
                    ploidy = as.numeric(sDf$ploidy),
                    altStat = (sDf[, gene]))
    )
}
plot_df$altStat <- as.factor(ifelse(plot_df$altStat == "1", "Altered", "Wild"))
a_plot_df <- na.omit(data.frame(aneuploidyScore = plot_df$aneuploidyScore
                , altStat = plot_df$altStat, gene = plot_df$gene))
p_plot_df <- na.omit(data.frame(ploidy = plot_df$ploidy
                , altStat = plot_df$altStat, gene = plot_df$gene))
print(max(plot_df$ploidy, na.rm = TRUE))
print(max(plot_df$aneuploidyScore, na.rm = TRUE))

for (i in c("Complex", "NSMCE2")) {
  pa[[i]] <- ggplot(a_plot_df[a_plot_df$gene == i, ], aes(x=altStat, y= aneuploidyScore
              , fill = altStat, colour = altStat))+
    geom_flat_violin(position = position_nudge(x = .25, y = 0)
          , adjust = 2, trim = TRUE)+
    geom_point(position = position_jitter(width = .15), size = .25)+
    geom_boxplot(aes(x = as.numeric(altStat)+0.25, y = aneuploidyScore)
          ,outlier.shape = NA, alpha = 0.3, width = .1, colour = "BLACK") +
    ylab("Aneuploidy score") + xlab("") +
    ylim(0, 30) +
    ggtitle(i) +
    scale_colour_manual(values = colors) + scale_fill_manual(values = colors) +
        stat_compare_means(method = "wilcox.test") +
    theme_cowplot() +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "plain"))
  # ploidy
  pp[[i]] <- ggplot(p_plot_df[p_plot_df$gene == i, ], aes(x=altStat, y= ploidy
              , fill = altStat, colour = altStat))+
    geom_flat_violin(position = position_nudge(x = .25, y = 0)
          , adjust = 2, trim = TRUE) +
    geom_point(position = position_jitter(width = .15), size = .25)+
    geom_boxplot(aes(x = as.numeric(altStat)+0.25, y = ploidy)
          , outlier.shape = NA, alpha = 0.3, width = .1, colour = "BLACK") +
    ylab("Ploidy score") + xlab("") +
    ylim(1, 6) +
    ggtitle(i) +
    scale_colour_manual(values = colors) + scale_fill_manual(values = colors) +
    theme_cowplot() +
        stat_compare_means(method = "wilcox.test") +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "plain"))

  if (i == "Complex") {
    pa[[i]] <- pa[[i]] + guides(fill = "none", color = "none")
    pp[[i]] <- pp[[i]] + guides(fill = "none", color = "none")
  } else {
    pa[[i]] <- pa[[i]] + guides(color = "none", fill=guide_legend(title="Status"))
    pp[[i]] <- pp[[i]] + guides(color = "none", fill=guide_legend(title="Status"))
  }
}

pa <- ggarrange(pa[[1]], pa[[2]], ncol = 2, nrow = 1, widths = c(0.4, 0.6))
pp <- ggarrange(pp[[1]], pp[[2]], ncol = 2, nrow = 1, widths = c(0.4, 0.6))

###################
###### PART 2 #####
###################
studies <- c("Prostate Adenocarcinoma (TCGA, PanCancer Atlas)")
# columns used for analysis
cols <- c(isalt = "Complex", NSMCE2 = "NSMCE2")
# color of the plot
ref_df <- as.data.frame(readRDS(sprintf("./data/cbioportal/formatted.rds")))

sur_dfs <- list()
titles <- list()
for (current_study in studies) {
    print(current_study)

    df <- ref_df[ref_df$study == current_study, ]
    for (l in 1:length(cols)) {
        col <- names(cols[l])
        namecol <- as.character(cols[l])
        print(namecol)
        sur_df <- data.frame(study = df$study, tissue = df$tissue, OVT= df$OVT,
                        OVS = df$OVS, group = df[, col])
        sur_df$OVT <- as.numeric(sur_df$OVT)

        sur_df$OVS <- as.numeric(sur_df$OVS)
        sur_df$group <- as.numeric(sur_df$group)
        sur_df$group <- ifelse(sur_df$group == 1, "Altered", "Wild") #namecol
        # The top row are for breast cancer, the bottom row is for prostate cancer. #
        sur_dfs[[paste(current_study, namecol)]] <- na.omit(sur_df)
        titles[[paste(current_study, namecol)]] <- namecol

    }
}

pdf("./reviewer-addressing/figS4ab.pdf", width = 13, height = 5)
plot(ggarrange(pa, pp, nrow = 1, ncol = 2, align = "hv"))
dev.off()

print("done")
