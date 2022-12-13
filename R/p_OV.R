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
source('./R/R_rainclouds.R')
dir <- getwd()
source("./R/routine_tasks.R")
genes <- c("NSMCE2", "Complex")
tissue <- "Ovarian"
colors <- c("#CB3814", "#3361BD")

df <- readRDS("./results/cbioportal_alt_all.rds")
df$Complex <- df$isalt

pa <- list()
pp <- list()
plot_df <- data.frame(gene = NA, aneuploidyScore = NA
                , ploidy = NA, altStat = NA)
for (gene in genes) {
    sDf <- df[
        grep(tissue, df$tissue, ignore.case = TRUE)
        ,c(gene, "aneuploidyScore", "ploidy")]
    plot_df <- rbind(plot_df, 
                data.frame(gene = gene, 
                    aneuploidyScore = as.numeric(sDf$aneuploidyScore), 
                    ploidy = as.numeric(sDf$ploidy),
                    altStat = (sDf[,gene]))
    )
}
plot_df$altStat <- as.factor(ifelse(plot_df$altStat == "1", "Altered", "Wild"))
a_plot_df <- na.omit(data.frame(aneuploidyScore = plot_df$aneuploidyScore
                , altStat = plot_df$altStat, gene = plot_df$gene))
p_plot_df <- na.omit(data.frame(ploidy = plot_df$ploidy
                , altStat = plot_df$altStat, gene = plot_df$gene))
print(max(plot_df$ploidy, na.rm = TRUE))
print(max(plot_df$aneuploidyScore, na.rm = TRUE))

# p1 <- bquote("P = " ~ '1' ~ 'x' ~ 10^-7)
# p2 <- bquote("P = " ~ '3' ~ 'x' ~ 10^-7)
# p3 <- bquote("P = " ~ '2' ~ 'x' ~ 10^-5)
# p4 <- bquote("P = " ~ '8' ~ 'x' ~ 10^-3)
for(i in c('Complex', 'NSMCE2')) {
  pa[[i]] <- ggplot(a_plot_df[a_plot_df$gene == i,], aes(x=altStat, y= aneuploidyScore
              , fill = altStat, colour = altStat))+
    geom_flat_violin(position = position_nudge(x = .25, y = 0)
          ,adjust = 2, trim = TRUE)+
    geom_point(position = position_jitter(width = .15), size = .25)+
    geom_boxplot(aes(x = as.numeric(altStat)+0.25, y = aneuploidyScore)
          ,outlier.shape = NA, alpha = 0.3, width = .1, colour = "BLACK") +
    ylab('Aneuploidy score') + xlab('') +
    ylim(0,30) +
    ggtitle(i) +
    scale_colour_manual(values = colors) + scale_fill_manual(values = colors) +
    stat_compare_means(method = "wilcox.test") + 
    theme_cowplot() +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'plain'))
  # ploidy
  pp[[i]] <- ggplot(p_plot_df[p_plot_df$gene == i,], aes(x=altStat, y= ploidy
              , fill = altStat, colour = altStat))+
    geom_flat_violin(position = position_nudge(x = .25, y = 0)
          , adjust = 2, trim = TRUE) +
    geom_point(position = position_jitter(width = .15), size = .25)+
    geom_boxplot(aes(x = as.numeric(altStat)+0.25, y = ploidy)
          , outlier.shape = NA, alpha = 0.3, width = .1, colour = "BLACK") +
    ylab('Ploidy score') + xlab('') +
    ylim(1,6) +
    ggtitle(i) +
    scale_colour_manual(values = colors) + scale_fill_manual(values = colors) + 
    theme_cowplot() +
    stat_compare_means(method = "wilcox.test") + 
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'plain'))

  if(i == 'Complex') {
    #pa[[i]] <- pa[[i]] + annotate(geom="text", x=1.5, y=28, color="black", size = 4,label=p1)
    #pp[[i]] <- pp[[i]] + annotate(geom="text", x=1.5, y=5.7, color="black", size = 4,label=p3)
    pa[[i]] <- pa[[i]] + guides(fill = 'none', color = 'none')
    pp[[i]] <- pp[[i]] + guides(fill = 'none', color = 'none')
  } else {
    #pa[[i]] <- pa[[i]] + annotate(geom="text", x=1.5, y=28, color="black", size = 4,label=p2)
    #pp[[i]] <- pp[[i]] + annotate(geom="text", x=1.5, y=5.7, color="black", size = 4,label=p4)
    pa[[i]] <- pa[[i]] + guides(color = 'none', fill=guide_legend(title="Status")) 
    pp[[i]] <- pp[[i]] + guides(color = 'none', fill=guide_legend(title="Status")) 
  }
}

pa <- ggarrange(pa[[1]], pa[[2]], ncol = 2, nrow = 1, widths = c(0.4, 0.6))
pp <- ggarrange(pp[[1]], pp[[2]], ncol = 2, nrow = 1, widths = c(0.4, 0.6))

pdf("./figures/p_OV.pdf", width = 13, height = 5)
plot(ggarrange(pa, pp, nrow = 1, ncol = 2, align = 'hv'))
dev.off()
print('done')