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
tissue <- c("Breast")
colors <- c("#CB3814", "#3361BD")

KM_survival_plot <- function(sur_df, colors, title, xlab = TRUE, ylab = TRUE, strata = FALSE) {
  gc()
  #censoring function
  censor <- 60
  sur_df$OVS[sur_df$OVT > censor] <- 0
  sur_df$OVT[sur_df$OVT > censor] <- censor

  diff <- survdiff(Surv(OVT, OVS)~ group, data = sur_df)
  p.val <- 1 - pchisq(diff$chisq, length(diff$n) - 1)

  title <- title
  fit <- NA
  fit <- survfit(Surv(OVT, OVS)~ group, data = sur_df)
  # https://stat.ethz.ch/pipermail/r-help/2007-April/130676.html
  # survival plot
  ggsurv <- 
    ggsurvplot(fit
        , conf.int = TRUE
        , pval = FALSE
        , palette = colors
        , xlab = '' 
        , ylab = ifelse(ylab, "Probability of overall survival",'')
        , title = title
        , legend.title='Status'
        , legend = c(.85,.4)
        , legend.labs = c('Altered', 'Wild')
        , risk.table = TRUE
        , axes.offset = FALSE
        , risk.table.height = 0.22
  ) 
  ggsurv$table <- ggrisktable(fit
        , data = sur_df
        , ylab = ''
        , xlab = ifelse(xlab, "Time (Months)", '')
        , risk.table.title = ''
        , palette = colors
        , color = 'strata'
        , legend = 'none'
        , axes.offset = TRUE
        , tables.theme = theme_classic()
        , fontsize = 3.25
        , risk.table.col = "strata"
        # NEEDED TO GET THE RISK TABLE CORRECTLY
        , break.time.by = 10
  ) + scale_y_discrete(labels = c('Wild', 'Altered')) +
      theme(axis.text.x = element_text(size = 12),
            axis.title.x = element_text(size = 12),
            legend.text = element_blank(),
            legend.title = element_blank())
  ggsurv$plot <- ggsurv$plot + theme(plot.title = element_text(hjust = 0.5, size = 12)) 
  ggsurv$plot <- ggsurv$plot + theme(axis.text.x = element_blank())
  ggsurv$plot <- ggsurv$plot + theme(plot.title = element_text(hjust = 0.5, size = 12)) 
  ggsurv$plot <- ggsurv$plot + theme(plot.margin = unit(c(5, 5, 0, 5), "points"))
  ggsurv$table <- ggsurv$table + theme(plot.margin = unit(c(5, 5, 0, 5), "points"))
  ggsurv$plot <- ggsurv$plot + scale_x_continuous(limits = c(0,65), breaks = seq(0, 60, 10))
  ggsurv$table <- ggsurv$table + scale_x_continuous(limits = c(0,65), breaks = seq(0, 60, 10))
   if(signif(p.val,1) == 1e-6) {
     print(p.val)
     p.val_text <- bquote("P = " ~ '1' ~ 'x' ~ 10^-6)
   }
   if(signif(p.val,1) == 5e-4) {
     print(p.val)
     p.val_text <- bquote("P = " ~ '5' ~ 'x' ~ 10^-4)
   }
  ggsurv$plot <- ggsurv$plot + annotate(
      geom="text", x=10, y=0.1,
              color="black", size = 5,
              label=p.val_text
              )
  # this is to bring risk table up
  rm(fit, sur_df)
  fit <- NA
  return(ggsurv)
}

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

p1 <- bquote("P = " ~ '3' ~ 'x' ~ 10^-4)
p2 <- bquote("P = " ~ '4' ~ 'x' ~ 10^-4)
p3 <- bquote("P = " ~ '1' ~ 'x' ~ 10^-6)
p4 <- bquote("P = " ~ '4' ~ 'x' ~ 10^-5)
for(i in c('Complex', 'NSMCE2')) {
  pa[[i]] <- ggplot(a_plot_df[a_plot_df$gene == i,], aes(x=altStat, y= aneuploidyScore
              , fill = altStat, colour = altStat))+
    geom_flat_violin(position = position_nudge(x = .25, y = 0)
          ,adjust = 2, trim = TRUE)+
    geom_point(position = position_jitter(width = .15), size = .25)+
    geom_boxplot(aes(x = as.numeric(altStat)+0.25, y = aneuploidyScore)
          ,outlier.shape = NA, alpha = 0.3, width = .1, colour = "BLACK") +
    ylab('Aneuploidy score') + xlab('') +
    ylim(0,40) +
    ggtitle(i) +
    scale_colour_manual(values = colors) + scale_fill_manual(values = colors) +
    theme_cowplot() +
    # stat_compare_means(method = "wilcox.test") + 
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
    ylim(1,7) +
    ggtitle(i) +
    scale_colour_manual(values = colors) + scale_fill_manual(values = colors) + 
    theme_cowplot() +
    # stat_compare_means(method = "wilcox.test") + 
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'plain'))

  if(i == 'Complex') {
    pa[[i]] <- pa[[i]] + annotate(geom="text", x=1.5, y=38, color="black", size = 4,label=p1)
    pp[[i]] <- pp[[i]] + annotate(geom="text", x=1.5, y=6.7, color="black", size = 4,label=p3)
    pa[[i]] <- pa[[i]] + guides(fill = 'none', color = 'none')
    pp[[i]] <- pp[[i]] + guides(fill = 'none', color = 'none')
  } else {
    pa[[i]] <- pa[[i]] + annotate(geom="text", x=1.5, y=38, color="black", size = 4,label=p2)
    pp[[i]] <- pp[[i]] + annotate(geom="text", x=1.5, y=6.7, color="black", size = 4,label=p4)
    pa[[i]] <- pa[[i]] + guides(color = 'none', fill=guide_legend(title="Status")) 
    pp[[i]] <- pp[[i]] + guides(color = 'none', fill=guide_legend(title="Status")) 
  }
}

pa <- ggarrange(pa[[1]], pa[[2]], ncol = 2, nrow = 1, widths = c(0.4, 0.6))
pp <- ggarrange(pp[[1]], pp[[2]], ncol = 2, nrow = 1, widths = c(0.4, 0.6))

###################
## Configuration ##
###################
studies <- c('Breast Cancer (METABRIC, Nature 2012 & Nat Commun 2016)'
    )
# columns used for analysis
cols <- c(isalt = "Complex", NSMCE2 = "NSMCE2")
# color of the plot
ref_df <- as.data.frame(readRDS(sprintf("./results/cbioportal_alt_all.rds")))

sur_dfs <- list()
titles <- list()
for (current_study in studies) {
    print(current_study)

    df <- ref_df[ref_df$study == current_study,]
    for (l in 1:length(cols)) {
      col <- names(cols[l])
      namecol <- as.character(cols[l])
      print(namecol)
      sur_df <- data.frame(study = df$study, tissue = df$tissue,OVT= df$OVT, 
                      OVS = df$OVS, group = df[, col])
      sur_df$OVT <- as.numeric(sur_df$OVT)

      sur_df$OVS <- as.numeric(sur_df$OVS)
      sur_df$group <- as.numeric(sur_df$group)
      sur_df$group <- ifelse(sur_df$group == 1, 'Altered', "Wild") #namecol
      sur_dfs[[paste(current_study,namecol)]] <- na.omit(sur_df)
      titles[[paste(current_study, namecol)]] <- namecol

      }
}

pdf("./figures/fig3row1.pdf", width = 13, height = 5)
plot(ggarrange(pa, pp, nrow = 1, ncol = 2, align = 'hv'))
dev.off()

pdf("./figures/fig3row2.pdf", width = 10, height = 5)
p <- list((KM_survival_plot(sur_dfs[[1]], color = colors, title = titles[[1]]))
            , (KM_survival_plot(sur_dfs[[2]], color = colors, title = titles[[2]])))
arrange_ggsurvplots(p,
  print = TRUE,
  ncol = 2,
  nrow = 1
)
dev.off()
print('done')