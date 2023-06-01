###################
## Configuration ##
###################
studies <- c("Breast Cancer (METABRIC, Nature 2012 & Nat Commun 2016)")
# columns used for analysis
cols <- c(isalt = "Complex", NSMCE2 = "NSMCE2")
# color of the plot
ref_df <- readRDS("./data/cbioportal/format_exOther.rds")
# code taken from:
# https://stackoverflow.com/questions/3483203/create-a-boxplot-in-r-that-labels-a-box-with-the-sample-size-n
# rain cloud plots from
# https://github.com/RainCloudPlots/RainCloudPlots
# libs
library(ComplexHeatmap)
library(ggplot2)
library(survminer)
library(survival)

KM_survival_plot <- function(sur_df, title, xlab = TRUE, ylab = TRUE, strata = FALSE, censor = TRUE) {
  colors <- c("#DE3B1C","#707176")
  gc()
  #censoring function
  if(censor == TRUE) {
    sur_df$OVS[sur_df$OVT > 60] <- 0
    sur_df$OVT[sur_df$OVT > 60] <- 60 
  }

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
        , xlab = "" 
        , ylab = ifelse(ylab, "Probability of overall survival","")
        , title = title
        , legend.title="Status"
        , legend = c(.85,.4)
        , legend.labs = c("Altered", "Wild")
        , risk.table = TRUE
        , axes.offset = FALSE
        , risk.table.height = 0.22
  ) 
  ggsurv$table <- ggrisktable(fit
        , data = sur_df
        , ylab = ""
        , xlab = ifelse(xlab, "Time (Months)", "")
        , risk.table.title = ""
        , palette = colors
        , color = "strata"
        , legend = "none"
        , axes.offset = TRUE
        , tables.theme = theme_classic()
        , fontsize = 3.25
        , risk.table.col = "strata"
        # NEEDED TO GET THE RISK TABLE CORRECTLY
        , break.time.by = ifelse(censor,10,50)
  ) + scale_y_discrete(labels = c("Wild", "Altered")) +
      theme(axis.text.x = element_text(size = 12),
            axis.title.x = element_text(size = 12, face = "bold"),
            axis.title.y = element_text(size = 12, face = "bold"),
            legend.text = element_blank(),
            legend.title = element_blank())
  ggsurv$plot <- ggsurv$plot + theme(plot.title = element_text(hjust = 0.5, size = 12)) 
  ggsurv$plot <- ggsurv$plot + theme(axis.text.x = element_blank())
  ggsurv$plot <- ggsurv$plot + theme(plot.title = element_text(hjust = 0.5, size = 12)) 
  ggsurv$plot <- ggsurv$plot + theme(plot.margin = unit(c(5, 5, 0, 5), "points"))
  ggsurv$table <- ggsurv$table + theme(plot.margin = unit(c(5, 5, 0, 5), "points"))
  if(censor == TRUE) {
    ggsurv$plot <- ggsurv$plot + scale_x_continuous(limits = c(0,65), breaks = seq(0, 60, 10))
    ggsurv$table <- ggsurv$table + scale_x_continuous(limits = c(0,65), breaks = seq(0, 60, 10))
  } else {
    ggsurv$table <- ggsurv$table + scale_x_continuous(limits = c(0,360), breaks = seq(0, 350, 50))
    ggsurv$plot <- ggsurv$plot + scale_x_continuous(limits = c(0,360), breaks = seq(0, 350, 50))
  }
  if(signif(p.val,1) == 1e-6) {
    print(p.val)
    p.val_text <- bquote("p = " ~ "1" ~ "x" ~ 10^-6)
  }
  if(signif(p.val,1) == 5e-4) {
    print(p.val)
    p.val_text <- bquote("p = " ~ "5" ~ "x" ~ 10^-4)
  }
  if(signif(p.val,1) == 1e-5) {
    print(p.val)
    p.val_text <- bquote("p = " ~ "1" ~ "x" ~ 10^-5)
  }
  if(signif(p.val,1) == 6e-4) {
    print(p.val)
    p.val_text <- bquote("p = " ~ "6" ~ "x" ~ 10^-4)
  }

  ggsurv$plot <- ggsurv$plot + annotate(
      geom="text", x=ifelse(censor,10,50), y=0.1,
              color="black", size = 3,
              label=p.val_text, face = "plain"
              )
  # this is to bring risk table up
  rm(fit, sur_df)
  fit <- NA
  return(ggsurv)
}



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
      sur_df$group <- ifelse(sur_df$group == 1, "Altered", "Wild") #namecol
      sur_dfs[[paste(current_study,namecol)]] <- na.omit(sur_df)
      titles[[paste(current_study, namecol)]] <- namecol

      }
}


pdf("./figures/fig3cd.pdf", width = 10, height = 5)
p <- list((KM_survival_plot(sur_dfs[[1]], title = "SMC5/6 complex", censor = FALSE))
            , (KM_survival_plot(sur_dfs[[2]], title = "NSMCE2", censor = FALSE)))
arrange_ggsurvplots(p,
  print = TRUE,
  ncol = 2,
  nrow = 1
)
dev.off()
print("done")
