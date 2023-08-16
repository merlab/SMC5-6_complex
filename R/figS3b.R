library(ComplexHeatmap)
library(ggplot2)
library(survminer)
library(survival)
KM_survival_plot <- function(sur_df, title, xlab = TRUE, ylab = TRUE, strata = FALSE, censor = TRUE) {
  gc()
  # censoring function
  if (censor == TRUE) {
    sur_df$OVS[sur_df$OVT > 60] <- 0
    sur_df$OVT[sur_df$OVT > 60] <- 60
  }

  diff <- survdiff(Surv(OVT, OVS) ~ group, data = sur_df)
  p.val <- 1 - pchisq(diff$chisq, length(diff$n) - 1)
  p.val_text <- paste("P = ", signif(p.val, 2))

  title <- title
  fit <- NA
  fit <- survfit(Surv(OVT, OVS) ~ group, data = sur_df)
  # https://stat.ethz.ch/pipermail/r-help/2007-April/130676.html
  # survival plot
  x <- unique(as.character(sur_df$group))
  x <- na.omit(x)
  x <- length(x)
  if (x == 2) {
    legend.labs <- c("Wild-type", "Mutation")
    colors <- c("#707176", "#377eb8")
  }
  if (x == 3) {
    legend.labs <- c("Wild-type", "Amplification", "Mutation")

    colors <- c("#707176", "#DE3B1C", "#377eb8")
  }

  ggsurv <-
    ggsurvplot(fit,
      conf.int = TRUE,
      pval = FALSE,
      censor = FALSE,
      palette = colors,
      xlab = "",
      ylab = ifelse(ylab, "Probability of overall survival", ""),
      title = title,
      legend.title = "Status",
      legend = c(.85, .2),
      legend.labs = legend.labs,
      risk.table = TRUE,
      axes.offset = FALSE,
      risk.table.height = 0.22
    )
  ggsurv$table <- ggrisktable(fit,
    data = sur_df,
    ylab = "",
    xlab = ifelse(xlab, "Time (Months)", ""),
    palette = colors,
    risk.table.title = "",
    axes.offset = TRUE,
    tables.theme = theme_classic(),
    fontsize = 3.25,
    # NEEDED TO GET THE RISK TABLE CORRECTLY
    break.time.by = ifelse(censor, 10, 50)
  ) + scale_y_discrete(labels = rev(legend.labs)) +
    theme(
      axis.text.x = element_text(size = 12),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      legend.text = element_blank(),
      legend.title = element_blank()
    )
  ggsurv$plot <- ggsurv$plot + theme(axis.text.x = element_blank())
  ggsurv$plot <- ggsurv$plot + theme(plot.title = element_text(hjust = 0.5, size = 12))
  ggsurv$plot <- ggsurv$plot + theme(plot.margin = unit(c(5, 5, 0, 5), "points"))
  ggsurv$table <- ggsurv$table + theme(plot.margin = unit(c(5, 5, 0, 5), "points"))
  if (censor == TRUE) {
    ggsurv$table <- ggsurv$table + scale_x_continuous(limits = c(0, 65), breaks = seq(0, 60, 10))
    ggsurv$plot <- ggsurv$plot + scale_x_continuous(limits = c(0, 65), breaks = seq(0, 60, 10))
  } else {
    ggsurv$table <- ggsurv$table + scale_x_continuous(limits = c(0, 360), breaks = seq(0, 350, 50))
    ggsurv$plot <- ggsurv$plot + scale_x_continuous(limits = c(0, 360), breaks = seq(0, 350, 50))
  }

  if (signif(p.val, 2) == 2.7e-7) {
    print(p.val)
    p.val_text <- bquote("P = " ~ "2.7" ~ "x" ~ 10^"-7")
  }

  if (signif(p.val, 2) == 8.7e-9) {
    print(p.val)
    p.val_text <- bquote("P = " ~ "8.7" ~ "x" ~ 10^"-9")
  }

  ggsurv$plot <- ggsurv$plot + annotate(
    geom = "text", x = ifelse(censor, 10, 50), y = 0.1,
    color = "black",
    # size = 3,
    size = 4,
    label = p.val_text
  )
  # this is to bring risk table up
  rm(fit, sur_df)
  fit <- NA
  return(ggsurv)
}
###################
## Configuration ##
###################
###
complexGenes <- c("SMC5", "SMC6", "NSMCE1", "NSMCE2", "NSMCE3", "NSMCE4A", "EID3", "SLF1", "SLF2")
complexGenesDet <- paste0(complexGenes, "_det")
studies <- c("TCGA")
# columns used for analysis
cols <- c(isalt_det = "Complex_det")
# color of the plot
ref_df <- readRDS("./data/cbioportal/cbpdDataWInst.rds")

isalt_det <- apply(ref_df, 1, function(x) {
  v <- x[complexGenesDet]
  v <- gsub(" \\(putative passenger\\)", "", v)
  v <- tolower(v)
  v <- v[!duplicated(v)]
  v <- paste(v, collapse = ";")
  v <- gsub(";;", ";", v)
  v <- gsub(";;", ";", v)
  v <- gsub(";;", ";", v)
  v <- gsub(";;", ";", v)
  v <- gsub("^;", "", v)
  v <- gsub(";$", "", v)
  v <- gsub("amplification;amplification", "amplification", v)
  v <- gsub(" \\(putative passenger\\)", "", v)
})

ref_df$isalt_det <- isalt_det


ref_df$group <- NA
ref_df$group[-grep("amplification", ref_df$isalt_det, ignore.case = TRUE)] <- "Mutation"
ref_df$group[ref_df$isalt_det == ""] <- "Wild-type"
ref_df$group[ref_df$isalt_det == "amplification"] <- "Amplification"
ref_df$group <- factor(ref_df$group, levels = c("Wild-type", "Amplification", "Mutation"))

pdf("./figures/figS3b.pdf", width = 5, height = 5)


df <- ref_df
df <- df[df$major == "Breast", ]
print(table(df$study))
title <- ""

tryCatch(expr = {
  sur_df <- df
  print(KM_survival_plot(sur_df = sur_df, title = title))
}, error = function(cond) message(cond))

dev.off()
print("done")
