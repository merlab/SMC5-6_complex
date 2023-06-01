library(survminer)
library(survival)
library(ggpubr)
library(data.table)
source("./R/routine_tasks.R")
KM_plot_PFS <- function(df, title) {
    df <- df[!is.na(df$PFT), ]
    df <- df[!is.na(df$PFS), ]
    fitisalt <- survfit(Surv(PFT, PFS) ~ isalt, data = df)
    fitnsmce2 <- survfit(Surv(PFT, PFS) ~ NSMCE2, data = df)
    # palette <- rev(c("#CB3814", "#3361BD"))
    palette <- c("#707176", "#DE3B1C")
    pisalt <- surv_pvalue(fitisalt)$pval
    pnsmce2 <- surv_pvalue(fitnsmce2)$pval
    print(pisalt)
    print(pnsmce2)
    pl <- list()
    pl[[1]] <- ggsurvplot(fitisalt,
        conf.int = TRUE,
        risk.table = FALSE, palette = palette, legend = c(0.8, 0.95),
        pval = FALSE # TRUE
        , xlab = "Progression free survival (Months)",
        ylab = "Probability of Progression free survival",
        title = "Complex",
        legend.title = "Complex",
        legend.labs = c("Wild", "Altered")
    )
    pl[[2]] <- ggsurvplot(fitnsmce2,
        conf.int = TRUE,
        risk.table = FALSE, palette = palette, legend = c(0.8, 0.95),
        pval = FALSE # TRUE
        , xlab = "Progression free survival (Months)",
        ylab = "Probability of Progression free survival",
        title = "NSMCE2",
        legend.title = "NSMCE2",
        legend.labs = c("Wild", "Altered")
    )
    for (i in 1:length(pl)) {
        if (i == 1) p.val <- pisalt
        if (i == 2) p.val <- pnsmce2
        print(p.val)
        if (signif(p.val, 2) == 0.57) p.val_text <- bquote("P = " ~ "0.57")
        if (signif(p.val, 2) == 0.38) p.val_text <- bquote("P = " ~ "0.38")
        if (signif(p.val, 2) == 0.0098) p.val_text <- bquote("P = " ~ "0.001")
        if (signif(p.val, 2) == 0.047) p.val_text <- bquote("P = " ~ "0.05")
        pl[[i]]$plot <- pl[[i]]$plot + annotate(
            geom = "text", x = max(df$PFT) / 10, y = 0.05,
            color = "black", size = 5,
            label = p.val_text
        )

        pl[[i]]$plot <- pl[[i]]$plot + theme(plot.title = element_text(hjust = 0.5, size = 12))
    }
    arrange_ggsurvplots(pl, print = TRUE) # , title = paste(title, 'PFS'))
}

### formatting
cl <- readRDS("./data/tcga-ov/clinical.rds")
df <- as.data.frame(data.table::fread("./data/tcga-ov-2011.tsv", header = TRUE))
cohort2011 <- colnames(df)[-c(1, 2)]
cohort2011 <- cohort2011[cohort2011 %in% rownames(cl)]
cl$PFS <- as.numeric(cl$PFS)
cl$PFT <- as.numeric(cl$PFT)
cl$DFS <- as.numeric(cl$DFS)
cl$DFT <- as.numeric(cl$DFT)
cl$OVS <- as.numeric(cl$OVS)
cl$OVT <- as.numeric(cl$OVT)

# NOTE: all data
pdf("./figures/figS3.pdf", width = 12, height = 6, onefile = TRUE)
# naming
m <- c("")
n <- c("")
q <- c("")
v <- c(m, n, q)
v <- c(v, "all ages")
name <- paste(c(v), collapse = "-")

# data formatting
db <- cl
df <- db
KM_plot_PFS(df, name)

# NOTE: SUBSET
# naming
m <- "TP53 mutation"
n <- "No radiotherapy"
q <- "2011 cohort"
age1 <- 55
age2 <- "<"
v <- c(m, n, q)
v <- c(v, paste0("age ", age2, age1))
v <- v[v != ""]
name <- paste(c(v), collapse = "-")

# data formatting
db <- cl
db <- db[db$rtherapy == "No", ]
db <- db[cohort2011, ]
db <- db[db$TP53 == 1, ]
db$age <- as.factor(ifelse(db$age < age1, "<", ">="))
db <- db[db$age == age2, ]
df <- db
KM_plot_PFS(df, name)
dev.off()

