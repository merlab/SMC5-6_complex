row <- 'isalt'
library(limma)
library(ggplot2)
library(ggpubr)
library(writexl)
source('./R/routine_tasks.R')
expmat <- readRDS('./results/metabric-brca/metagx_microarray.rds')
samples <- colnames(expmat)
cbpd <- readRDS("./results/cbioportal_alt_all.rds")
rownames(cbpd) <- cbpd$name
cbpd <- cbpd[cbpd$name %in% samples, ]
expmat <- expmat[, colnames(expmat) %in% cbpd$name]
cbpd <- cbpd[colnames(expmat),]
expmat <- expmat[, cbpd$name]

cbpd[, row] <- (ifelse(cbpd[,row]== "1", "Mutated", "Wild"))
cbpd[, row] <- factor(cbpd[,row], levels = c("Wild", 'Mutated'))
er <- as.numeric(cbpd[,row])
design <- model.matrix(~ er)
design <- as.matrix(design)
row.names(design) <- cbpd$name
fit <- lmFit(expmat, design)
fitted.ebayes <- eBayes(fit)
  
out <- topTable(fitted.ebayes, number = Inf)

out$logFDR <- -log10(out$adj.P.Val)
out$Significance <- "Not Significant"
out$Significance[out$adj.P.Val < 0.01] <- "Significant"
out$gene <- rownames(out)
out <- out[order(abs(out$logFC), decreasing = TRUE),]

saveRDS(out, './results/metabric-brca/limma.rds')
df <- out
df <- df[order(abs(df$logFC), decreasing = TRUE), ]
df$FDR <- df$adj.P.Val
df$adj.P.Val <- NULL
df$Significance <- NULL
df$t <- NULL
df$B <- NULL
df$name <- NULL
df <- df[, c('gene', 'logFC', 'AveExpr', 'P.Value', 'FDR')]
write_xlsx(df, './results/metabric-brca/differential_gene_expr_brca_metabric.xlsx')
print("done")
