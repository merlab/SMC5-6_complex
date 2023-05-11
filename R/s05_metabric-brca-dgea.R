row <- "isalt"
library(limma)
library(ggplot2)
library(ggpubr)
source("./R/routine_tasks.R")
expmat <- readRDS("./data/metabric-brca/microarray-metagx.rds")
samples <- colnames(expmat)
cbpd <- readRDS("./data/cbioportal/format_exOther.rds")
rownames(cbpd) <- cbpd$name
cbpd <- cbpd[cbpd$name %in% samples, ]
expmat <- expmat[, colnames(expmat) %in% cbpd$name]
cbpd <- cbpd[colnames(expmat), ]
expmat <- expmat[, cbpd$name]

cbpd[, row] <- (ifelse(cbpd[, row] == "1", "Mutated", "Wild"))
cbpd[, row] <- factor(cbpd[, row], levels = c("Wild", "Mutated"))
er <- as.numeric(cbpd[, row])
design <- model.matrix(~er)
design <- as.matrix(design)
row.names(design) <- cbpd$name
fit <- lmFit(expmat, design)
fitted.ebayes <- eBayes(fit)

out <- topTable(fitted.ebayes, number = Inf)

out$logFDR <- -log10(out$adj.P.Val)
out$Significance <- "Not Significant"
out$Significance[out$adj.P.Val < 0.01] <- "Significant"
out$gene <- rownames(out)
out <- out[order(abs(out$logFC), decreasing = TRUE), ]

folder_check("./data/metabric-brca")
saveRDS(out, "./data/metabric-brca/dgea-limma.rds")

out <- out[order(out$adj.P.Val, decreasing = FALSE), ]
signatureGenes <- out$gene[seq_len(10)]
write.table(signatureGenes, file = "./results/metabricSignature.txt", quote = FALSE, sep = "\n")
