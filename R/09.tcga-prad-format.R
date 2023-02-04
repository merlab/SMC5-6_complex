# purpose: prep the STAR-count of TCGA biolinker for edgeR analysis
# max percentage of patients with 0 count allowed for each genes:
source('./R/routine_tasks.R')
library(Biobase)
a <- readRDS('./data/tcga-prad-star-counts.rds')
print(dim(a))
a <- a[, colData(a)$sample_type == 'Primary Tumor']
a <- a[rowData(a)$gene_type == 'protein_coding',]
print(dim(a))
colnames(a) <- colData(a)$patient
colnames(a) <- gsub('\\.', '-', colnames(a))
print(dim(a))
rvar <- (apply(assays(a)[['unstranded']], 1, var))
rvar <- rvar[!is.na(rvar) & rvar != 0]
a <- a[names(rvar),]
print(dim(a))

# remove patients with 0 variance or NA count
cvar <- (apply(assays(a)[['unstranded']], 2, var))
cvar <- cvar[!is.na(cvar) & cvar != 0]
a <- a[ ,names(cvar)]
print(dim(a))

duplicatedGenes <- unique(rowData(a)$gene_name[duplicated(rowData(a)$gene_name)])
a <- a[!rowData(a)$gene_name %in% duplicatedGenes, ]
rownames(a) <- make.unique(rowData(a)$gene_name)
keep <- rownames(a) %in% rowData(a)$gene_name
a <- a[keep, ]
print(dim(a))

folder_check('./data/tcga-prad')
saveRDS(assays(a)[['unstranded']], './data/tcga-prad/rnaseq.rds')
print('done')