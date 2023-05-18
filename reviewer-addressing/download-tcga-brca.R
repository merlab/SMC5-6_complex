library(TCGAbiolinks)
library(DESeq2)
query <- GDCquery(
    project = "TCGA-BRCA",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
)

# this is to save them all in one place
wd <- getwd()
setwd("~/.cache/TCGAbiolinks/")
GDCdownload(
    query = query
)
eset <- GDCprepare(
    query = query,
    save = TRUE
)
setwd(wd)


eset <- eset[, colData(eset)$sample_type == "Primary Tumor"]
colnames(eset) <- gsub("\\.", "-", colnames(eset))
# remove non-protein coding and duplicated gene
eset <- eset[rowData(eset)$gene_type == "protein_coding", ]
eset <- eset[!duplicated(rowData(eset)$gene_name), ]
rownames(eset) <- rowData(eset)$gene_name
# colnames(eset) <- colData(eset)$patient
clinical <- colData(eset)
studyExpr <- assays(eset)[["unstranded"]]

rvar <- apply(studyExpr, 1, var)
studyExpr <- studyExpr[complete.cases(rvar), ]
rvar <- apply(studyExpr, 1, var)
studyExpr <- studyExpr[rvar != 0, ]
cvar <- apply(studyExpr, 2, var)
studyExpr <- studyExpr[, complete.cases(cvar)]
cvar <- apply(studyExpr, 2, var)
studyExpr <- studyExpr[, cvar != 0]

dds <- DESeqDataSetFromMatrix(countData = studyExpr, colData = clinical, design = ~1)
dds <- estimateSizeFactors(dds)
studyExpr <- counts(dds, normalized = TRUE)

saveRDS(studyExpr, "./reviewer-addressing/tcga/BRCA-expr.rds")
