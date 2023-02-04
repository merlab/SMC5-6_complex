row <- 'isalt'
# row <- 'NSMCE2'
# NOTE: guide used:
# https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html
# https://programmer.ink/think/rna-4-edger-based-on-tcga-differential-expression-in-sci-articles.html
library(edgeR)
library(Biobase)
library(writexl)
# library(SummarizedExperiment)

TCGA <- readRDS('./data/tcga-prad/rnaseq.rds')
#dds <- edgeR::calcNormFactors(dds)
# remove genes with NA values

cbpd <- readRDS('./data/cbioportal/formatted.rds')
rownames(cbpd) <- cbpd$name
sharedPatients <- intersect(cbpd$name, colnames(TCGA))
cbpd <- cbpd[sharedPatients, ]
countData <- TCGA[, sharedPatients]

v <- as.factor(ifelse(cbpd[, row] == '1', 'Yes', 'No'))
v <- factor(v, levels = c("No", 'Yes'))
alterationData <- setNames(v, cbpd$name)

d <- DGEList(counts = countData, group = alterationData)
print(dim(d))

d.full <- d

keep <- rowSums(cpm(d)>100) >= 2
d <- d[keep,]
print(dim(d))

d <- calcNormFactors(d)

# plotMDS(d, method="bcv", col=as.numeric(d$samples$group))
# legend("bottomleft", as.character(unique(d$samples$group)), col=1:3, pch=20)

d1 <- estimateCommonDisp(d, verbose = TRUE)
d1 <- estimateTagwiseDisp(d1, verbose = TRUE)
plotBCV(d1)



design.mat <- model.matrix(~ 0 + d$samples$group)
colnames(design.mat) <- levels(d$samples$group)
d2 <- estimateGLMCommonDisp(d,design.mat)
# You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
# The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
# NOTE: try each model to see which one fits better to your model!
d2 <- estimateGLMTrendedDisp(d2,design.mat, method = "auto")
plotBCV(d2)
d2 <- estimateGLMTagwiseDisp(d2,design.mat)
plotBCV(d2)

et12 <- exactTest(d1, pair=c(1,2)) # compare groups 1 and 2
print(et12)
print(topTags(et12, n=10))
de1 <- decideTestsDGE(et12, adjust.method="BH", p.value=0.05)
summary(de1)

fit <- glmFit(d2, design.mat)
# compare (group 1 - group 2) to 0:
# this is equivalent to comparing group 1 to group 2
lrt12 <- glmLRT(fit, contrast=c(-1,1))

model1 <- (topTags(et12, n = Inf))
model2 <- (topTags(lrt12, n = Inf))

# NOTE: i am not fully sure which model is better
print(head(model1))
print(head(model2))

df <- model2[[1]]


saveRDS(df, './data/tcga-prad/dgea-edgeR.rds')


print('done')
