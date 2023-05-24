# purpose: download maf data of all TCGA studies
source("./R/routine_tasks.R")
source("./R/R_rainclouds.R")
# complexGenes <- c("NSMCE2", "SMC6", "SMC5", "NSMCE1", "NSMCE3", "NSMCE4A", "EID3")
complexGenes <- c("SMC5", "SMC6", "NSMCE1", "NSMCE2", "NSMCE3", "NSMCE4A", "EID3")
folder_check("./reviewer-addressing/tcga/")
library(TCGAbiolinks)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(viridis)
library(rstatix)

# allTypes <- c(
#     "LAML", "ACC", "BLCA", "LGG", "BRCA", "CESC", "CHOL", "LCML", "COAD",
#     "CNTL", "ESCA", "FPPP", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LIHC",
#     "LUAD", "LUSC", "DLBC", "MESO", "MISC", "OV", "PAAD", "PCPG", "PRAD",
#     "READ", "SARC", "SKCM", "STAD", "TGCT", "THYM", "THCA", "UCS", "UCEC",
#     "UVM"
# )
# for (type in allTypes) {
projects <- getGDCprojects()$id
projects <- grep("TCGA-", projects, value = TRUE)


for (project in projects) {
    tryCatch(
        {
            f <- sprintf("./reviewer-addressing/tcga/%s-mut.rds", project)
            if (file.exists(f)) next()
            query <- GDCquery(
                project = project,
                data.category = "Simple Nucleotide Variation",
                access = "open",
                # legacy = FALSE,
                data.type = "Masked Somatic Mutation"
                # workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
            )
            wd <- getwd()
            setwd("~/.cache/TCGAbiolinks/")
            GDCdownload(
                query = query
            )
            dataPrep <- GDCprepare(
                query = query,
                save = TRUE
            )
            setwd(wd)
            saveRDS(dataPrep, f)
        },
        error = function(cond) {
            message(cond)
        }
    )
}
print("done")
