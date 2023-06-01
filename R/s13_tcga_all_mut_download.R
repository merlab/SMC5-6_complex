# purpose: download maf data of all TCGA studies
source("./R/routine_tasks.R")
source("./R/R_rainclouds.R")
folder_check("./data/tcga")
complexGenes <- c("SMC5", "SMC6", "NSMCE1", "NSMCE2", "NSMCE3", "NSMCE4A", "EID3")
folder_check("./reviewer-addressing/tcga-mut/")
library(TCGAbiolinks)
library(ggpubr)
library(ggplot2)
library(dplyr)
library(viridis)
library(rstatix)

projects <- getGDCprojects()$id
projects <- grep("TCGA-", projects, value = TRUE)


for (project in projects) {
    tryCatch(
        {
            f <- sprintf("./data/tcga-mut/%s-mut.rds", project)
            if (file.exists(f)) next()
            query <- GDCquery(
                project = project,
                data.category = "Simple Nucleotide Variation",
                access = "open",
                data.type = "Masked Somatic Mutation"
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
