# purpose: format the data for analysis
library(writexl)
complexGenes <- c("SMC5", "SMC6", "NSMCE1", "NSMCE2", "NSMCE3", "NSMCE4A", "EID3", "SLF1", "SLF2")
instabilityGenes <- c(
    "TP53", "BRCA1", "BRCA2", "NBN", "TTK", "AURKA", "PLK1",
    "CHEK2", "CCNE1", "RB1", "RECQL4", "BLM"
)
source("./R/routine_tasks.R")
folder_check("./results/")
cbpd <- readRDS("./reviewer-addressing/cbioportal/deduplicated.rds")
# only keep samples that have been profiled for mutation in the complex
print(dim(cbpd))
cbpd <- cbpd[cbpd$profiledmut == "Yes" | cbpd$profiledsv == "Yes" | cbpd$profiledcna == "Yes", ]
print(dim(cbpd))
# only keep studies with +20 cohort size
studies <- table(cbpd$study)
for (i in names(studies[studies < 20])) {
    cbpd <- cbpd[cbpd$study != i, ]
}
print(dim(cbpd))
# remove pediatric studies
cbpd <- cbpd[-grep("Pediatric", cbpd$study, ignore.case = TRUE), ]
print(dim(cbpd))

# find the major cancer type of the samples
tissue_types <- c(
    "Breast", "Prostate", "Melanoma", "Ovarian", "Endometrial", "Lung", "Pancreatic", "Bladder",
    "Hepatobiliary", "Esophagogastric"
)
cbpd$major <- "Other"
for (i in tissue_types) {
    cbpd$major[grep(i, cbpd$tissue, ignore.case = TRUE)] <- i
    cbpd$major[grep(i, cbpd$tissue_detailed, ignore.case = TRUE)] <- i
}

# convert some important columns to numeric
cbpd$PFT <- as.numeric(cbpd$PFT)
cbpd$PFS <- as.numeric(cbpd$PFS)
cbpd$age <- as.numeric(cbpd$age)
cbpd$DFS <- as.numeric(cbpd$DFS)
cbpd$DFT <- as.numeric(cbpd$DFT)
cbpd$OVT <- as.numeric(cbpd$OVT)
cbpd$OVS <- as.numeric(cbpd$OVS)
cbpd$race <- as.factor(cbpd$race)
cbpd$major <- factor(cbpd$major)
cbpd$age <- as.numeric(cbpd$age)
for (i in c("isalt", instabilityGenes, complexGenes)) {
    cbpd[, i] <- as.numeric(cbpd[, i])
}
# save the formatted results
print(table(cbpd$major))
saveRDS(cbpd, "./reviewer-addressing/cbioportal/formatted.rds")
# part 2 remove others
cbpd <- cbpd[!cbpd$major %in% "Other", ]
saveRDS(cbpd, "./reviewer-addressing/cbioportal/format_exOther.rds")
cbpd <- cbpd[cbpd$isalt == 1, ]
cbpd$major <- factor(cbpd$major)
tissuef <- table(cbpd$major)
tissuef <- tissuef[order(tissuef, decreasing = TRUE)]
cbpd <- cbpd[cbpd$major %in% names(tissuef)[1:10], ]
saveRDS(cbpd, "./reviewer-addressing/cbioportal/top10_mut.rds")
print("done")
