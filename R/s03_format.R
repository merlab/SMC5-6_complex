# purpose: format the data for analysis
library(writexl)
source("./R/routine_tasks.R")
cbpd <- readRDS("./data/cbioportal/deduplicated.rds")
# only keep samples that have been profiled for mutation in the complex
print(dim(cbpd))
cbpd <- cbpd[cbpd$profiledmut == 'Yes' | cbpd$profiledsv == 'Yes' | cbpd$profiledcna == 'Yes', ]
print(dim(cbpd))
# only keep studies with +20 cohort size
studies <- table(cbpd$study)
for(i in names(studies[studies < 20])) {
  cbpd <- cbpd[cbpd$study != i, ]
}
print(dim(cbpd))
# remove pediatric studies
cbpd <- cbpd[-grep("Pediatric", cbpd$study, ignore.case = TRUE), ]
print(dim(cbpd))

# find the major cancer type of the samples
tissue_types <- c("Breast", "Prostate", "Melanoma", "Ovarian", "Endometrial", "Lung", 'Pancreatic', 'Bladder'
  , 'Hepatobiliary', 'Esophagogastric')
cbpd$major <- 'Other'
for(i in tissue_types) {
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
for(i in c('isalt', genes)) {
    cbpd[,i] <- as.numeric(cbpd[,i])
}
# save the formatted results
print(table(cbpd$major))
saveRDS(cbpd, "./data/cbioportal/formatted.rds")
write_xlsx(cbpd, './results/cbioportal_formatted.xlsx')
# part 2 remove others
cbpd <- cbpd[!cbpd$major %in% "Other",]
saveRDS(cbpd, "./data/cbioportal/format_exOther.rds")
# get specific columns we want
df <- cbpd[,c('name','tissue','study','major', grep('_det$',colnames(cbpd), value = TRUE))]
df$major <- factor(df$major
  , levels = c("Breast", "Prostate", "Lung", "Melanoma"
             , "Ovarian", "Esophagogastric", "Hepatobiliary"
             , "Endometrial", "Pancreatic", "Bladder"))
df <- df[order(df$major),]
df$major <- NULL
colnames(df)[3] <- 'Study name'
colnames(df) <- gsub('_det', '', colnames(df))
write_xlsx(df, './results/Supplementary-Data-1_cbioportal_oncoprint.xlsx')
# part 3 keep only mut
cbpd <- cbpd[cbpd$isalt == 1,]
cbpd$major <- factor(cbpd$major)
tissuef <- table(cbpd$major)
tissuef <- tissuef[order(tissuef, decreasing = TRUE)]
cbpd <- cbpd[cbpd$major %in% names(tissuef)[1:10],]
saveRDS(cbpd, "./data/cbioportal/top10_mut.rds")
print('done')
