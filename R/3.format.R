# purpose: format the data for analysis
library(writexl)
source("./R/routine_tasks.R")
cbpd <- readRDS("./data/cbioportal/deduplicated.rds")
# only keep samples that have been profiled for mutation in the complex
print(dim(cbpd))
cbpd <- cbpd[cbpd$profiledmut == 'Yes' | cbpd$profiledsv == 'Yes' | cbpd$profiledcna == 'Yes', ]
print(dim(cbpd))
# print(table(cbpd$tissue))
# only keep studies with +20 cohort size
studies <- table(cbpd$study)
for(i in names(studies[studies < 20])) {
  cbpd <- cbpd[cbpd$study != i, ]
}
print(dim(cbpd))
# remove pediatric studies
cbpd <- cbpd[-grep("Pediatric", cbpd$study, ignore.case = TRUE), ]
print(dim(cbpd))
cbpd <- cbpd[-grep("embryonal", cbpd$tissue, ignore.case = TRUE), ]
#cbpd <- cbpd[-grep("embryo", cbpd$tissue, ignore.case = TRUE), ]
#cbpd <- cbpd[-grep("cell line", cbpd$tissue, ignore.case = TRUE), ]
#cbpd <- cbpd[-grep("cell-line", cbpd$tissue, ignore.case = TRUE), ]
print(dim(cbpd))

# find the major cancer type of the samples
tissue_types <- list(
    "Adrenocortical" = c('adrenocortical', 'adrenal')
  , "Ampullary" = c('ampullary', 'ampulla')
  , "Biliary Tract" = c("gallbladder", 'bile duct', 'cholangiocarcinoma', 'biliary')
  , 'Bladder' = c('bladder', 'urothelial')
  , "Bone" = c('bone', 'ewing', 'osteo')
  , "Bowel" = c('colorectal', 'rectal', 'colon', 'bowel')
  , "Breast" = c('breast')
  # NOTE: we categorized CNS and PNS under one term
  , "CNS/PNS" = c(
                  # present in both
                  'nerve', 'nervous', 'glioma'
                  # embryonic
                  # related to cns
                  , 'glioblastoma', 'medulloblastoma', 'astrocytoma'
                  , 'meningioma', 'CNS', 'paraganglioma', 'brain'
                  # related to pns
                  , 'PNS', 'spine','spinal', 'peripheral nerve'
                  , 'neuroblastoma', 'peripheral nervous system'
                  )
  , "Cervical" = c('cervix', 'cervical')
  , 'Esophagogastric/Gastrointestinal' = c('esophagogastric', 'esophageal', 'gastric', 'stomach', 'esophagus', 'gastrointestinal')
  , "Eye" = c('eye', 'optical', 'uvea', 'uveal', 'retinoblastoma', 'retina')
  , "Head/Neck/Salivary" = c('head', 'neck', 'oral', 'nasopharyngeal'
                          , 'nasopharynx', 'adenoid', 'nasal'
                          , 'salivary', 'saliva'
                          )
  , "Kidney" = c('kidney', 'renal', 'nephron', 'rhabdoid', 'wilms')
  , "Liver" = c('hepatocellular', 'liver', 'hepato', 'hepatic', 'heptocellular')
  , "Lung" = c('lung', 'thoracic')
  # NOTE: we mixed lympoid and myloid cacner into blood
  , "Blood" = c(
                  # general
                  'leukemia', 'blood'
                  # lymphatic
                  , 'lymphoid', 'thymus', 'lymphoblastic', "b cell", 'b-cell'
                  , 'lymphoma', 'lymphocytic'
                  # myloid  
                  , 'myeloid', 'myelodysplastic', 'myelodysplasia', 'histiocytosis'
                  , 'myelodysplastic', 'myeloproliferative', 'thrombocythemia'
                  )
  , "Ovarian" = c('ovary', 'ovarian')
  , 'Pancreatic' = c('pancreatic', 'pancreas', 'acinar')
  , "Pleura" = c('pleura', 'pleural','mesothelioma')
  , "Prostate" = c('prostate')
  , "Melanoma/Skin" = c('skin', 'melanoma', 'acral')
  , "Soft tissue" = c('soft tissue','sarcoma', 'angiosarcoma', 'sarcomal', 'rhabdomyosarcoma')
  , "Testis" = c('testicular', 'germ cell', 'testis')
  , "Thymus" = c('thymus', 'thymic', 'thymoma')
  , "Thyroid" = c('thyroid')
  , "Uterus/Endometrial" = c('endometrial', 'uterine', 'uterus')
  , "Vulva/Vagina" = c('vulva', 'vagina', 'vaginal')
  , "Complex/Unknow" = c("unknown", 'mixed')
  )
cbpd$major <- 'Other'
for(i in 1:length(tissue_types)) {
  for(j in tissue_types[[i]]) {
    cbpd$major[grep(j, cbpd$tissue, ignore.case = TRUE)] <- names(tissue_types)[i]
  }
  #cbpd$major[grep(i, cbpd$tissue, ignore.case = TRUE)] <- names(i)
    #cbpd$major[grep(i, cbpd$tissue_detailed, ignore.case = TRUE)] <- i
}
# x <- (table(cbpd$major))
# x <- x[order(x, decreasing = TRUE)]
# print(x)
# y <- (table(cbpd$isalt, cbpd$major))
# print(y)

# convert some important values to numberic
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
cbpd$sex <- as.factor(cbpd$sex)
cbpd$ageCat <- ifelse(cbpd$age <= 60, '<=60', '60+')
for(i in c('isalt', genes)) {
    cbpd[,i] <- as.numeric(cbpd[,i])
}
print(table(cbpd$major))
saveRDS(cbpd, "./data/cbioportal/formatted.rds")
write_xlsx(cbpd, './results/cbioportal_formatted.xlsx')
cbpd <- cbpd[!cbpd$major %in% "Other",]
saveRDS(cbpd, "./data/cbioportal/format_exOther.rds")
cbpd <- cbpd[cbpd$isalt == 1,]
cbpd$major <- factor(cbpd$major)
tissuef <- table(cbpd$major)
tissuef <- tissuef[order(tissuef, decreasing = TRUE)]
cbpd <- cbpd[cbpd$major %in% names(tissuef)[1:10],]
saveRDS(cbpd, "./data/cbioportal/top10_mut.rds")
print('done')
# cbpd <- cbpd[cbpd$major %in% c("Breast", "Prostate", "Melanoma"
#                               , "Ovarian", "Endometrial", "Lung", 'Pancreatic'
#                               , 'Bladder', 'Hepatobiliary', 'Esophagogastric')
#             ,]