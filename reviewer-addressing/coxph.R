library(ggplot2)
library(survminer)
library(survival)
cbpd <- readRDS("./data/cbioportal/format_exOther.rds")
cbpd$instabilityScore <- as.numeric(cbpd$instabilityScore)
cbpd$instabilityAbove1 <- cbpd$instabilityScore > 1

#  [1] "name"             "study"
#   "tissue"
#  [4] "tissue_detailed"  "age"
#   "sex"
#  [7] "OVS"              "OVT"
#   "PFS"
# [10] "PFT"              "DFS"
#   "DFT"
# [13] "histology"        "grade"
#   "stage"
# [16] "isalt"            "SMC5"
#   "SMC6"
# [19] "NSMCE1"           "NSMCE2"
#   "NSMCE3"
# [22] "NSMCE4A"          "EID3"
#   "SLF1"
# [25] "SLF2"             "SMC5_det"
#   "SMC6_det"
# [28] "NSMCE1_det"       "NSMCE2_det"
#   "NSMCE3_det"
# [31] "NSMCE4A_det"      "EID3_det"
#   "SLF1_det"
# [34] "SLF2_det"         "instabilityScore
# " "TP53"
# [37] "BRCA1"            "BRCA2"
#   "NBN"
# [40] "TTK"              "AURKA"
#   "PLK1"
# [43] "CHEK2"            "CCNE1"
#   "RB1"
# [46] "RECQL4"           "BLM"
#   "TP53_det"
# [49] "BRCA1_det"        "BRCA2_det"
#   "NBN_det"
# [52] "TTK_det"          "AURKA_det"
#   "PLK1_det"
# [55] "CHEK2_det"        "CCNE1_det"
#   "RB1_det"
# [58] "RECQL4_det"       "BLM_det"
#   "aneuploidyScore"
# [61] "ploidy"           "ntherapy"
#   "rtherapy"
# [64] "race"             "profiledmut"
#   "profiledsv"
# [67] "profiledcna"      "noSample"
#   "major"

pdf("coxana.pdf", height = 8, width = 8)
#### breast analysis
breast <- cbpd[cbpd$major == "Breast" & cbpd$study == "Breast Cancer (METABRIC, Nature 2012 & Nat Commun 2016)", ]
cox <- coxph(Surv(OVT, OVS) ~ grade + histology + instabilityScore + isalt, data = breast)
plot(ggforest(cox))
###
prostate <- cbpd[cbpd$major == "Prostate" &
    cbpd$study %in% c(
        "Pan-cancer analysis of whole genomes (ICGC/TCGA, Nature 2020)",
        "Prostate Adenocarcinoma (TCGA, PanCancer Atlas)"
    ), ] #  & cbpd$study == "", ]
cox <- coxph(Surv(OVT, OVS) ~ age + instabilityScore + isalt, data = prostate)
plot(ggforest(cox))
####
dev.off()
print("done")
