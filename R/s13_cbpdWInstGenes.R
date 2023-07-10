# purpose: makes the dataframe for the alteration analysis of the complexes in the genes + instability genes
complexGenes <- c("SMC5", "SMC6", "NSMCE1", "NSMCE2", "NSMCE3", "NSMCE4A", "EID3", "SLF1", "SLF2")
complexGenesDet <- paste0(complexGenes, "_det")
instabilityGenes <- c("TP53", "BRCA1", "BRCA2", "NBN", "TTK", "AURKA", "PLK1",
                      "CHEK2", "CCNE1", "RB1", "RECQL4", "BLM")
instabilityGenesDet <- paste0(instabilityGenes, "_det")
source("./R/routine_tasks.R")
folder_check("./reviewer-addressing/cbioportal")

enable_parallelization(detectCores() %/% 2 + 1)
library(data.table)

###################
## Configuration ##
###################
# what to subset the oncoprint file to


df <- as.data.frame(fread("./data/new_onco/PATIENT_DATA_oncoprint-2023-04-17.tsv",
                                    header = TRUE))
df2 <- as.data.frame(fread("./data/new_onco/instability-genes-2023-04-17.tsv",
                                    header = TRUE))
df3 <- as.data.frame(fread("./data/new_onco/myc-2023-07-03.tsv",
                                    header = TRUE))
# output_alt <- as.data.frame(matrix(nrow = 0, ncol = (27 + 2*length(complexGenes) + 2*length(instabilityGenes))))
output_alt <- as.data.frame(matrix(nrow = 0, ncol = (29 + 2*length(complexGenes) + 2*length(instabilityGenes))))
colnames(output_alt) <- c(
             "name", "study", "tissue", "tissue_detailed",
             "age", "sex", "OVS", "OVT", "PFS", "PFT", "DFS", "DFT",
             "histology", "grade", "stage",
             "isalt", complexGenes, complexGenesDet,
             "instabilityScore", instabilityGenes, instabilityGenesDet,
             "MYC", "MYC_det",
             "aneuploidyScore", "ploidy",
             "ntherapy", "rtherapy", "race",
             "profiledmut", "profiledsv", "profiledcna",
             "noSample"
             )
# foreach here combines the results of the comutation
# (outputed as return() function) into a row in the output
# dataframe
output_alt <- foreach(i = 3:ncol(df), .combine = "rbind") %dopar% {
# for diagnosis purposes
#for (i in 3:ncol(df)) {
    # patient metadata
    name <- colnames(df)[i]
    #
    profiledmut <- df[df$track_name == "Profiled for mutations", c(i)]
    if (profiledmut == "") {
        profiledmut <- "NO"
    }
    #
    profiledsv <- df[df$track_name == "Profiled for structural variants", c(i)]
    if (profiledsv == "") {
        profiledsv <- "NO"
    }
    #
    profiledcna <- df[df$track_name == "Profiled for copy number alterations", c(i)]
    if (profiledcna == "") {
        profiledcna <- "NO"
    }
    #
    study <- df[df$track_name == "Study of origin", c(i)]
    if (study == "") {
        study <- NA
    }
    #
    tissue <- df[df$track_name == "Cancer Type", c(i)]
    if (tissue == "") {
        tissue <- NA
    }
    #
    tissue_detailed <- df[df$track_name == "Cancer Type Detailed", c(i)]
    if (tissue_detailed == "") {
        tissue_detailed <- NA
    }
    #
    noSample <- df[df$track_name == "# Samples per Patient", c(i)]
    if (noSample == "") {
        noSample <- NA
    }
    #
    age <- df[df$track_name == "Diagnosis Age", c(i)]
    if (age == "") {
        age <- NA
    }
    #
    sex <- df[df$track_name == "Sex", c(i)]
    sex <- tolower(sex)
    if (sex == "") {
        sex <- NA
    } else if (sex == "male") {
        sex <- "Male"
    } else if (sex == "female") {
        sex <- "Female"
    }
    #
    histology <- df[df$track_name == "Cancer Type Detailed", c(i)]
    if (histology == "") {
        histology<- NA
    }
    #
    ntherapy <- df[df$track_name == "Neoadjuvant Therapy Type Administered Prior To Resection Text", c(i)]
    if (ntherapy == "") {
        ntherapy <- NA
    }
    #
    rtherapy <- df[df$track_name == "Radiation Therapy", c(i)]
    if (rtherapy == "") {
        rtherapy <- NA
    }
    #
    race <- df[df$track_name == "Race Category", c(i)]
    if (race== "") {
        race <- NA
    }
    #
    aneuploidyScore <- df[df$track_name == "Aneuploidy Score", c(i)]
    if (aneuploidyScore == "") {
        aneuploidyScore <- NA
    }
    #
    ploidy <- df[df$track_name == "Ploidy", c(i)]
    if (ploidy == "") {
        ploidy <- NA
    }
    #
    stage <- df[df$track_name == "Neoplasm American Joint Committee on Cancer Clinical Group Stage",name]
    if (stage == "") {
        stage <- NA
    }
    #
    grade <- df[df$track_name == "Neoplasm Histologic Grade", c(i)]
    if (grade == "") {
        grade<- NA
    }
    # OV
    OVS <- df[df$track_name == "Overall Survival Status", i]
    if (OVS == "0:LIVING") {
        OVS <- 0
    } else if (OVS == "1:DECEASED") {
        OVS <- 1
    } else {
        OVS <- NA
    }
    OVT <- df[df$track_name == "Overall Survival (Months)", c(i)]
    if (OVT == "") {
        OVT <- NA
    }
    # PFS
    PFS <- df[df$track_name == "Progression Free Status", i]
    if (PFS == "0:CENSORED") {
        PFS <- 0
    } else if (PFS == "1:PROGRESSION") {
        PFS <- 1
    } else {
        PFS <- NA
    }
    PFT <- df[df$track_name == "Progress Free Survival (Months)", c(i)]
    if (PFT == "") {
        PFT <- NA
    }
    # DF
    DFS <- df[df$track_name == "Disease Free Status", i]
    if (DFS == "0:DiseaseFree") {
        DFS <- 0
    } else if (DFS == "1:Recurred/Progressed") {
        DFS <- 1
    } else {
        DFS <- NA
    }
    DFT <- df[df$track_name == "Disease Free (Months)", c(i)]
    if (DFT == "") {
        DFT <- NA
    }
    #################
    # complex genes #
    #################
    isaltGenes <- setNames(rep(0, length(complexGenes)), complexGenes)
    isaltGenesDet <- setNames(rep(0, length(complexGenesDet)), complexGenesDet)
    for (k in seq_along(complexGenes)) {
        # this takes all the available alteration info for a gene
        # if anything is altered it would consider the gene altered
        # NOTE: this was modified to remove teh mRNA and protein variations in cbioportal data
        foo <- df[
                df$track_name == complexGenes[k] & df$track_type %in% c("CNA", "MUTATIONS", "STRUCTURAL_VARIANT")
                , i]
        foo <- foo[foo != ""]
        if (length(foo) == 0) {
            foo <- ""
        }
        foo <- paste(as.vector(foo), collapse = ";", sep = ";")
        if (foo %in% c(";;", ";", "")) foo <- ""
        foo <- gsub(";;", ";", foo)
        isaltGenes[k] <- ifelse(foo == "", 0, 1)
        isaltGenesDet[k] <- foo
    }
    # is at least one of the genes altered?
    isalt <- ifelse(sum(isaltGenes) > 0, 1, 0)
    ######################
    # instability genes #
    ######################
    instabilitySum <- setNames(rep(0, length(instabilityGenes)), instabilityGenes)
    instabilitySumDet <- setNames(rep(0, length(instabilityGenesDet)), instabilityGenesDet)
    for (k in seq_along(instabilityGenes)) {
        # this takes all the available alteration info for a gene
        # if anything is altered it would consider the gene altered
        # NOTE: this was modified to remove teh mRNA and protein variations in cbioportal data
        foo <- df2[
                df2$track_name == instabilityGenes[k] & df2$track_type %in% c("CNA", "MUTATIONS", "STRUCTURAL_VARIANT")
                , i]

        foo <- foo[foo != ""]
        if (length(foo) == 0) {
            foo <- ""
        }
        foo <- paste(as.vector(foo), collapse = ";", sep = ";")
        # foo <- gsub("splice", "", foo)
        foo <- gsub(";;", ";", foo)
        if (foo %in% c(";;", ";", "")) foo <- ""
        foo <- gsub(";;", ";", foo)
        if (foo %in% c(";;", ";", "")) foo <- ""
        instabilitySum[k] <- ifelse(foo == "", 0, 1)
        instabilitySumDet[k] <- foo
    }
    instabilityScore <- sum(instabilitySum)

    ### MYC gene

        # this takes all the available alteration info for a gene
        # if anything is altered it would consider the gene altered
        # NOTE: this was modified to remove teh mRNA and protein variations in cbioportal data
    MYC <- 0
    MYC_det <- ""
    foo <- df3[
            df3$track_name == "MYC" & df3$track_type %in% c("CNA", "MUTATIONS", "STRUCTURAL_VARIANT")
            , i]
    foo <- foo[foo != ""]
    if (length(foo) == 0) {
        foo <- ""
    }
    foo <- paste(as.vector(foo), collapse = ";", sep = ";")
    if (foo %in% c(";;", ";", "")) foo <- ""
    foo <- gsub(";;", ";", foo)
    MYC <- ifelse(foo == "", 0, 1)
    MYC_det <- foo


    # progress bar
    if (i %% 1000 == 0) {
        print(paste(i / ncol(df) * 100, "done"))
    }

    # this is the vector that is returned which consists of
    # alteration status and patient clincal data
    # the rbind command adds these to a dataframe
    return(c(name = name,
             study = study,
             tissue = tissue,
             tissue_detailed = tissue_detailed,
             age = age,
             sex = sex,
             OVS = OVS,
             OVT = OVT,
             PFS = PFS,
             PFT = PFT,
             DFS = DFS,
             DFT = DFT,
             histology = histology,
             grade = grade,
             stage = stage,
             # complex genes
             isalt = isalt,
             isaltGenes,
             isaltGenesDet,
             # instability markers
             instabilityScore = instabilityScore,
             instabilitySum,
             instabilitySumDet,
             MYC = MYC,
             MYC_det = MYC_det,
             # rest of stuff
             aneuploidyScore = aneuploidyScore,
             ploidy = ploidy,
             ntherapy = ntherapy,
             rtherapy = rtherapy,
             race = race,
             profiledmut = profiledmut,
             profiledsv = profiledsv,
             profiledcna = profiledcna,
             noSample = noSample
            )
    )
}

output_alt <- as.data.frame(output_alt)

# saveRDS(output_alt, "./reviewer-addressing/cbioportal/raw.rds")
# print("done")
# 
# 
# df <- readRDS("./reviewer-addressing/cbioportal/raw.rds")
df <- output_alt
repeated_samples <- df$name[duplicated(df$name)]
print(head(repeated_samples))
print(length(repeated_samples))
print(dim(df))
print("checking for non-equivalent replicates")
for (i in repeated_samples) {
    temp_df <- df[df$name == i, ]
    if (nrow(temp_df) == 1) next()
    torm <- c()
    for (j in 2:nrow(temp_df)) {
        j_in_jmin1 <- sum(temp_df[j, ] %in% temp_df[j - 1, ])
        jmin1_in_j <- sum(temp_df[j - 1, ] %in% temp_df[j, ])
        if (j_in_jmin1 > jmin1_in_j) {
            torm <- c(torm, j)
        } else {
            torm <- c(torm, j - 1)
        }
    }
    temp_df <- temp_df[-torm, ]
    if (nrow(temp_df) > 1) stop()
    df <- df[df$name != i, ]
    df <- rbind(df, temp_df[1, ])
}
print("deduplication done")
print(dim(df))
rownames(df) <- df$name
#
# purpose: format the data for analysis
library(writexl)
folder_check("./results/")
cbpd <- df
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
for (i in c("isalt", instabilityGenes, complexGenes, "MYC")) {
    cbpd[, i] <- as.numeric(cbpd[, i])
}
saveRDS(cbpd, "./data/cbioportal/cbpdDataWInst.rds")
print("done")
