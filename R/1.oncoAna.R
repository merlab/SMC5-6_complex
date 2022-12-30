# purpose: makes the dataframe for the alteration analysis of the complexes in the genes
genes <- c("SMC5", "SMC6", "NSMCE1", "NSMCE2", "NSMCE3", "NSMCE4A", "EID3")
genesdet <- paste0(genes, '_det')
source("./R/routine_tasks.R")
folder_check('./data/cbioportal')
enable_parallelization(detectCores() %/% 2)
library(data.table)

###################
## Configuration ##
###################
# what to subset the oncoprint file to


df <- as.data.frame(data.table::fread("./data/PATIENT_DATA_oncoprint.tsv",header = TRUE))

output_alt <- as.data.frame(matrix(nrow = 0, ncol = (26 + 2*length(genes))))
colnames(output_alt) <- c(
             'name', 'study', 'tissue', 'tissue_detailed', 
             'age', 'sex', 'OVS', 'OVT', 'PFS', 'PFT', 'DFS', 'DFT',
             'histology', 'grade', 'stage',
             'isalt', genes, genesdet, 'TP53',
             'aneuploidyScore', 'ploidy',
             'ntherapy', 'rtherapy', 'race',
             'profiledmut', 'profiledsv', 'profiledcna',
             'noSample'
             )
# foreach here combines the results of the comutation 
# (outputed as return() function) into a row in the output
# dataframe
output_alt <- foreach (i = 3:ncol(df), .combine = "rbind") %dopar% {
# for diagnosis purposes
#for (i in 3:ncol(df)) {
    # patient metadata
    name <- colnames(df)[i]
    profiledmut <- df[df$track_name == "Profiled for mutations", c(i)]; if (profiledmut == "") profiledmut <- 'NO'
    profiledsv <- df[df$track_name == "Profiled for structural variants", c(i)]; if (profiledsv == "") profiledsv <- 'NO'
    profiledcna <- df[df$track_name == "Profiled for copy number alterations", c(i)]; if (profiledcna == "") profiledcna <- 'NO'
    study <- df[df$track_name == "Study of origin", c(i)]; if (study == "") study <- NA
    tissue <- df[df$track_name == "Cancer Type", c(i)]; if (tissue == "") tissue <- NA
    tissue_detailed <- df[df$track_name == "Cancer Type Detailed", c(i)]; if (tissue_detailed == "") tissue_detailed <- NA
    noSample <- df[df$track_name == "# Samples per Patient", c(i)]; if (noSample == "") noSample <- NA
    age <- df[df$track_name == "Diagnosis Age", c(i)]; if (age == "") age <- NA
    sex <- df[df$track_name == "Sex",c(i)]
    sex <- tolower(sex)
    if (sex == "") sex <- NA
    else if(sex == 'male') sex <- "Male"
    else if(sex == 'female') sex <- "Female"
    histology <- df[df$track_name == "Cancer Type Detailed",c(i)]; if (histology == "") histology<- NA
    ntherapy <- df[df$track_name == "Neoadjuvant Therapy Type Administered Prior To Resection Text",c(i)]; if (ntherapy == "") ntherapy <- NA
    rtherapy <- df[df$track_name == "Radiation Therapy",c(i)]; if (rtherapy == "") rtherapy <- NA
    race <- df[df$track_name == "Race Category",c(i)]; if (race== "") race <- NA
    aneuploidyScore <- df[df$track_name == "Aneuploidy Score",c(i)]; if (aneuploidyScore == "") aneuploidyScore <- NA
    ploidy <- df[df$track_name == "Ploidy",c(i)]; if (ploidy == "") ploidy <- NA
    stage <- df[df$track_name == "Neoplasm American Joint Committee on Cancer Clinical Group Stage",name]; if (stage == "") stage <- NA
    grade <- df[df$track_name == "Neoplasm Histologic Grade",c(i)]; if (grade == "") grade<- NA

    # OV
    OVS <- df[df$track_name == "Overall Survival Status",i]
    if (OVS == "0:LIVING") {
        OVS <- 0
    } else if (OVS == "1:DECEASED") {
        OVS <- 1
    } else OVS <- NA 
    OVT <- df[df$track_name == "Overall Survival (Months)",c(i)]
    if (OVT == "") OVT <- NA
    # PFS
    PFS <- df[df$track_name == "Progression Free Status",i]
    if (PFS == "0:CENSORED") {
        PFS <- 0
    } else if (PFS == "1:PROGRESSION") {
        PFS <- 1
    } else PFS <- NA 
    PFT <- df[df$track_name == "Progress Free Survival (Months)",c(i)]
    if (PFT == "") PFT <- NA
    # DF
    DFS <- df[df$track_name == "Disease Free Status",i]
    if (DFS == "0:DiseaseFree") {
        DFS <- 0
    } else if (DFS == "1:Recurred/Progressed") {
        DFS <- 1
    } else DFS <- NA 
    DFT <- df[df$track_name == "Disease Free (Months)",c(i)]
    if (DFT == "") DFT <- NA

    # is the gene altered
    isaltgenes <- setNames(rep(0, length(genes)), genes)
    isaltgenesdet <- setNames(rep(0, length(genesdet)), genesdet)
    for (k in 1:length(genes)) {
        # this takes all the available alteration info for a gene 
        # if anything is altered it would consider the gene altered
        # NOTE: this was modified to remove teh mRNA and protein variations in cbioportal data
        foo <- df[df$track_name == genes[k] & df$track_type %in% c('CNA', 'MUTATIONS', 'STRUCTURAL_VARIANT')
                , i]
        foo <- foo[foo != '']
        if(length(foo) == 0) foo <- ''
        foo <- paste(as.vector(foo), collapse = ';', sep = ';')
        if(foo %in% c(';;', ';', '') ) foo <- ''
        foo <- gsub(';;', ';', foo)
        isaltgenes[k] <- ifelse(foo == "", 0, 1)
        isaltgenesdet[k] <- foo
    }

    # is at least one of the genes altered?
    isalt <- ifelse(sum(isaltgenes) > 0, 1, 0)

    TP53 <- df[df$track_name == 'TP53' & df$track_type %in% c('MUTATIONS', 'STRUCTURAL_VARIANT'), i]
    TP53 <- TP53[TP53 != '']
    if(length(TP53) == 0) TP53 <- ''
    TP53 <- paste(as.vector(TP53), collapse = ';', sep = ';')
    if(TP53 %in% c(';;', ';', '') ) TP53 <- ''
    TP53 <- ifelse(TP53 == "", 0, 1)


    # progress bar
    if (i %% 1000 == 0) print(paste(i / ncol(df) * 100, "done"))

    # this is the vector that is returned which consists of 
    # alteration status and patient clincal data
    # the rbind command adds these to a dataframe
    return(c(name = name, study = study, tissue = tissue, tissue_detailed = tissue_detailed, 
             age = age, sex = sex, OVS = OVS, OVT = OVT, PFS = PFS, PFT = PFT, DFS = DFS, DFT = DFT,
             histology = histology, grade = grade, stage = stage,
             isalt = isalt, isaltgenes, isaltgenesdet, TP53 = TP53,
             aneuploidyScore = aneuploidyScore, ploidy = ploidy,
             ntherapy = ntherapy, rtherapy = rtherapy, race = race,
             profiledmut = profiledmut, profiledsv = profiledsv, profiledcna = profiledcna,
             noSample = noSample
            )
    )
}

output_alt <- as.data.frame(output_alt)

saveRDS(output_alt, "./data/cbioportal/raw.rds")
print("done")