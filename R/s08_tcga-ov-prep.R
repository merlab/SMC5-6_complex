# purpose: makes the dataframe for the alteration analysis of the complexes in the genes
genes <- c("SMC5", "SMC6", "NSMCE1", "NSMCE2", "NSMCE3", "NSMCE4A", "EID3")
genesdet <- paste0(genes, "_det")
dir <- getwd()
source(sprintf("./R/routine_tasks.R"))
enable_parallelization()
library(foreach)
library(data.table)
mutcheck <- function(track, row, df) {
  if (track == "TP53") {
    out <- df[df$track_name == track & df$track_type %in% c("MUTATIONS"), row]
  } else {
    out <- df[df$track_name == track & df$track_type %in% c("MUTATIONS", "STRUCTURAL_VARIANT", "CNA"), row]
  }
  out <- out[out != ""]
  if (length(out) == 0) out <- ""
  out <- paste(as.vector(out), collapse = ";", sep = ";")
  if (out %in% c(";;", ";", "")) out <- ""
  out <- ifelse(out == "", 0, 1)
  return(out)
}

###################
## Configuration ##
###################
df <- as.data.frame(
  data.table::fread("./data/tcga-ov-atlas.tsv", header = TRUE)
)
# we need this for stage data
df2 <- as.data.frame(
  data.table::fread("./data/tcga-ov-firehouse.tsv", header = TRUE)
)


print("df imported")

print(df[, 1])

output <- as.data.frame(matrix(nrow = 0, ncol = (19 + length(genes))))
# output_alt <- as.data.frame(matrix(nrow = ncol(df)-2, ncol = (15 + 2*length(genes))))
colnames(output) <- c(
  "tissue_detailed", "name",
  "age", "OVS", "OVT", "PFS", "PFT", "DFS", "DFT",
  "histology", "grade", "stage",
  "isalt", genes, "TP53", "BRCA1", "BRCA2",
  "rtherapy", "ntherapy", "race"
)
# foreach here combines the results of the comutation
# (outputed as return() function) into a row in the output
# dataframe
output <- foreach(i = 3:ncol(df), .combine = "rbind") %dopar% {
  # for(i in 3:ncol(df)) {
  tissue_detailed <- df[df$track_name == "Cancer Type Detailed", c(i)]
  if (tissue_detailed == "") tissue_detailed <- NA
  name <- colnames(df)[i]
  age <- df[df$track_name == "Diagnosis Age", c(i)]
  if (age == "") age <- NA
  sex <- df[df$track_name == "Sex", c(i)]
  if (sex == "") sex <- NA
  histology <- df[df$track_name == "Cancer Type Detailed", c(i)]
  if (histology == "") histology <- NA
  ntherapy <- df[df$track_name == "Neoadjuvant Therapy Type Administered Prior To Resection Text", c(i)]
  if (ntherapy == "") ntherapy <- NA
  rtherapy <- df[df$track_name == "Radiation Therapy", c(i)]
  if (rtherapy == "") rtherapy <- NA
  race <- df[df$track_name == "Race Category", c(i)]
  if (race == "") race <- NA
  if (name %in% colnames(df2)) {
    stage <- df2[df2$track_name == "Neoplasm American Joint Committee on Cancer Clinical Group Stage", name]
    if (stage == "") stage <- NA
  } else {
    stage <- NA
  }
  grade <- df[df$track_name == "Neoplasm Histologic Grade", c(i)]
  if (grade == "") grade <- NA
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
  if (OVT == "") OVT <- NA
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
  if (PFT == "") PFT <- NA
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
  if (DFT == "") DFT <- NA

  TP53 <- mutcheck("TP53", i, df)
  BRCA1 <- mutcheck("BRCA1", i, df)
  BRCA2 <- mutcheck("BRCA2", i, df)

  # is the gene altered
  isaltgenes <- setNames(rep(0, length(genes)), genes)
  for (k in 1:length(genes)) {
    foo <- df[
      df$track_name == genes[k] & df$track_type %in% c("CNA", "MUTATIONS", "STRUCTURAL_VARIANT"),
      i
    ]
    foo <- foo[foo != ""]
    if (length(foo) == 0) foo <- ""
    foo <- paste(as.vector(foo), collapse = ";", sep = ";")
    if (foo %in% c(";;", ";", "")) foo <- ""
    foo <- gsub(";;", ";", foo)
    isaltgenes[k] <- ifelse(foo == "", 0, 1)
  }
  # TP53 <- df[df$track_name == 'TP53' & df$track_type %in% c('CNA', 'MUTATIONS', 'STRUCTURAL_VARIANT'), i]
  # is at least one of the genes altered?
  isalt <- ifelse(sum(isaltgenes) > 0, 1, 0)
  return(c(
    tissue_detailed = tissue_detailed, name = name,
    age = age, sex = sex, OVS = OVS, OVT = OVT, PFS = PFS, PFT = PFT, DFS = DFS, DFT = DFT,
    histology = histology, grade = grade, stage = stage,
    isalt = isalt, isaltgenes, TP53 = TP53, BRCA1 = BRCA1, BRCA2 = BRCA2,
    ntherapy = ntherapy, rtherapy = rtherapy,
    race = race
  ))
}
# cols <- c('OVT','DFT','PFT','age', genes, 'TP53')
# output[,cols] <- as.numeric(output[,cols])
output <- as.data.frame(output)
rownames(output) <- output$name
output$PFT <- as.numeric(output$PFT)
output$PFS <- as.numeric(output$PFS)
output$age <- as.numeric(output$age)
output$DFS <- as.numeric(output$DFS)
output$DFT <- as.numeric(output$DFT)
output$OVT <- as.numeric(output$OVT)
output$OVS <- as.numeric(output$OVS)
output$race <- as.factor(output$race)
output$isalt <- as.numeric(output$isalt)
output$NSMCE2 <- as.numeric(output$NSMCE2)
print(head(output))
dim(df)
dim(output)
folder_check("./data/tcga-ov/")
saveRDS(output, "./data/tcga-ov/clinical.rds")
print("done")

