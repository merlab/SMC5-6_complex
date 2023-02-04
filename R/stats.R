# purpse: make the statistics generated for the paper
source('./R/routine_tasks.R')
cbpd <- readRDS("./data/cbioportal/formatted.rds")
print(dim(cbpd))
cbpd$age <- as.numeric(cbpd$age)
cbpd$ageCat <- ifelse(cbpd$age <= 60, '<=60', '60+')
cbpd$sex <- as.factor(cbpd$sex)
cbpd$isalt <- as.factor(ifelse(cbpd$isalt == 1, "Altered", "Wild-type"))
cbpd <- cbpd[cbpd$major != 'Other',]

# complex and sex
print('sex stats')
tbl <- table(cbpd$sex, cbpd$isalt)
print(tbl)
print(chisq.test(tbl))

x <- as.data.frame(tbl)

print('agecat stats')
tbl <- table(cbpd$ageCat, cbpd$isalt)
print(tbl)
print(chisq.test(tbl))

print('breast alteration percentage:')
breast <- cbpd[cbpd$major == 'Breast',]
n_breast <- nrow(breast)
n_alt_breast <- sum(breast$isalt == 'Altered')
print(n_alt_breast / n_breast * 100)


print('prostate alteration precentage:')
prostate <- cbpd[cbpd$major == 'Prostate',]
n_prostate <- nrow(prostate)
n_alt_prostate <- sum(prostate$isalt == 'Altered')
print(n_alt_prostate / n_prostate * 100)


print('co-mutation vs all population')
n_all <-  nrow(cbpd)
genes <- c("SMC5", "SMC6", "NSMCE1", "NSMCE2", "NSMCE3", "NSMCE4A", "EID3")
n_muts <- rowSums(cbpd[,genes])
n_comut <- sum(n_muts > 1)
print(n_comut / n_all * 100)
