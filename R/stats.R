source('./R/routine_tasks.R')
cbpd <- readRDS("./data/cbioportal/formatted.rds")
print(dim(cbpd))
cbpd$age <- as.numeric(cbpd$age)
#df$sex[df$sex == "MALE"] <- 'Male'
cbpd$sex <- as.factor(cbpd$sex)
cbpd$isalt <- as.factor(ifelse(cbpd$isalt == 1, "Altered", "Wild-type"))
cbpd$ageCat <- ifelse(cbpd$age <= 60, '<=60', '60+')

keep <- c()
for(i in unique(cbpd$major)) {
        s <- na.omit(cbpd$sex[cbpd$major == i])
        if(length(s) == 0) next()
        mr <- sum(s == 'Male') / length(s)
        print(i)
        print(mr)
        if(mr >= 0.25 & mr <= 0.75) keep <- c(keep, i)
}
print(length(keep))
print(dim(cbpd[cbpd$major %in% keep,]))
tissue <- cbpd[cbpd$major %in% keep,]
print(dim(tissue))

# complex and sex
print('sex')
print(table(tissue$sex, tissue$major, useNA = 'a'))
tbl <- table(tissue$sex, tissue$isalt)
print(chisq.test(tbl))
print(fisher.test(tbl))
x <- as.data.frame(tbl)
colnames(x) <- c('sex', 'isalt', 'Freq')
x$sex <- factor(x$sex, levels = c("Male", "Female"))
x$Freq[x$isalt == "Altered"] <- x$Freq[x$isalt == "Altered"] / sum(x$Freq[x$isalt == "Altered"])
x$Freq[x$isalt == "Wild-type"] <- x$Freq[x$isalt == "Wild-type"] / sum(x$Freq[x$isalt == "Wild-type"])

df <- tissue
tbl <- table(df[,c('isalt','ageCat')])
print(chisq.test(tbl))
print(fisher.test(tbl))
print(wilcox.test(df$age[df$isalt == "Altered"], df$age[df$isalt == "Wild-type"]))
print(t.test(df$age[df$isalt == "Altered"], df$age[df$isalt == "Wild-type"]))
