
library(ggplot2)
source('./R/routine_tasks.R')
cbpd <- readRDS("./data/cbioportal_data.rds")
print(dim(cbpd))
cbpd$age <- as.numeric(cbpd$age)
#df$sex[df$sex == "MALE"] <- 'Male'
cbpd$sex <- as.factor(cbpd$sex)
cbpd$isalt <- as.factor(ifelse(cbpd$isalt == 1, "Altered", "Wild-type"))
cbpd$ageCat <- ifelse(cbpd$age <= 60, '<=60', '60+')
#print(table(cbpd$sex, cbpd$tissue))
keep <- c()
for(i in unique(cbpd$tissue)) {
        s <- na.omit(cbpd$sex[cbpd$tissue == i])
        if(length(s) == 0) next()
        mr <- sum(s == 'Male') / length(s)
        print(i)
        print(mr)
        if(mr >= 0.25 & mr <= 0.75) keep <- c(keep, i)
}
print(length(keep))
print(dim(cbpd[cbpd$tissue %in% keep,]))
#tissue <- cbpd[cbpd$major != "Other",]
tissue <- cbpd[cbpd$tissue %in% keep,]
#tissue <- cbpd
print(dim(tissue))
df <- tissue[,c('isalt', 'age', 'sex', 'major')]

# complex and sex
print('sex')
print(table(df$major))
print(table(df$sex, useNA = 'a'))
t <- df[#!df$major %in% c("Breast", "Prostate", "Ovarian") 
                ,c('isalt','sex')] # , 'major'
tbl <- table(t)
print(chisq.test(tbl))
print(fisher.test(tbl))
x <- as.data.frame(tbl)
x$sex <- factor(x$sex, levels = c("Male", "Female"))
x$Freq[x$isalt == "Altered"] <- x$Freq[x$isalt == "Altered"] / sum(x$Freq[x$isalt == "Altered"])
x$Freq[x$isalt == "Wild-type"] <- x$Freq[x$isalt == "Wild-type"] / sum(x$Freq[x$isalt == "Wild-type"])
p <-ggplot(x, aes(x = isalt, y = Freq, fill = sex)) +
        #geom_boxplot(outlier.color = NA) +
        geom_bar(position="stack", stat="identity", color = 'black', width = 1) +
        #geom_violin() +
        #stat_compare_means(method = "wilcox.test") + 
        scale_fill_manual(values = c('cyan', 'pink')) +
        xlab('')+
        ylab('Sex proportions') +
        scale_y_continuous(expand = c(0,0), limits = c(0,1), breaks = seq(0,1,0.25)) +
        guides(fill = 'none') +
        coord_flip() +
        theme_classic()
p <- rmbg(p)
pdf('./sex.pdf', height = 2, width = 4)
plot(p)
dev.off()
# chisq$observed
# chisq$expected
# chisq$p.value

# # complex and age
print('age')
# pdf('./age.pdf', width = 3, height = 4)
# p <-ggplot(df, aes(x = isalt, y = age, fill = isalt)) +
#         #geom_boxplot(outlier.color = NA) +
#         geom_violin() +
#         #stat_compare_means(method = "wilcox.test") + 
#         scale_fill_manual(values = c('#DE3B1C','#707176')) +
#         xlab('')+
#         ylab('Age (yr)') +
#         guides(fill = 'none') +
#         theme_classic()
# p <- rmbg(p)
# plot(p)
# dev.off()
#df <- cbpd
#tbl <- table(df[,c('isalt','ageCat')])
# print(chisq.test(tbl))
# print(fisher.test(tbl))
# print(t.test(df$age[df$isalt == "Altered"], df$age[df$isalt == "Wild-type"]))
