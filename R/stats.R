library(ggplot2)
source('./R/routine_tasks.R')
raw <- readRDS('./results/cbioportal_alt_all_raw.rds')
print(dim(raw))
cbpd <- readRDS('./results/cbioportal_alt_all.rds')
print(dim(cbpd))
tissue <- cbpd[cbpd$major != "Other",]
print(dim(tissue))
stop()
df <- tissue[,c('isalt', 'age', 'sex')]
df$age <- as.numeric(df$age)
df$sex[df$sex == "MALE"] <- 'Male'
df$sex <- as.factor(df$sex)
df$isalt <- as.factor(ifelse(df$isalt == 1, "Altered", "Wild-type"))
df$ageCat <- ifelse(df$age <= 60, '<=60', '60+')

# complex and age
print('age')
pdf('./age.pdf', width = 3, height = 4)
p <-ggplot(df, aes(x = isalt, y = age, fill = isalt)) +
        #geom_boxplot(outlier.color = NA) +
        geom_violin() +
        #stat_compare_means(method = "wilcox.test") + 
        scale_fill_manual(values = c('#DE3B1C','#707176')) +
        xlab('')+
        ylab('Age (yr)') +
        guides(fill = 'none') +
        theme_classic()
p <- rmbg(p)
plot(p)
dev.off()
tbl <- table(df[,c('isalt','ageCat')])
print(chisq.test(tbl))
print(t.test(df$age[df$isalt == "Altered"], df$age[df$isalt == "Wild-type"]))

# complex and sex
print('sex')
tbl <- table(df[,c('isalt','sex')])
print(chisq.test(tbl))

pdf('./sex.pdf', height = 2, width = 4)
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
plot(p)
dev.off()
# chisq$observed
# chisq$expected
# chisq$p.value