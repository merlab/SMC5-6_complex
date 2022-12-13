library(writexl)
library(ComplexHeatmap)
library(circlize)
source("./R/routine_tasks.R")
# add plots for genes
add_points <- function(colname,y, color, cur_sec
              , cbpd, cex = 0.001, pch = 0, xlen = 1) {
    cbpd <- cbpd[cbpd$sec == cur_sec, colname]
    cbpd[is.na(cbpd)] <- 0
    x <- as.numeric(cbpd) * rev(seq(from = 1, to = length(cbpd), by = 1))
    for(a in 1:length(x)) {
      if(x[a] == 0) next()
      circos.rect(xleft = x[a], ybottom = y - 0.4 # 0.5
                , xright = x[a] + xlen, ytop = y + 0.4 # 0.5
                , sector.index = cur_sec
                , col = color, border = color 
                , lwd = .03
                )
    }
}
# https://colorbrewer2.org/#
colors <- rev(c('#e31a1c','#b2df8a','#1f78b4','#33a02c', '#a6cee3','#ff7f00','#cab2d6','#b15928', '#6a3d9a', '#fdbf6f'))
cbpd <- readRDS('./results/cbioportal_alt_all.rds')
# NOTE: analyze only those with genomic alteration
cbpd <- cbpd[cbpd$isalt == 1, ]
cbpd <- cbpd[cbpd$major != 'Other', ]

# NOTE: we rmeoved other section
cbpd$major <- factor(cbpd$major)
a <- table(cbpd$major)
a <- a[order(a, decreasing = TRUE)]
a <- a[!names(a) %in% 'Other']
# levels <- c('Other', rev(names(a)))
levels <- c(rev(names(a)))
cbpd$major <- factor(cbpd$major, levels = levels)
cbpd <- cbpd[order(cbpd$major), ]
cbpd$sector <- factor(cbpd$major)
# OVS
cbpd$alive[cbpd$OVT > 60 | cbpd$OVS == 1] <- 1
cbpd$dead[cbpd$OVT <= 60 | cbpd$OVS == 0] <- 1
# gender
cbpd$male[cbpd$sex %in% c('MALE', 'Male')] <- 1
cbpd$female[cbpd$sex %in% c('FEMALE', 'Female')] <- 1
# ploidy
cbpd$ploidy <- ifelse(is.na(cbpd$ploidy), 0, 1)
cbpd$aneuploidyScore <- ifelse(is.na(cbpd$aneuploidyScore), 0, 1)
cbpd$one <- 1
# age
cbpd[,c('to60', 'to200')] <- 0
cbpd$to60[cbpd$age <= 60] <- 1
cbpd$to200[cbpd$age > 60] <- 1
# study
cbpd$TCGA <- 0
cbpd$TCGA[grep("TCGA",cbpd$study, ignore.case = TRUE)] <-  1
cbpd$metabric <- 0
cbpd$metabric[grep("metabric",cbpd$study, ignore.case = TRUE)] <-  1
cbpd$otherstudy <- 0
cbpd$otherstudy[cbpd$metabric * cbpd$TCGA == 0] <-  1
# https://jokergoo.github.io/circlize_book/book/circos-heatmap.html
col_fun1 <- colorRamp2(c(0, 1), c("white", "red"))
sectors <- cbpd$sector

colors <- setNames(colors, unique(sectors))
pdf('./figures/fig2.pdf', height = 8.5, width = 8.5)
n <- length(unique(sectors))
deg <- 60
circos.par(gap.degree = c(rep(1, n-1), deg)
  , track.height = 0.4
  , track.margin = c(0, 0), cell.padding = c(0,0,0,0)
  , start.degree = 90 - deg
  )

circos.initialize(sectors
  , xlim = matrix(c(rep(1, n), (table(sectors) + 1)), ncol = 2))
# NOTE: plots
ymax <- length(genes) + 4
ylabels <- c('Cancer type', 'Study type', 'Sex', 'Age', rev(genes))

circos.track(ylim = c(0.5, ymax + 0.5), bg.border = 'white')

# add data for each sector
for(i in unique(sectors)) {
  # add grey background color
  for(j in 1:ymax) {
    x <- as.numeric(table(sectors)[i])
    circos.rect(xleft = 1, ybottom = j - 0.4 # 0.5
                , xright = x + 1, ytop = j + 0.4 # 0.5
                , sector.index = i
                , col = 'grey90', border = 'grey90'
                , lwd = .03
                )
  }

  # NOTE: order and sort the samples to remove the stray stuff
  cbpdtmp <- cbpd[cbpd$sec == i,]
  for(k in rev(genes)) {
    cbpdtmp <- cbpdtmp[order(cbpdtmp[, k], decreasing = TRUE), ]
  }
  # NOTE: plot the gene lanes
  for(k in length(genes):1) {
    gene <- genes[k]
    add_points(gene, ymax - k + 1, colors[i], i, cbpdtmp)
  }
  # NOTE: plot the metadata lanes
  add_points('male', 3, 'cyan', i, cbpd)
  add_points('female', 3, 'pink', i, cbpd)
  add_points('one', 1, colors[i], i, cbpd)
  # add age
  c <- c('#1b9e77','#d95f02','#7570b3','#e7298a')
  add_points('to60', 4, c[1], i, cbpd)
  add_points('to200', 4, c[4], i, cbpd)
  # add study
  add_points('TCGA' , 2, "#d95f02", i, cbpd)
  add_points('metabric', 2, '#7570b3', i, cbpd)
  # add lines
  for(j in 1:(ymax + 1)) {
    circos.lines(x = c(1, length(sectors[sectors %in% i]))
               , y = c(rep(j - 0.5, 2))
               , col = 'black'
               , sector.index = i
               , lwd = 0.5
                )
  }
}
sex_leg <- Legend(at = c('Male', "Female"), type = 'points'
  , legend_gp = gpar(col = c('cyan', 'pink'), fontsize = 18)
  , title_position = "topleft", title = "Sex")
age_leg <- Legend(at = c('<=60', '>60'), type = 'points'
  , legend_gp = gpar(col = c[c(1,4)], fontsize = 18)
  , title_position = "topleft", title = "Age range")
study_leg <- Legend(at = c('TCGA', "Metabric"), type = 'points'
  , legend_gp = gpar(col = c('#d95f02', '#7570b3'), fontsize = 18)
  , title_position = "topleft", title = "Study type")
lgd_list_vertical <- packLegend(age_leg, sex_leg, study_leg)
draw(lgd_list_vertical, x = unit(6, 'inch'), y = unit(7,'inch'))


# y axis
circos.yaxis(labels = ylabels, at = seq(from = 1, to = ymax, by = 1), side = 'right'
             , col = 'white'
             , labels.cex = 0.7
              )
# add sector names and labels
for(i in unique(sectors)) {
  circos.text(x = length(sectors[sectors %in% i]) / 2
            , y = ymax + 2.1
            , sector.index = i
            , labels = paste0(i, "\n", " (n=", length(sectors[sectors %in% i]), ")")
            , cex = 1
            )
}
circos.clear()
dev.off()
print('done')