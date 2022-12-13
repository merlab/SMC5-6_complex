# code inspired from:
# https://jokergoo.github.io/circlize_book/book/circos-heatmap.html
# https://colorbrewer2.org/
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
      circos.rect(xleft = x[a], ybottom = y - 0.4 
                , xright = x[a] + xlen, ytop = y + 0.4
                , sector.index = cur_sec
                , col = color, border = color 
                , lwd = .03
                )
    }
}
# color for each cancer type
# colors <- rev(c('#e31a1c','#b2df8a','#1f78b4','#33a02c', '#a6cee3','#ff7f00','#cab2d6','#b15928', '#6a3d9a', '#fdbf6f'))
colors <- c('#FEBF70', '#563D9F', '#BD5928', '#CAB9DB', '#FD7D27', '#A3D7E8', '#3CA644', '#257EBB', '#B0DB8E', '#F92F1B')
over60c <- '#F53E98'
under60c <- '#33A679'
malec <- '#70D5E1'
femalec <- '#FAB7CB'
studytcgac <- '#EA6025'
studymetabricc <- '#6673BA'
cbpd <- readRDS("./data/cbioportal_data.rds")
# NOTE: analyze only those with genomic alteration
cbpd <- cbpd[cbpd$isalt == 1, ]
# remove other cancer types
cbpd <- cbpd[cbpd$major != 'Other', ]

# NOTE: we rmeoved other section
cbpd$major <- factor(cbpd$major)
a <- table(cbpd$major)
a <- a[order(a)]
levels <- names(a)
cbpd$major <- factor(cbpd$major, levels = levels)
# make the sectores based on the major cancer types
cbpd <- cbpd[order(cbpd$major), ]
cbpd$sector <- factor(cbpd$major)
# cancer type row
cbpd$one <- 1
# gender
cbpd$male[cbpd$sex %in% c('MALE', 'Male')] <- 1
cbpd$female[cbpd$sex %in% c('FEMALE', 'Female')] <- 1
# age
cbpd[,c('under60', 'over60')] <- 0
cbpd$under60[cbpd$age < 60] <- 1
cbpd$over60[cbpd$age >= 60] <- 1
# study
cbpd$TCGA <- 0
cbpd$TCGA[grep("TCGA",cbpd$study, ignore.case = TRUE)] <-  1
cbpd$metabric <- 0
cbpd$metabric[grep("metabric",cbpd$study, ignore.case = TRUE)] <-  1
cbpd$otherstudy <- 0
cbpd$otherstudy[cbpd$metabric * cbpd$TCGA == 0] <-  1
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
    circos.rect(xleft = 1, ybottom = j - 0.4 
                , xright = x + 1, ytop = j + 0.4 
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
  add_points('male', 3, malec, i, cbpd)
  add_points('female', 3, femalec, i, cbpd)
  add_points('one', 1, colors[i], i, cbpd)
  # add age
  add_points('under60', 4, under60c, i, cbpd)
  add_points('over60', 4, over60c, i, cbpd)
  # add study
  add_points('TCGA' , 2, studytcgac, i, cbpd)
  add_points('metabric', 2, studymetabricc, i, cbpd)
  # add genomic alteration
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
  , legend_gp = gpar(col = c(malec, femalec), fontsize = 18)
  , title_position = "topleft", title = "Sex")
age_leg <- Legend(at = c('<60', '>=60'), type = 'points'
  , legend_gp = gpar(col = c(under60c, over60c), fontsize = 18)
  , title_position = "topleft", title = "Age range")
study_leg <- Legend(at = c('TCGA', "Metabric"), type = 'points'
  , legend_gp = gpar(col = c(studytcgac, studymetabricc), fontsize = 18)
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