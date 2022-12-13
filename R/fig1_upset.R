# https://krassowski.github.io/complex-upset/articles/Examples_R.html
source("./R/routine_tasks.R")
library(ComplexHeatmap)
cbpd <- readRDS('./results/cbioportal_alt_all.rds')
cbpd <- cbpd[cbpd$isalt == 1, ]
cbpd <- cbpd[cbpd$major != 'Other',]
cbpd$major <- factor(cbpd$major)

makeUpset <- function(sdf, title, cut = 2) {
    colors <- c('#66c2a5','#fc8d62','#8da0cb')
    m <- make_comb_mat(sdf, mode = 'distinct')
    m <- m[comb_size(m) >= cut]
    print(m)
    # data at the top of teh upset plot
    breaks <- c(1, 3, 10, 20, 50, 100, 250, 500, 1000, 2000)
    top_ha <- HeatmapAnnotation(
        'Number of \n patients' = anno_barplot(
                log10(comb_size(m))
                , bar_width = 0.8
                , gp = gpar(fill = 'black')#colors[1])#"black")
                , border = FALSE
                , height = unit(5, "cm")
                , ylim = c(1,log10(2000))
                , axis_param = list(at = log10(breaks), labels = breaks)
                                )
        , annotation_name_side = "left"
        , annotation_name_rot = 0
    )
    # data at the bottom of the upset plot
    right_ha = rowAnnotation(
        'Set size' = anno_barplot(
                log10(set_size(m))
                , bar_width = 0.8
                , gp = gpar(fill = 'black')#colors[3])#"black")
                , border = FALSE
                , width = unit(5, "cm")
                , ylim = c(1, log10(2000))
                , axis_param = (list(at = log10(breaks), labels = breaks))
                )
        , annotation_name_side = "bottom"
        )
    UpSet(m, comb_order = rev(order(comb_size(m)))
        , top_annotation = top_ha
        , row_names_side = 'left'
        , left_annotation = right_ha
        )


}

pdf('./figures/fig1_upset.pdf', onefile = TRUE, width = 8, height = 5)
plot(makeUpset(cbpd[, genes], 'All', 10))
dev.off()

print('done')