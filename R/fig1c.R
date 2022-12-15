# https://krassowski.github.io/complex-upset/articles/Examples_R.html
#source("./R/routine_tasks.R")
library(ComplexHeatmap)
genes <- c("SMC5", "SMC6", "NSMCE1", "NSMCE2", "NSMCE3", "NSMCE4A", "EID3")
cbpd <- readRDS("./data/cbioportal_curated.rds")
#cbpd <- cbpd[cbpd$major != 'Other',]
cbpd <- cbpd[cbpd$isalt == 1, ]
cbpd$major <- factor(cbpd$major)

makeUpset <- function(sdf, title, cut = 2) {
    m <- make_comb_mat(sdf, mode = 'distinct')
    # we remove any combination under 2 samples
    m <- m[comb_size(m) >= cut]
    # data at the top of teh upset plot
    breaks <- c(1, 10, 100, 1000, 2000)
    # top annotation bar plot
    top_ha <- HeatmapAnnotation(
        'Number of patients \n (log10 scale)' = anno_barplot(
                log10(comb_size(m))
                , bar_width = 0.8
                , gp = gpar(fill = '#3C3D41', color = '#3C3D41', fontface = 'bold')
                , border = FALSE
                , height = unit(5, "cm")
                , ylim = c(1,log10(2000))
                , axis_param = list(at = log10(breaks), labels = breaks))
        , annotation_name_side = "left"
        , annotation_name_rot = 90
    )
    # right annotation bar plot
    right_ha = rowAnnotation(
        'Set size \n (log10 scale)' = anno_barplot(
                log10(set_size(m))
                , bar_width = 0.8
                , gp = gpar(fill = '#3C3D41', color = '#3C3D41', fontface = 'bold')
                , border = FALSE
                , width = unit(5, "cm")
                , ylim = c(1, log10(2000))
                , axis_param = list(at = log10(breaks), labels = breaks)
                )
        , annotation_name_side = "bottom"
        )
    # main upset plot
    UpSet(m, comb_order = rev(order(comb_size(m)))
        , top_annotation = top_ha
        , row_names_side = 'left'
        , left_annotation = right_ha
        )
}

pdf('./figures/fig1c.pdf', onefile = TRUE, width = 8, height = 5)
plot(makeUpset(cbpd[, genes], 'All', 10))
dev.off()

print('done')