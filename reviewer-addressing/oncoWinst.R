# purpose: makes a traditional oncoprint
genes <- c("SMC5", "SMC6", "NSMCE1", "NSMCE2", "NSMCE3", "NSMCE4A", "EID3")
dir <- getwd()
source("./R/routine_tasks.R")
library(ggplot2)
library(ggpubr)
suppressPackageStartupMessages(library(ComplexHeatmap))

###################
## Configuration ##
###################
# what to subset the oncoprint file to
db <- readRDS("./reviewer-addressing/cbioportal/format_exOther.rds")
rownames(db) <- db$name
db <- db[, c(grep("_det", colnames(db), value = TRUE), "major")]
db <- db[, -1]
colnames(db) <- gsub("_det", "", colnames(db))
db <- apply(db[, 1:ncol(db)], 2, function(x) {
    x <- gsub("\\(putative passenger\\)", "", x)
    x <- gsub("sv", "Structural variation", x)
    x <- gsub("splice", "Splice", x)
    x <- gsub("_rec", "", x)
    x <- gsub("\\(putative driver\\)", "", x)
    x <- gsub("amp", "Amplification", x)
    x <- gsub("homdel", "Homologous Deletion", x)
    x <- gsub(" $", "", x)
})

tissue_types <- c(
    "Breast", "Prostate", "Melanoma", "Ovarian", "Endometrial", "Lung", "Pancreatic", "Bladder",
    "Hepatobiliary", "Esophagogastric", "Bladder"
)
tissueName <- c(
    All = "All",
    Ovarian = "Ovarian cancer",
    Breast = "Breast cancer",
    Prostate = "Prostate cancer"
    # Endometrial = "Endometrial cancer",
    # Melanoma = "Melanoma",
    # Esophagogastric = "Esophagogastric cancer",
    # Bladder = "Bladder cancer",
    # Hepatobiliary = "Hepatobiliary cancer",
    # Pancreatic = "Pancreatic cancer",
    # Lung = "Lung cancer"
)

pdf("./reviewer-addressing/oncoWinst.pdf", height = 6, width = 11, onefile = TRUE)
for (i in names(tissueName)) {
    tissue <- i
    if (tissue == "All") {
        mat <- db[db[, "major"] != "Other", ]
    } else {
        mat <- db[db[, "major"] == i, ]
    }
    mat <- mat[, -ncol(mat)]
    print(table(unlist(mat)))
    mat <- t(mat)
    names <- c(
        "Amplification", "Structural variation", "Deep Deletion",
        "Inframe Mutation", "Missense Mutation", "Truncating mutation",
        "Splice"
    )
    colors <- c("blue", "red", "#008000", "#984ea3", "#ff7f00", "#f781bf", "#a6761d", "black")
    print(i)
    print(dim(mat))
    col <- setNames(colors, names)
    alter_fun <- list(
        background = alter_graphic("rect", fill = "#CCCCCC"),
        "Amplification" = alter_graphic("rect", fill = col["Amplification"]),
        "Deep Deletion" = alter_graphic("rect", fill = col["Deep Deletion"]),
        "Inframe Mutation" = alter_graphic("rect", fill = col["Inframe Mutation"]),
        "Truncating mutation" = alter_graphic("rect", fill = col["Truncating mutation"]),
        "Splice" = alter_graphic("rect", height = 0.33, fill = col["Splice"]),
        "Structural variation" = alter_graphic("rect", height = 0.33, fill = col["Structural variation"]),
        "Missense Mutation" = alter_graphic("rect", height = 0.33, fill = col["Missense Mutation"]),
        "Homologous Deletion" = alter_graphic("rect", height = 0.33, fill = col["Homologous Deletion"])
    )
    column_title <- tissueName[i]

    if (tissue == "All") {
        heatmap_legend_param <- list(
            title = "Genomic variations",
            at = names,
            labels = names
        )
        plot(oncoPrint(mat,
            alter_fun = alter_fun, col = col,
            column_title = column_title,
            heatmap_legend_param = heatmap_legend_param,
            alter_fun_is_vectorized = FALSE,
            column_names_gp = gpar(fontsize = 10),
        ))
    } else {
        plot(oncoPrint(mat,
            alter_fun = alter_fun, col = col,
            column_title = column_title,
            alter_fun_is_vectorized = FALSE,
            column_names_gp = gpar(fontsize = 10),
        ))
    }
}
dev.off()
