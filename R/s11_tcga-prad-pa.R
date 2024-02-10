library(fgsea)
make_pathway_db <- function() {
  pathway_db_dir <- "./data"
  pathway_db <- NA
  pathway_db_names <- c("c2.cp.v7.5.1.symbols.gmt")
  for (i in pathway_db_names) {
    pathway_db <- c(
      pathway_db,
      gmtPathways(sprintf("%s/%s", pathway_db_dir, i))
    )
    print(paste("added", i))
  }
  pathway_db <- pathway_db[!duplicated(pathway_db)]
  # print(paste("number of pathways:", length(pathway_db)))
  return(pathway_db)
}
# perform genedrug association given a genelist


df <- readRDS("./data/tcga-prad/dgea-edgeR.rds")
gene_list <- na.omit(setNames(df$logFC, rownames(df)))

gene_list <- gene_list[order(gene_list, decreasing = TRUE)]
vector <- gene_list
genes <- c("SMC5", "SMC6", "NSMCE1", "NSMCE2", "NSMCE3", "NSMCE4A", "EID3")
o <- fgsea(
  pathways = make_pathway_db(),
  stats = gene_list,
  minSize = 35,
  maxSize = 145
)
o <- o[order(o$padj), ]
saveRDS(o, "./data/tcga-prad/pathway-gsea.rds")
print("done")

