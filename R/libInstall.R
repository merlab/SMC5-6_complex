packages <- c('devtools', 'Rcpp','sf', 
  'cowplot', 'readr',  'dplyr', 'lavaan', 'Hmisc', 'tidyr', 'igraph',
  'ggplot2','ggpubr', 'ggraph', 'ggrepel', 'circlize', 'network', 
  'networkD3', 'gridExtra','BiocManager', 'readxl', 'writexl', 'remotes', 
  'foreach', 'doParallel', 'ndtv', 'parallel', 'data.table',
  'survminer', 'survival', 'scales', 'sna', 'threejs', 'visNetwork'
    )
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))
}

biocpackages <- c(
  'fgsea', 'edgeR', 'limma', 'DOSE',
  'MetaGxBreast', 'TCGAbiolinks', 'Biobase',
  'clusterProfiler', 'enrichplot', 'ComplexHeatmap'

)
if (length(setdiff(biocpackages, rownames(installed.packages()))) > 0) {
  BiocManager::install(setdiff(biocpackages, rownames(installed.packages())))
}

