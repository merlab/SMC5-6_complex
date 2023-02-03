source('./R/routine_tasks.R')
#library(MetaGxBreast)
library(Biobase)
enable_parallelization(cores = detectCores() %/% 4 + 1)

esets <- readRDS('./data/metabric-brca-metagx-download.rds')

eset <- esets[[1]]
expmat <- assayData(eset)[['exprs']]
colnames(expmat) <- gsub('_', '-', colnames(expmat))

deduplicate_genes <- function(features, expmat) {
    duplicated_genes <- unique(as.character(features[duplicated(features$gene),'gene']))
    if(length(duplicated_genes) == 0) return(expmat)
    for(i in duplicated_genes) {
        x <- (rownames(features)[features$gene == i])
        subexpmat <- expmat[x,]
        rvar <- (apply(subexpmat, 1, var))
        rowstorm <- names((which(rvar != max(rvar))))
        if(length(rowstorm) == 0) next()
        expmat <- expmat[!rownames(expmat) %in% rowstorm,]
    }
    return(expmat)
}

# remove genes with NA values
rvar <- apply(expmat, 1, var)
rvar <- rvar[!is.na(rvar)]
expmat <- expmat[names(rvar),]
features <- as(featureData(eset), 'data.frame')
features <- features[rownames(expmat),]

expmat <- expmat[order(features$gene),]
features <- features[rownames(expmat),]

expmattemp <- expmat
chunks <- c(seq(25,1,-6))
for(x in 1:length(chunks)) {
    n <- nrow(expmattemp)
    print(n)
    l <- n %/% chunks[x]
    if(chunks[x] == 1) {
      expmattemp <- deduplicate_genes(features[1:n,], expmattemp[1:n,])
      features <- features[rownames(expmattemp), ]
      break()
    } else {
      expmattemp <- foreach(i = 1:(chunks[x]+1), .combine = 'rbind') %dopar% {
        # subsetting work relatively equally for all cores
          from <- (l*(i-1))+1
          to <- l*i
          if(i > chunks[x]) to <- n
          print(paste(from, '-', to))
          return(deduplicate_genes(features[from:to,], expmattemp[from:to,]))
      }
    }
    features <- features[rownames(expmattemp), ]
    if(length(unique(as.character(features[duplicated(features$gene),'gene']))) == 0) break()
}
rownames(expmattemp) <- features$gene
folder_check('./data/metabric-brca')
saveRDS(expmattemp, './data/metabric-brca/microarray_metagx.rds')
print('done')