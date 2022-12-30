# purpose: makes the dataframe for the alteration analysis of the complexes in the genes
source("./R/routine_tasks.R")
enable_parallelization()
# start
df <- readRDS("./data/cbioportal/raw.rds")

repeated_samples <- df$name[duplicated(df$name)]
print(head(repeated_samples))
print(length(repeated_samples))
print(dim(df))
print("checking for non-equivalent replicates")
for (i in repeated_samples) {
        temp_df <- df[df$name == i,]
        if(nrow(temp_df) == 1) next()
        for (j in 2:nrow(temp_df)) {
                if(sum(temp_df[j,] != temp_df[j-1,], na.rm = TRUE) != 0) stop()
        }
        df <- df[df$name != i, ]
        df <- rbind(df, temp_df[1, ])
}
print("deduplication done")
print(dim(df))
rownames(df) <- df$name
saveRDS(df, "./data/cbioportal/deduplicated.rds")
print("done")