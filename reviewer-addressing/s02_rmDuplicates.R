# a purpose: makes the dataframe for the alteration analysis of the complexes in the genes
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
    temp_df <- df[df$name == i, ]
    if (nrow(temp_df) == 1) next()
    torm <- c()
    for (j in 2:nrow(temp_df)) {
        # if they are not the same, we look at the one with more info
        # if (sum(temp_df[j, ] != temp_df[j - 1, ], na.rm = TRUE) != 0) {
        j_in_jmin1 <- sum(temp_df[j, ] %in% temp_df[j - 1, ])
        jmin1_in_j <- sum(temp_df[j - 1, ] %in% temp_df[j, ])
        # if j has more info in j - 1
        if (j_in_jmin1 > jmin1_in_j) {
            torm <- c(torm, j)
        } else {
            torm <- c(torm, j - 1)
        }
        # }
    }
    temp_df <- temp_df[-torm, ]
    if (nrow(temp_df) > 1) stop()
    df <- df[df$name != i, ]
    df <- rbind(df, temp_df[1, ])
}
print("deduplication done")
print(dim(df))
rownames(df) <- df$name
saveRDS(df, "./data/cbioportal/deduplicated.rds")
print("done")
