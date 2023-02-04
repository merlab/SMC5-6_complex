# config
trans <- 0.3
min_path <- 35
max_path <- 145
show_n_common_gene <- 20
rm_if_p_common <- 1
# program
library(gridExtra)
library(tidyr)
library(igraph) 
library(network) 
library(sna)
library(ggraph)
library(visNetwork)
library(threejs)
library(networkD3)
source('./R/routine_tasks.R')
set.seed(123)
colors <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#cab2d6','#6a3d9a','#ffff99','#b15928')
colors <- adjustcolor(colors, alpha = trans)

# NOTE: making edge & node data.fames
# TODO: fix names. never use PDF
pathDf <- as.data.frame(readRDS('./data/metabric-brca/pathway-gsea.rds'))
rownames(pathDf) <- pathDf$pathway
pathDf <- pathDf[pathDf$padj < 0.01 & pathDf$pval < 0.01, ]
nodes <- pathDf
nodes$label <- nodes$pathway
nodes <- nodes[nodes$size >= min_path, ]
nodes <- nodes[nodes$size <= max_path, ]
# NOTE: sentimen analysis
s <- pathway_sentiment_ana(nodes$label)
nodes$type  <- as.factor(s$category)
# NOTE: pathway removed due to lack of data
nodes <- nodes[nodes$type != 'Other', ]
nodes <- nodes[nodes$type != 'MuscleCardiac', ]
nodes$label <- make_label(nodes$label)
nodes$color <- colors[as.factor(nodes$type)]
# removing too small or large pathways
pathDf <- pathDf[nodes$pathway, ]

pathways_to_rm <- c()
edges <- data.frame()
for(i in 1:(nrow(pathDf)-1)) {
  from <- rownames(pathDf)[i]
  for(j in i:nrow(pathDf)) {
    to <- rownames(pathDf)[j]
    if(from == to) next()
    n <- find_n_sim_genes(pathDf$leadingEdge[i], pathDf$leadingEdge[j])
    if(n / min(pathDf$size[i]) >= rm_if_p_common) pathways_to_rm <- c(pathways_to_rm, from)
    if(n / min(pathDf$size[j]) >= rm_if_p_common) pathways_to_rm <- c(pathways_to_rm, to)
    w <- 0
    if(n > show_n_common_gene) w <- (n / 1000)
    # this is to ensure that pathways in one file are very close to each other
    #   as it improves attraction between members of the same group (line widht not changed)
    if(nodes$type[i] == nodes$type[j]) w <- w + 0.3 

    edges <- rbind(edges, c(from = from, to = to, width = n, weight = w))
  }
}
pathways_to_rm <- unique(pathways_to_rm)
colnames(edges) <- c('from', 'to', 'width', 'weight')
edges$width <- as.numeric(edges$width)
edges$weight <- as.numeric(edges$weight)
edges <- edges[edges$weight > 0,]

# NOTE: rmeove pathways that share above %T of genes with other pathways
nodes <- nodes[!nodes$pathway %in% pathways_to_rm, ]
edges <- edges[!edges$from %in% pathways_to_rm, ]
edges <- edges[!edges$to %in% pathways_to_rm, ]

net <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
E(net)$arrow.size <- 0
E(net)$edge.color <- nodes$type
x <- as.numeric(E(net)$width)
E(net)$width <- (x / 50) + 0.5
V(net)$size <- as.numeric(V(net)$size) / 20 + 0.3
V(net)$label <- as.character(nodes$type)

set.seed(123)
# different layout styles
fr2 <- layout_with_fr(net, niter = 5000, dim = 2)
kk <- layout_with_kk(net)
cs <- layout_with_graphopt(net, charge = 0.1, mass = 0.2, spring.length = 2, spring.constant = 0.1)

pdf('./figures/fig4c.pdf', width = 6, height = 6)
par(mar = c(2,2,2,2))
clp <- cluster_label_prop(net)
clp$membership <- nodes$type

new_cols <- colors[membership(clp)]
plot(clp
  , net
  , col='grey70'
  , mark.border="grey70"
  , mark.col = unique(new_cols)
  , edge.color = 'grey50'
  , vertex.label.cex = 0.1
  , layout = layout_nicely
  , remove.multiple = TRUE
  , remove.loops = TRUE
  , vertex.label=NA
  , vertex.frame.color=NA
  )

# Add labels
types <- unique(as.factor(nodes$type))
types <- as.character(types[order(as.numeric(types))])
colors <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#cab2d6')
colors <- adjustcolor(colors, alpha = trans)
types <- c(
  DNA = 'DNA replication & repair'
, Disease = 'Cancer-specific'
, CellCycle = 'Cell Cyle-related'
, Metabolism = 'Metabolic'
, Signaling = 'Signaling'
, ECM = 'Extracellular Matrix'
, Immunity = 'Immune reaction'
)
legend(x = -1.2, y = 1.1
     , types
     , col = colors
     , pt.bg = colors
     , pt.cex = 0.8, cex = 0.8, bty = "n", ncol = 1, pch = 21
     , title = 'Pathway Category', title.cex = 1
     )
dev.off()
print('fig4 igraph done')
x <- V(net)
df <- data.frame(names = names(x), type = x$label)
V(net)$label <- c(1:nrow(df))

set.seed(123)
# different layout styles
fr2 <- layout_with_fr(net, niter = 5000, dim = 2)
kk <- layout_with_kk(net)
cs <- layout_with_graphopt(net, charge = 0.1, mass = 0.2, spring.length = 2, spring.constant = 0.1)
pdf('./figures/figS2a.pdf', width = 6, height = 6)
par(mar = c(2,2,2,2))
clp <- cluster_label_prop(net)
clp$membership <- nodes$type

new_cols <- colors[membership(clp)]
plot(clp
  , net
  , col='grey70'
  , mark.border="grey70"
  , mark.col = unique(new_cols)
  , edge.color = 'grey50'
  , vertex.label.cex = 0.5
  , vertex.label.color = 'black'
  , layout = layout_nicely
  , remove.multiple = TRUE
  , remove.loops = TRUE
  ,vertex.frame.color=NA
  )

# Add labels
types <- unique(as.factor(nodes$type))
types <- as.character(types[order(as.numeric(types))])
colors <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#cab2d6')
colors <- adjustcolor(colors, alpha = trans)

types <- c(
  DNA = 'DNA replication & repair'
, Disease = 'Cancer-specific'
, CellCycle = 'Cell Cyle-related'
, Metabolism = 'Metabolic'
, Signaling = 'Signaling'
, ECM = 'Extracellular Matrix'
, Immunity = 'Immune reaction'
)
legend(x = -1.2, y = 1.1
     , types
     , col = colors
     , pt.bg = colors
     , pt.cex = 0.8, cex = 0.8, bty = "n", ncol = 1, pch = 21
     , title = 'Pathway Category', title.cex = 1
     )
dev.off()



pdf('./figures/figS2b.pdf', width = 20, height = 20)

types <- c(
  DNA = 'DNA replication & repair'
, Disease = 'Cancer-specific'
, CellCycle = 'Cell Cyle-related'
, Metabolism = 'Metabolic'
, Signaling = 'Signaling'
, ECM = 'Extracellular Matrix'
, Immunity = 'Immune reaction'
)
# format type column name
df$'Pathway Name' <- df$names
df$'Pathway Type' <- types[df$type]
df$'Enrichment Score' <- pathDf[df$names, 'ES']
df$'-log10FDR' <- -log10(pathDf[df$names, 'padj'])
df$'Pathway Name' <- gsub('_', ' ', df$'Pathway Name')
df$names <- NULL; df$type <- NULL
tt2 <- ttheme_default(core=list(fg_params=list(hjust=0, x=0)),
                      rowhead=list(fg_params=list(hjust=0, x=0)))
grid.table(df, theme = tt2)
dev.off()

print('done')