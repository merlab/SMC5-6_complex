# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8670522/
library(fgsea)
library(ggpubr)
library(ggplot2)
library(ggrepel)
set.seed(123)
source('./R/routine_tasks.R')
colors <- c('#d95f02', '#7570b3', '#e7298a')

library(fgsea)
library(ggplot2)
library(ggrepel)
library(ggpubr)
source('./R/routine_tasks.R')
set.seed(123)
n <- 1
rm_if_p_common <- 0.90

make_label_local <- function(v) {
  for(i in 1:length(v)) {
    x <- strsplit(v[i], '_')[[1]]
    x <- x[-1]
    x <- paste(x, collapse = ' ')
    print(x)
    v[i] <- x
  }
  # abbreviation
  v <- gsub('_', ' ', v)
  v <- gsub("NUCLEAR PORE COMPLEX", "NPC", v)
  v <- gsub("NPC NPC", "NPC", v)
  v <- gsub("DOUBLE STRANDED BREAK", "DSB", v)
  v <- gsub("DOUBLE STRAND BREAK", "DSB", v)
  # shortening of documents
  v <- gsub("MOLECULES", "MOLEC.", v)
  v <- gsub("MOLECULAR", "MOLEC.", v)
  v <- gsub("MOLECULE", "MOLEC.", v)
  v <- gsub("CYTOCHROME P450", "CYP450", v)
  v <- gsub("EXTENSION", "EXT.", v)
  v <- gsub("PROGRESSION", "PROG.", v)
  v <- gsub("TELOMERES", "telomere", v)
  # lexical simplication
  v <- gsub(" PATHWAY ", " ", v)
  v <- gsub(" AND ", " & ", v)
  return(v)
}

plotEnrichment <- function (pathway, stats, gseaParam = 1, ticksSize = 0.2, title = NA
    , p = 0, fdr = 0, es = 0, nes = 0, color = 'green') {
    rnk <- rank(-stats)
    ord <- order(rnk)
    statsAdj <- stats[ord]
    statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
    statsAdj <- statsAdj/max(abs(statsAdj))
    pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
    pathway <- sort(pathway)
    gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, 
        returnAllExtremes = TRUE)
    bottoms <- gseaRes$bottoms
    print(min(bottoms))
    tops <- gseaRes$tops
    n <- length(statsAdj)
    xs <- as.vector(rbind(pathway - 1, pathway))
    ys <- as.vector(rbind(bottoms, tops))
    toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
    x = y = NULL
    yannot <- min(bottoms)
    xannot <- 8000#18000 
    size <- 2
    g1 <- ggplot(toPlot, aes(x = x, y = y)) + 
        geom_point(color = color, size = 0.1) + 
        geom_line(color = color) + 
        ylab('Enrichment score') +
        xlab('') +
        theme_classic() + 
        scale_x_continuous(expand = c(0,0), limits = c(0,length(stats)), position="bottom") +
        theme(axis.line.x = element_blank(), axis.text.x = element_blank()
            , axis.ticks.x = element_blank(),
              axis.title.y = element_text(size = 8))
    g1 <- rmbg(g1)

    g2 <- ggplot(toPlot, aes(x = x, y = y)) +  
        geom_segment(data = data.frame(x = pathway), mapping = aes(x = x, 
            y = -0.25, xend = x, yend = 0.25), size = ticksSize) + 
        ylab('') +
        xlab('Gene rank') +
        scale_x_continuous(expand = c(0,0), limits = c(0,length(stats)), position="bottom",
          breaks = seq(0, 25000, 5000)) +
        scale_y_continuous(limits = c(-0.25,0.25)) +
        ylim(-0.25,0.25) +
        geom_hline(yintercept = 0.25, color = 'black') +#, size = 1) +
        geom_vline(xintercept = n, color = 'black') +#, size = 1) +
        theme_classic() +
        theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
              axis.ticks.y = element_blank(), axis.text.y = element_blank(),
              ) 
    g2 <- rmbg(g2)
    ghigh <- ggplot() +
        annotate(geom = 'text', hjust = 0.5, label = title, color = 'black'
                , x = 1, y = 1, size = 4, fontface = 2) +
        xlab('') + ylab('') +
        theme_void()
    ghigh <- rmbg(ghigh)
    glow <- ggarrange(g1, NULL, g2, nrow = 3, ncol = 1, heights= c(0.6, -0.08, 0.3), align = 'hv')
    glow <- rmbg(glow)
    g <- ggarrange(ghigh, NULL, glow, nrow = 3, ncol = 1, heights = c(0.2, -0.05, 0.9), align = 'hv')
    g <- rmbg(g)
    g <- g + border()
    return(g)
}

df <- readRDS('./results/metabric-brca/limma.rds')
rank <- setNames(df$logFC, rownames(df))
rank <- rank[order(rank, decreasing = TRUE)]


pathways <- make_pathway_db()
o <- fgsea(pathways = pathways,
                stats = rank,
                minSize = 25,
                maxSize = 150)

p <- list()
chosen_pathways <- c(
  'REACTOME_PEPTIDE_HORMONE_METABOLISM', 
  'REACTOME_ECM_PROTEOGLYCANS',
  "REACTOME_DNA_REPLICATION",
  'REACTOME_ACTIVATION_OF_ATR_IN_RESPONSE_TO_REPLICATION_STRESS'
  )
for(i in chosen_pathways) {
    p[[i]] <- plotEnrichment(pathways[[i]],rank, title = make_label_local(i)
      , p = o$pval[o$pathway == i], fdr = o$padj[o$pathway == i], es = o$ES[o$pathway == i], nes = o$NES[o$pathway == i]
      , color = colors[3]) 
}
enrichment <- ggarrange(plotlist = p, nrow = 4, ncol = 1, align = 'hv')
pdf('./figures/fig4row2.pdf', width = 5.5, height = 12)
x <- rmbg(enrichment)
plot(x)
dev.off()
print('done')