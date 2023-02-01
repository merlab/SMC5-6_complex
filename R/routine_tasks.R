genes <- c("SMC5", "SMC6", "NSMCE1", "NSMCE2", "NSMCE3", "NSMCE4A", "EID3")
assign("genes", genes, envir = .GlobalEnv)
tissue_types <- c("Breast", "Prostate", "Melanoma", "Ovarian", "Endometrial", "Lung", 'Pancreatic', 'Bladder'
  , 'Hepatobiliary', 'Esophagogastric', 'Bladder')
top_genes <- c('RAD54B', 'CCNE2', 'RECQL4', 'MCM4', 'AURKA', #'NSMCE2',
  'KIF13B', 'PPP2R2A', 'PARP3', 'TP53BP1'
  )
assign("top_genes", top_genes, envir = .GlobalEnv)
#colors <- c('#ef8a62','#67a9cf')
#color <- c("cornflowerblue", "firebrick")
#colors <- c('orange', "purple")
# gene names and column names set by default
# like the cut command in GNU cor utils
# takes the input, sepeartes by delemiter, gives the field you specify
#b_rank <- na.omit(setNames(b$logFC, rownames(b)))

make_pathway_db <- function() {
      pathway_db <- gmtPathways("./data/c2.cp.v7.5.1.symbols.gmt")
      pathway_db <- pathway_db[!duplicated(pathway_db)]
      return(pathway_db)
}
pathway_sentiment_ana <- function(input) {
  x <- input
  for(i in 1:length(x)) {
    t <- x[i]
    t <- strsplit(t, '_')[[1]]
    t <- paste(t[-1], collapse = ' ')
    x[i] <- t
  }
  x <- gsub('_', ' ', x)
  x <- gsub("THE ROLE OF", "", x)
  x <- gsub("PATHWAY", "", x)
  x <- gsub(" THE ", " ", x)
  x <- gsub(" AND ", " ", x)
  x <- gsub(" WITH ", " ", x)
  x <- gsub("MEDIATED ", " ", x)
  x <- gsub("RESPONSE ", " ", x)
  x <- tolower(x)
  all <- unique(unlist(strsplit(x, ' ')))
  all <- sapply(all, function(i) ifelse(nchar(i) > 1, i, NA))
  all <- unname(all)
  all <- na.omit(all)
  tdm <- matrix(0, nrow = length(all), ncol = length(x))
  rownames(tdm) <- all
  colnames(tdm) <- x
  for(i in all) {
    for(j in x) {
      if(grepl(i,j, fixed = TRUE)) tdm[i,j] <- tdm[i,j] + 1
    }
  }
  colnames(tdm) <- input
  v <- rowSums(tdm)
  v <- v[order(v, decreasing = TRUE), drop = FALSE]
  categories <- list(
    Signaling = c('receptor', 'signaling', 'activation', 'homeostasis', 'binding', 'uptake', 'interactions', 'mapk', 'activates', 'ion', 'secretion'),
    ECM = c('tight', 'junction', 'plasma', 'ecm', 'integrin', 'basement', 'proteoglycans', 'glycoproteins', 'matrix', 'extracellular', 'membranes', 'collagens', 'fibrils'),
    Immunity = c('infection', 'parasite', 'viral', 'inflammatory', 'immunoregulatory', 'antigen', 'lymphoid', 'il10', 'nf', 'kb', 'lymphoid', 'interferon', 'interleukin'),
    MuscleCardiac = c('muscle', 'contraction', 'striated', 'cardiac', 'myocarditis', 'cardiomyopathy'),
    Metabolism = c('glycosylation', 'tca', 'glycolysis', 'metabolism', 'gluconeogenesis', 'phosphorylation', 'mitochondrial', 'break', 'assembly', 'cytochrome', 'viral', 'degradation', 'p450', 'insulin'),
    CellCycle = c('checkpoints', 'cellular', 'replication', 'nuclear', 'checkpoint', 'g2', 'mitotic'),
    DNA = c('dna', 'mrna', 'genes', 'nucleus', 'transcription', 'transcriptional', 'replication'),
    Disease = c('cancer', 'telomeres','telomere', 'disease') # Hormonal = c('insulin'),
  )
  # cancer, p53, telomere, formation, disease
  # print(v)
  type <- data.frame(matrix(0, nrow = length(x), ncol = 2 + length(categories)))
  rownames(type) <- input
  colnames(type) <- c('pathway', 'category', names(categories))
  type$pathway <- x
  for(i in objects(categories)) {
    for(j in categories[[i]]) {
      if(!j %in% rownames(tdm)) next()
      print(j)
      type[,i] <- type[,i] + ifelse(tdm[j,] > 0, 1, 0)
    }
  }
  for(i in 1:nrow(type)) {
    lbl <- (names(categories)[(match(max(type[i, 3:ncol(type)]), type[i,])) - 2])
    if(length(lbl) == 0) lbl = 'Other'
    type$category[i] <- lbl
  }
  return(type)

}


cut_name <- function(input, delemiter = ".rds", field = 1) {
    return(strsplit(input, delemiter)[[1]][field])
}

# Enabled parallelization by mclapply and foreach
enable_parallelization <- function(cores = 0, ignore_error = FALSE) {
    library(parallel)
    library(foreach)
    library(doParallel)
    if (cores == 0) {
        registerDoParallel(cores = detectCores())
    } else if (cores <= detectCores() | ignore_error == TRUE) {
            registerDoParallel(cores = cores)
    } else if (ignore_error == FALSE) {
            return(warning("More cores chosen than what computer has!
                            Cores set to max the computer has"))
            registerDoParallel(cores = detectCores())
    } else {
        return(warning("Parallelization not enabled!"))
    }
}

# computes pearson, spearman, kendall correlation
# input must be a list containing two items:
# list(metadata, df with 2 rows)
# df is used for correlation computation
# returns the output in the following format:
# vactor(metadata, peasonCor, pearsonPval, spearmanCor, spearmanPval,
#       kendallCor, kendallPvak)
all_cors <- function(x) {
    df <- x[[2]]
    df <- df[!is.infinite(df$V1), ]
    df <- df[!is.infinite(df$V2), ]

    pl <- (cor.test(df$V1, df$V2, method = "pearson", use = "complete.obs"))
    s <- (cor.test(df$V1, df$V2, method = "spearman", use = "complete.obs"))
    k <- (cor.test(df$V1, df$V2, method = "kendall", use = "complete.obs"))

    output <- c(x[[1]], p$estimate, p$p.value,
                                 s$estimate, s$p.value,
                                 k$estimate, k$p.value)
    return(output)
}

# checks if folder exists - makes the folder if not present
folder_check <- function(inputdir) {
    if (file.exists(inputdir)) {
    } else {
        dir.create(inputdir)
    }
}


remove_special_chr <- function(x) {
    for(i in 1:length(x)) {
        x[i] <- gsub("[[:punct:]]", "", x[i])
    }
    return(x)
}


censor_data <- function(cl) {
    censor <- 60
    cl$DFS[cl$DFT > censor] <- 0
    cl$PFS[cl$PFT > censor] <- 0
    cl$OVS[cl$OVT > censor] <- 0
    cl$PFT[cl$PFT > censor] <- censor
    cl$OVT[cl$OVT > censor] <- censor
    cl$DFT[cl$DFT > censor] <- censor
    return(cl)
}

center_title <- function(l) {
  for(i in 1:length(l)) {
    l[[i]]$plot <- l[[i]]$plot + theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold')) 
  }
  return(l)
}
KM_plot_OVS <- function(fitisalt, fitnsmce2, title = NA, p = 0.05) {
  palette <- rev(c("#CB3814", "#3361BD"))
  pisalt <- surv_pvalue(fitisalt)$pval
  pnsmce2 <- surv_pvalue(fitnsmce2)$pval
  if(pisalt <= p | pnsmce2 <= p ) {
  #if(pisalt <= p & pnsmce2 <= p ) {
    pl <- list()
    print(pnsmce2)
    print(pisalt)
    pl[[1]] <- ggsurvplot(fitisalt, conf.int = TRUE, pval = TRUE
                    , risk.table = TRUE, palette = palette
                    , xlab = "Overall Survival (Months)"
                    , ylab = "Probability of Overall Survival"
                    , legend = c(0.8,0.95)
                    , risk.table.height = 0.20
                    , title = 'Complex'
                    , legend.labs = c('Wild', 'Altered')
                  )
    pl[[2]] <- ggsurvplot(fitnsmce2, conf.int = TRUE, pval = TRUE
                    , risk.table = TRUE, palette = palette
                    , xlab = "Overall Survival (Months)"
                    , ylab = "Probability of Overall Survival"
                    , legend = c(0.8,0.95)
                    , risk.table.height = 0.20
                    , title = 'NSMCE2' 
                    , legend.labs = c('Wild', 'Altered')
                  )
    arrange_ggsurvplots(center_title(pl), print = TRUE, title = paste(title, 'OVS'))
  }
}
KM_plot_PFS <- function(fitisalt, fitnsmce2, title = NA, p = 0.05) {
  palette <- rev(c("#CB3814", "#3361BD"))
  pisalt <- surv_pvalue(fitisalt)$pval
  pnsmce2 <- surv_pvalue(fitnsmce2)$pval
  #print(pisalt)
  #print(pnsmce2)
  if(pisalt <= p | pnsmce2 <= p ) {
  #if(pisalt <= p & pnsmce2 <= p ) {
    pl <- list()
    print(pnsmce2)
    print(pisalt)
    pl[[1]] <- ggsurvplot(fitisalt, conf.int = TRUE, pval = TRUE
                    , risk.table = TRUE, palette = palette
                    , xlab = "Progression free survival (Months)"
                    , ylab = "Probability of Progression free survival"
                    , legend = c(0.8,0.95)
                    , risk.table.height = 0.20
                    , title = 'Complex'
                    , legend.labs = c('Wild', 'Altered')
                  )
    pl[[2]] <- ggsurvplot(fitnsmce2, conf.int = TRUE, pval = TRUE
                    , risk.table = TRUE, palette = palette
                    , xlab = "Progression free survival (Months)"
                    , ylab = "Probability of Progression free survival"
                    , legend = c(0.8,0.95)
                    , risk.table.height = 0.20
                    , title = 'NSMCE2' 
                    , legend.labs = c('Wild', 'Altered')
                  )
    arrange_ggsurvplots(center_title(pl), print = TRUE, title = paste(title, 'PFS'))
  }
}
KM_plot_DFS <- function(fitisalt, fitnsmce2, title = NA, p = 0.05) {
  palette <- rev(c("#CB3814", "#3361BD"))
  pisalt <- surv_pvalue(fitisalt)$pval
  pnsmce2 <- surv_pvalue(fitnsmce2)$pval
  if(pisalt <= p | pnsmce2 <= p ) {
  #if(pisalt <= p & pnsmce2 <= p ) {
    pl <- list()
    print(pnsmce2)
    print(pisalt)
    pl[[1]] <- ggsurvplot(fitisalt, conf.int = TRUE, pval = TRUE
                    , risk.table = TRUE, palette = palette
                    , xlab = "Disease free status (Months)"
                    , ylab = "Probability of Disease free status"
                    , legend = c(0.8,0.95)
                    , risk.table.height = 0.20
                    , title = 'Complex'
                    , legend.labs = c('Wild', 'Altered')
                  )

    pl[[2]] <- ggsurvplot(fitnsmce2, conf.int = TRUE, pval = TRUE
                    , risk.table = TRUE, palette = palette
                    , xlab = "Disease free status (Months)"
                    , ylab = "Probability of Disease free status"
                    , legend = c(0.8,0.95)
                    , risk.table.height = 0.20
                    , title = 'NSMCE2' 
                    , legend.labs = c('Wild', 'Altered')
                  )
    arrange_ggsurvplots(center_title(pl), print = TRUE, title = paste(title, 'DFS'))
  }
}

make_label <- function(v) {
  for(i in 1:length(v)) {
    x <- v[i]
    x <- strsplit(x, '_')[[1]]
    x <- paste(x[-1], collapse = ' ')
    v[i] <- x
  }
  # abbreviation
  v <- gsub("NUCLEAR PORE COMPLEX", "NPC", v)
  v <- gsub("NPC NPC", "NPC", v)
  v <- gsub("DOUBLE STRANDED BREAK", "DSB", v)
  v <- gsub("DOUBLE STRAND BREAK", "DSB", v)
  # shortening of documents
  v <- gsub("MECHANISM", "MECH.", v)
  v <- gsub("REGULATION OF", "REG.", v)
  v <- gsub("CELLULAR", "CELL", v)
  v <- gsub("NEGATIVE", "NEG.", v)
  v <- gsub("SYSTEM", "SYS.", v)
  v <- gsub("UBIQUITINATION", "UB.", v)
  v <- gsub("REMOVAL", "RM.", v)
  v <- gsub("REMOVE", "RM.", v)
  v <- gsub("EXPRESSION", "EXP.", v)
  v <- gsub("ACTIVITY", "ACTY.", v)
  v <- gsub("ACTIVATION", "ACTY.", v)
  v <- gsub("FORMATION", "FORMN.", v)
  v <- gsub("REPLICATION", "RPL.", v)
  v <- gsub("REPLICATIVE", "RPL.", v)
  v <- gsub("DEGRADATION", "DEG.", v)
  v <- gsub("INITIATION", "INIT.", v)
  v <- gsub("MOLECULES", "MOLEC.", v)
  v <- gsub("MOLECULAR", "MOLEC.", v)
  v <- gsub("MOLECULE", "MOLEC.", v)
  v <- gsub("ASSOCIATED", "ASSOC.", v)
  v <- gsub("ASSOCIATION", "ASSOC.", v)
  v <- gsub("RESPONSE", "RESP.", v)
  v <- gsub("SIGNALING", "SIG.", v)
  v <- gsub("DAMAGE", "DMG.", v)
  v <- gsub("RECOGNITION", "RECOG.", v)
  v <- gsub("REPROGRAMMING", "REPROG.", v)
  v <- gsub("PROGRAMMING", "PROG.", v)
  v <- gsub("INTERACTIONS", "INT.", v)
  v <- gsub("INTERACTION", "INT.", v)
  v <- gsub("CYTOCHROME P450", "CYP450", v)
  v <- gsub("NETWORK", "NTW.", v)
  v <- gsub("METABOLIC", "METAB.", v)
  v <- gsub("METABOLISM", "METAB.", v)
  v <- gsub("CHECKPOINTS", "CKPT.", v)
  v <- gsub("CHECKPOINT", "CKPT.", v)
  v <- gsub("MAINTENANCE", "MAINT.", v)
  v <- gsub("SYNTHESIS", "SYN.", v)
  v <- gsub("COMPLEX", "CMPLX.", v)
  v <- gsub("EXTENSION", "EXT.", v)
  v <- gsub("PROGRESSION", "PROG.", v)
  v <- gsub("PHOSPHORYLATION", "PY.", v)
  v <- gsub("STABILITY", "STB.", v)
  v <- gsub("TARGETTED", "TARGET", v)
  v <- gsub("TELOMERES", "telomere", v)
  # lexical simplication
  v <- gsub("THE ROLE OF", "ROLE OF", v)
  v <- gsub(" PATHWAY ", " ", v)
  v <- gsub(" THE ", " ", v)
  v <- gsub(" AND ", " & ", v)
  v <- gsub(" WITH ", " w. ", v)
  v <- gsub(" THROUGH ", " via ", v)
  v <- gsub(" VIA ", " via ", v)
  v <- gsub(" MEDIATED ", " med. ", v)
  v <- gsub(" OF ", " of ", v)
  v <- gsub(" TO ", " to ", v)
  v <- gsub(" IN ", " in ", v)
  v <- gsub(" PRE ", " pre ", v)
  v <- gsub(" DERIVED FROM ", " from ", v)
  return(v)
    # x <- gsub("", "METAB.", x)
    # x <- gsub("DISEASE", "METAB.", x)
    # x <- gsub("BETA", "β", x)
    # x <- gsub("ALPHA", "α", x)
}
find_n_sim_genes <- function(v1,v2) {
  v1 <- v1[[1]]
  v2 <- v2[[1]]
  s1 <- v1[v1 %in% v2]
  s2 <- v2[v2 %in% v1]
  if(length(s1) != length(s2)) stop('something is wrong')
  return(length(s1))
}

# https://stackoverflow.com/questions/7455046/how-to-make-graphics-with-transparent-background-in-r-using-ggplot2
rmbg <- function(g) {
    return(g + theme(
          panel.background = element_rect(fill = "transparent", color = NA), # bg of the panel
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          legend.background = element_rect(fill = "transparent", color = NA), # get rid of legend bg
          legend.box.background = element_rect(fill = "transparent", color = NA), # get rid of legend panel bg
          panel.grid.major = element_blank(), # get rid of major grid
          panel.grid.minor = element_blank(), # get rid of minor grid
          panel.border = element_blank()
        ))
}
#rmallaxis <- function()
#rmalllegend <- function()
#addtxt <- function()
#absylim <- function(o, y1,y2) return(o + scale_y_continuous(limits = c(y1,y2))




rainCloudPlot <- function(mRNA, mutation, type) {
    row <- 'isalt'
    colors <- c("#CB3814", "#3361BD")
    df <- as.data.frame(t(mRNA))
    genes <- colnames(df)
    df <- cbind(df, mutation[,row])
    df <- as(tidyr::pivot_longer(df, genes), 'data.frame')
    df$value <- as.numeric(df$value)
    colnames(df)[1] <- 'altStat'
    df$altStat <- as.factor(ifelse(df$altStat == '1', 'Altered', 'Wild'))
    out <- list()
    for(i in unique(df$name)) {
      sdf <- df[df$name == i, ]
      p <- (wilcox.test(sdf$value[sdf$altStat == 'Wild'], sdf$value[sdf$altStat == 'Altered']))$p.value
      p <- signif(p, 1)
      # format p value
      if(p == 3e-12) p <- bquote('P' ~ '=' ~ '3' ~ 'x' ~ 10^-12)
      if(p >0.9e-70 & p < 1.1e-70) p <- bquote('P' ~ '=' ~ '1' ~ 'x' ~ 10^-70)
      if(p == 3e-61) p <- bquote('P' ~ '=' ~ '3' ~ 'x' ~ 10^-61)
      if(p > 1.9e-57 & p < 2.1e-57) p <- bquote('P' ~ '=' ~ '2' ~ 'x' ~ 10^-57)
      if(p == 4e-49) p <- bquote('P' ~ '=' ~ '4' ~ 'x' ~ 10^-49)
      if(p == 6e-30) p <- bquote('P' ~ '=' ~ '6' ~ 'x' ~ 10^-30)
      if(p == 3e-26) p <- bquote('P' ~ '=' ~ '3' ~ 'x' ~ 10^-26)
      if(p > 1.9e-29 & p < 2.1e-29) p <- bquote('P' ~ '=' ~ '2' ~ 'x' ~ 10^-29)
      if(p == 2e-20) p <- bquote('P' ~ '=' ~ '2' ~ 'x' ~ 10^-20)

      #print(head(sdf))
      plt <- ggplot(sdf, aes(x=altStat , y= value, fill = altStat, color = altStat)) +
        geom_boxplot(width = .6, color = 'black') +
        annotate(geom = 'text', x = 1.5, y = max(sdf$value) + 0.5, label = p) +
        xlab('') +
        guides(fill = FALSE, color = FALSE) +
        theme_cowplot()+ 
        ggtitle(i) +
        scale_colour_manual(values = colors) + scale_fill_manual(values = colors) +
        theme(plot.title=element_text(hjust=0.5, size = 11, face = 'bold'))
      
      print(type == 'microarray')

      if(type == 'microarray') plt <- plt + ylab('Gene expression')
      if(type == 'rnaseq') plt <- plt + ylab('Count')
      out[[i]] <- rmbg(plt)
    }
    return(out)
}


# give a mRNA expression matrix, the function finds the mutaiton dataset matching it
obtain_mut_from_mRNA <- function(mRNA) {
      samples <- colnames(mRNA) 
      cbpd <- readRDS("./data/cbioportal/format_exOther.rds")
      cbpd <- cbpd[cbpd$name %in% samples, ]
      rownames(cbpd) <- cbpd$name
      cbpd$ploidy <- cbpd$aneuploidyScore <- cbpd$study <- cbpd$tissue <- cbpd$tissue_detailed <- cbpd$name <- cbpd$noSample <- cbpd$sex <- NULL
      return(cbpd)
}


# perform genedrug association given a genelist
calc_gsea <- function(gene_list, rdsout, xlsxout, minSize = 35, maxSize = 145) {
      gene_list <- gene_list[order(gene_list, decreasing = TRUE)]
      #print(head(gene_list))
      vector <- gene_list
      genes <- c("SMC5", "SMC6", "NSMCE1", "NSMCE2", "NSMCE3", "NSMCE4A", "EID3")
      o <- fgsea(pathways = make_pathway_db(),
                      stats = gene_list,
                      minSize = minSize,
                      maxSize = maxSize)
      o <- o[order(o$padj), ]
      saveRDS(o, rdsout)
      o$leadingEdge <- NULL
      write_xlsx(o, xlsxout)
}

# clean the pathway names
clean_path_name <- function(vector) {
      vector <- gsub('GOBP_', '', vector)
      vector <- gsub('KEGG_', '', vector)
      vector <- gsub('HALLMARK_', '', vector)
      vector <- gsub('_', ' ', vector)
      return(vector)
}
