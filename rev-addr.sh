#!/bin/sh
# for no amp
Rscript ./reviewer-addressing/noAmp-Fig3ab.R
Rscript ./reviewer-addressing/noAmp-Fig3cd.R
Rscript ./reviewer-addressing/noAmp-FigS4ab.R
# for amp only
Rscript ./reviewer-addressing/yesAmp-Fig3ab.R
Rscript ./reviewer-addressing/yesAmp-Fig3cd.R
Rscript ./reviewer-addressing/yesAmp-FigS4ab.R
# lfc analysis between instability and smc5-6
Rscript ./reviewer-addressing/metabric-brca-lfc-comp.R
# co-mutation plots
Rscript ./reviewer-addressing/co-mutation-dot.R
# for transcriptational signature
Rscript ./reviewer-addressing/expr-Fig3cd.R
Rscript ./reviewer-addressing/expr-figS4cd.R
# new oncoprints
Rscript ./reviewer-addressing/oncoWinst.R
# new figure 2 w/o cancer type row
Rscript ./reviewer-addressing/newFig2.R
# heterogenity / homogenoity / bialeic analysis of TCGA data
Rscript ./reviewer-addressing/TCGA_hetero_ana.R
