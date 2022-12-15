#!/bin/sh
Rscript ./R/1.oncoAna.R
Rscript ./R/2.rmDuplicates.R
Rscript ./R/3.format.R
Rscript ./R/fig1b.R
Rscript ./R/fig1c.R
Rscript ./R/fig2.R
Rscript ./R/fig3ab.R
Rscript ./R/fig3cd.R
