#!/bin/sh
Rscript ./R/1.oncoAna.R
Rscript ./R/2.rmDuplicates.R
Rscript ./R/3.format.R & \
  Rscript ./R/4.metabric-brca-format.R
Rscript ./R/5.metabric-brca-dgea.R
Rscript ./R/6.metabric-brca-gsea.R
Rscript ./R/fig1b.R & \
  Rscript ./R/fig1c.R & \
  Rscript ./R/fig2.R & \
  Rscript ./R/fig3ab.R & \
  Rscript ./R/fig3cd.R & \
  Rscript ./R/fig4ab.R & \
  Rscript './R/fig4c&figS2.R' & \
  Rscript ./R/figS1.R & \
  Rscript ./R/figS3.R & \
  Rscript ./R/figS4ab.R & \
  Rscript ./R/figS4cd.R 

