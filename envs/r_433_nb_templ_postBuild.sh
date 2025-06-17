#!/bin/bash

source ~/.bashrc
conda activate r_433_nb

Rscript -e "devtools::install_github('tpq/balance')"
Rscript -e "devtools::install_github('malucalle/selbal')"
Rscript -e "devtools::install_bitbucket('knomics/nearestbalance')"
Rscript -e "devtools::install_github('leylabmpi/LeyLabRMisc')"

