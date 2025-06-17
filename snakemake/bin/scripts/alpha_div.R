suppressPackageStartupMessages(library(phyloseq))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(GenomeInfoDb))

# Load the .RData file
load(snakemake@input[[1]])

#Alpha diversity

values.plot=c("#e64b35b2","#4dbbd5b2","#00a087b2","#3c5488b2","#f39b7fb2","#8491b4b2")

# alpha diversity estimation
alpha_div=mutate(estimate_richness(psq, split = TRUE, measures = c("Observed","Chao1","Shannon","Simpson")))

# save alpha diversity in table
write.table(alpha_div, file = snakemake@output[[1]], sep = "\t", quote = FALSE, row.names = FALSE)
