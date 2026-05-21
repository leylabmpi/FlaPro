suppressPackageStartupMessages(library(phyloseq))
suppressPackageStartupMessages(library(dplyr))

# Load the .RData file
load(snakemake@input[[1]])

#Alpha diversity

# alpha diversity estimation
alpha_div=mutate(estimate_richness(psq, split = TRUE, measures = c("Observed","Chao1","Shannon","Simpson")))

# save alpha diversity in table
write.table(alpha_div, file = snakemake@output[[1]], sep = "\t", quote = FALSE, row.names = FALSE)
