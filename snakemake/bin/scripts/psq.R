
suppressPackageStartupMessages(library(phyloseq))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(dplyr))


#Importing counts
counts_df = read.table(snakemake@input[[1]], header = TRUE, sep = "\t")
print("Counts read.")

#Importing taxonomy
taxonomy_df = read.table(snakemake@input[[2]], sep = "\t", header = TRUE, fill = TRUE)
print("Taxonomy imported.")

#Importing metadata
#metadata_df = read.table(snakemake@input[[3]], sep = '\t', header = TRUE)

# Convert the dataframe to a tibble
counts_df <- as_tibble(counts_df)
taxonomy_df <- as_tibble(taxonomy_df)
#metadata_df <- as_tibble(metadata_df)
print("Converted to tibbles.")

# Change rownames
otu_table <- counts_df %>% 
  tibble::column_to_rownames(var = "Family")
print('Change rownames')

tax_table <- taxonomy_df %>% 
  tibble::column_to_rownames(var = "Flagellin_ID")

print("Rownames changed.")

#sample_table <- metadata_df %>% 
#  tibble::column_to_rownames(var = "SRR")

  # Convert to matrix for phyloseq
otu_table = as.matrix(otu_table)
tax_table = as.matrix(tax_table)

print("Converted to matrices.")
print(dim(otu_table))
print(dim(tax_table))

# Create phyloseq object
OTU = otu_table(otu_table, taxa_are_rows = TRUE)
TAX = tax_table(tax_table)
#samples = sample_data(sample_table)
psq = phyloseq(OTU, TAX)#, samples)

print("Created phyloseq object.")

psq <- prune_taxa(taxa_sums(psq) > 0, psq)
print("Taxa pruned.")

#psq <- subset_samples(psq, Diagnosis %in% c("HC", "CD"))
# Save the phyloseq object
print("Phyloseq object:")
print(psq)

save(psq, file = snakemake@output[[1]])
print("Phyloseq object saved.")