# Load necessary libraries
suppressPackageStartupMessages(library(dplyr))

# Read the input file (using snakemake's input directive)
df <- read.delim(snakemake@input[[1]], sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Check if the necessary columns 'Hits' and 'Count' exist
if (!all(c('Hits', 'Count') %in% colnames(df))) {
  stop("The input file does not contain 'Hits' and/or 'Count' columns.")
}

# Create the 'RealCount' column if it does not exist
if (!"RealCount" %in% colnames(df)) {
  df$RealCount <- 0
}

# Apply the condition to update the 'RealCount' column
df <- df %>% mutate(RealCount = ifelse(Count > 0, Hits, 0))

# Save the updated dataframe to the output file (using snakemake's output directive)
write.table(df, snakemake@output[[1]], sep = "\t", row.names = FALSE, quote = FALSE)