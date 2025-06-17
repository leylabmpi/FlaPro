# Load necessary libraries
library(dplyr)
library(readr)
library(purrr)

# List of file paths (replace with snakemake@input for Snakemake workflow)
file_names <- snakemake@input

# Initialize an empty list to store DataFrames
dataframes_rcount <- list()

# Loop through file names
for (file in file_names) {
  
  # Extract sample name from the file name
  sample_name <- strsplit(basename(file), "\\.", fixed = TRUE)[[1]][1]
  
  ### COUNTS
  # Read 'Family' and 'RealCount' columns
  df_realcount <- read_delim(file, delim = "\t", col_types = cols_only(Family = col_character(), RealCount = col_double()))
  
  # Rename 'RealCount' column to the sample name
  colnames(df_realcount)[colnames(df_realcount) == "RealCount"] <- sample_name
  
  # Add the dataframe to the list
  dataframes_rcount[[length(dataframes_rcount) + 1]] <- df_realcount
}

# Merge all DataFrames on 'Family' column
merged_df_rcount <- reduce(dataframes_rcount, full_join, by = "Family")

# Fill missing values with 0
merged_df_rcount[is.na(merged_df_rcount)] <- 0

# Save the result (replace with snakemake@output for Snakemake workflow)
write_delim(merged_df_rcount, snakemake@output[[1]], delim = "\t")