#!/usr/bin/env Rscript

# Load necessary libraries
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(purrr))

# Command-line arguments
args <- commandArgs(trailingOnly = TRUE)
input_files <- args[1:(length(args)-1)]  # All input files
output_file <- args[length(args)]  # Output file path

# Function to process and merge files
process_and_merge <- function(input_files, output_file) {
  
  # Read and process each file
  dataframes <- map(input_files, function(file) {
    
    # Extract sample name from file name (without extension)
    sample_name <- tools::file_path_sans_ext(basename(file))
    
    # Read 'Family' and 'RealCount' columns
    df <- read_delim(file, delim = "\t", col_types = cols_only(Family = col_character(), RealCount = col_double()))
    
    # Rename 'RealCount' to the sample name
    colnames(df)[colnames(df) == "RealCount"] <- sample_name
    
    return(df)
  })
  
  # Merge all dataframes on 'Family'
  merged_df <- reduce(dataframes, full_join, by = "Family")
  
  # Fill missing values with 0
  merged_df[is.na(merged_df)] <- 0
  
  # Write the result to the output file
  write_delim(merged_df, output_file, delim = "\t")
}

# Call the function to process files and merge
process_and_merge(input_files, output_file)