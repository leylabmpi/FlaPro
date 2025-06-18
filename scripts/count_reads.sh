#!/bin/bash

# Define the input tab-separated file
input_file=$1 #"samples.txt"

# Skip the first line (header) and process each line
tail -n +2 "$input_file" | while IFS=$'\t' read -r Sample Read1 Read2
do
  # Count the number of lines in Read1 and Read2
  count_read1=$(zcat "$Read1" | wc -l)

  # Divide the line counts by 4
  count_read1=$((count_read1 / 4))
  
  # Output the results
  echo "$Sample	$count_read1"
  
done

