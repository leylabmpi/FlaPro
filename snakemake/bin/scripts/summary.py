import pandas as pd
import sys

def process_real_counts(file_names, output_file):
    # Initialize an empty list to store DataFrames
    dataframes_rcount = []

    for file in file_names:
        # Extract sample name from the file name
        sample_name = file.split('/')[-1].split('.', 1)[0]
        
        ### COUNTS
        # Read 'Family' and 'RealCount' columns
        df_realcount = pd.read_csv(file, sep="\t", usecols=['Family', 'RealCount'])
        # Rename 'RealCount' column to indicate the sample name
        df_realcount.rename(columns={'RealCount': sample_name}, inplace=True)
        # Append the DataFrame to the list
        dataframes_rcount.append(df_realcount)

    # Merge all DataFrames on 'Family'
    merged_df_rcount = dataframes_rcount[0]
    for df_rcount in dataframes_rcount[1:]:
        merged_df_rcount = pd.merge(merged_df_rcount, df_rcount, on='Family', how='outer')

    # Fill missing values with 0 (assuming missing counts should be 0)
    merged_df_rcount.fillna(0, inplace=True)

    # Save
    merged_df_rcount.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    # Extract command-line arguments
    file_names = sys.argv[1:-1]  # all arguments except the last one are input files
    output_file = sys.argv[-1]  # last argument is the output file

    # Call the function with the provided command line arguments
    process_real_counts(file_names, output_file)