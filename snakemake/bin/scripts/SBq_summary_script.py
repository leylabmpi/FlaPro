import pandas as pd

# Folder path

# List of file paths
file_names = snakemake.input

# Initialize an empty list to store DataFrames
dataframes_rcount = []

for file in file_names:
    # Extract sample name from the file name
    sample_name = file.split('/')[-1].split('.', 1)[0]
    
    ### COUNTS
    # Read 'Family' and 'Count' columns
    df_realcount = pd.read_csv(file, sep="\t", usecols=['Family', 'RealCount'])
    # Rename 'Count' column to indicate the sample name
    df_realcount.rename(columns={'RealCount': f'{sample_name}'}, inplace=True)
    # Append the DataFrame to the list
    dataframes_rcount.append(df_realcount)
    
    

# Merge all DataFrames on 'Family' - COUNT
merged_df_rcount = dataframes_rcount[1]
for df_rcount in dataframes_rcount[1:]:
    merged_df_rcount = pd.merge(merged_df_rcount, df_rcount, on='Family', how='outer')

# Fill missing values with 0 (assuming missing counts should be 0)
merged_df_rcount.fillna(0, inplace=True)

# Save
merged_df_rcount.to_csv(snakemake.output[0], sep='\t', index=False)