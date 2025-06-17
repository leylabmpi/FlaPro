import pandas as pd

# Read the input file specified in the Snakemake rule
df = pd.read_csv(snakemake.input[0], sep='\t')

# Check if the necessary columns exist
if 'Hits' not in df.columns or 'Count' not in df.columns:
    raise KeyError("The input file does not contain 'Hits' and/or 'Count' columns.")

# Create the 'RealCount' column if it does not exist
if 'RealCount' not in df.columns:
    df['RealCount'] = 0

# Apply the condition to create/update the 'RealCount' column
df['RealCount'] = df.apply(lambda x: x['Hits'] if x['Count'] > 0 else 0, axis=1)

# Save the updated dataframe to the output file specified in the Snakemake rule
df.to_csv(snakemake.output[0], sep='\t', index=False)
