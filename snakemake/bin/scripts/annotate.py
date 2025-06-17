import pandas as pd

# Load the dataframes
tax_df = pd.read_csv(snakemake.input['taxonomy'], sep='\t')
phen_df = pd.read_csv(snakemake.input['phenotype'], sep=',', header=None, names=['Flagellin_ID'])
counts_df = pd.read_csv(snakemake.input['merged'], sep='\t')

# Rename the column in counts_df for consistent merging
counts_df.rename(columns={"Family": "Flagellin_ID"}, inplace=True)

# Identify potentially silent genes
fla_set = set(phen_df['Flagellin_ID'])

# Directly annotate 'Phenotype' in phen_df based on fla_set
phen_df['Phenotype'] = phen_df['Flagellin_ID'].apply(lambda x: 'potentially silent' if x in fla_set else 'non-silent')

# Merge tax_df with phen_df
# Note: Ensure that 'how='left'' to keep all records from tax_df
merged_df = pd.merge(tax_df, phen_df, on='Flagellin_ID', how='left')

# For flagellin IDs not matched in phen_df, mark as 'non-silent'
# Adjusted to avoid the FutureWarning
merged_df['Phenotype'] = merged_df['Phenotype'].fillna('non-silent')

# Merge the merged_phen_tax with counts_df
final_merged_df = pd.merge(merged_df, counts_df, on='Flagellin_ID', how='right')

# Save the final merged dataframe
final_merged_df.to_csv(snakemake.output[0], sep='\t', index=False)
