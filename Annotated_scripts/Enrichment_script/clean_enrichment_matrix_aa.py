import pandas as pd

# Function to clean the matrix based on amplicon ranges and remove unwanted rows
def clean_matrix(matrix, amplicon_columns, amplicon_ranges):
    # Extract the nucleotide positions from the mutation descriptions
    positions = matrix.index.to_series().str.extract(r'(\d+)', expand=False).astype(int)
    
    # Filter out rows outside the global range (removing everything before position 25)
    valid_positions_mask = positions >= 25
    matrix = matrix[valid_positions_mask]
    positions = positions[valid_positions_mask]

    # Create a cleaned matrix with zeroed values outside the amplicon-specific ranges
    cleaned_matrix = matrix.copy()
    for amplicon, (start, end) in amplicon_ranges.items():
        for col in amplicon_columns[amplicon]:
            mask = positions.between(start, end, inclusive='both')
            cleaned_matrix[col] = matrix[col].where(mask, 0)
    
    return cleaned_matrix

# Define the ranges for each amplicon
amplicon_ranges = {
    'Amplicon 1': (25, 81),
    'Amplicon 2': (82, 137),
    'Amplicon 3': (138, 195),
    'Amplicon 4': (196, 247)
}

# Define the columns for each amplicon in the new matrix
all_enrichment_amplicon_columns = {
    'Amplicon 1': ['13vs17', '13vs21', '17vs21'],
    'Amplicon 2': ['14vs18', '14vs22', '18vs22'],
    'Amplicon 3': ['15vs19', '15vs23', '19vs23'],
    'Amplicon 4': ['16vs20', '16vs24', '20vs24']
}

# Load the all enrichment matrix
all_enrichment_matrix = pd.read_csv('/home/sp877/Independent_project/enchrichment_count/multiple_mutations/all_multi_mutation_enrichment_matrix.tsv', sep='\t', index_col=0)

# Clean the matrix
cleaned_all_enrichment_matrix = clean_matrix(all_enrichment_matrix, all_enrichment_amplicon_columns, amplicon_ranges)

# Save the cleaned matrix
cleaned_all_enrichment_matrix.to_csv('cleaned_all_multi_mutation_enrichment_matrix.tsv', sep='\t')

print("Cleaned all enrichment matrix has been saved.")

