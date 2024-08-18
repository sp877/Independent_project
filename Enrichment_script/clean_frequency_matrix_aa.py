import pandas as pd

# Function to clean the matrix based on amplicon ranges and to remove unwanted rows
def clean_matrix(matrix, amplicon_columns, amplicon_ranges):
    # Extract the nucleotide positions from the mutation descriptions
    positions = matrix.index.to_series().str.extract(r'(\d+)', expand=False).astype(int)
    
    # Filter out rows outside the whole range (removing everything before position 25)
    valid_positions_mask = positions >= 25
    matrix = matrix[valid_positions_mask]
    positions = positions[valid_positions_mask]

    # Create a cleaned matrix with 0 values outside of the amplicon ranges
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

# Define the columns for each amplicon in the input data
plasmid_amplicon_columns = {
    'Amplicon 1': ['13'],
    'Amplicon 2': ['14'],
    'Amplicon 3': ['15'],
    'Amplicon 4': ['16']
}

mutating_amplicon_columns = {
    'Amplicon 1': ['17'],
    'Amplicon 2': ['18'],
    'Amplicon 3': ['19'],
    'Amplicon 4': ['20']
}

treated_amplicon_columns = {
    'Amplicon 1': ['21'],
    'Amplicon 2': ['22'],
    'Amplicon 3': ['23'],
    'Amplicon 4': ['24']
}

# Set the paths for each of the input matrix
plasmid_matrix = pd.read_csv('/home/sp877/Independent_project/enchrichment_count/multiple_mutations/plasmid_multi_mutation_frequency_matrix.tsv', sep='\t', index_col=0)
mutating_cells_matrix = pd.read_csv('/home/sp877/Independent_project/enchrichment_count/multiple_mutations/mutating_cells_multi_mutation_frequency_matrix.tsv', sep='\t', index_col=0)
treated_cells_matrix = pd.read_csv('/home/sp877/Independent_project/enchrichment_count/multiple_mutations/treated_cells_multi_mutation_frequency_matrix.tsv', sep='\t', index_col=0)

# Clean the matrices
cleaned_plasmid_matrix = clean_matrix(plasmid_matrix, plasmid_amplicon_columns, amplicon_ranges)
cleaned_mutating_cells_matrix = clean_matrix(mutating_cells_matrix, mutating_amplicon_columns, amplicon_ranges)
cleaned_treated_cells_matrix = clean_matrix(treated_cells_matrix, treated_amplicon_columns, amplicon_ranges)

# Save the cleaned matrices
cleaned_plasmid_matrix.to_csv('cleaned_plasmid__multi_mutation_frequency_matrix.tsv', sep='\t')
cleaned_mutating_cells_matrix.to_csv('cleaned_mutating_cells_multi_mutation_frequency_matrix.tsv', sep='\t')
cleaned_treated_cells_matrix.to_csv('cleaned_treated_cells_multi_mutation_frequency_matrix.tsv', sep='\t')

print("Cleaned frequency matrices have been saved.")

