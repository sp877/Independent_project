import pandas as pd
import numpy as np
import re

# Directory paths for ASMV outputs
base_path = '/home/sp877/Independent_project/enchrichment_count/ASMV_output_noumi'

# Function to read and process the dataset
def read_and_process_data(file_path):
    try:
        # Read file line by line
        with open(file_path, 'r') as file:
            lines = file.readlines()

        # Process each line
        # creates empty list to store processed lines
        processed_lines = []
        # Iterate over each line in the file
        for line in lines:
            split_line = line.strip().split('\t')
            # Checks if the line has less than 9 columns
            if len(split_line) < 9:
                # Adds empty strings to ensure the line has exactly 9 columns
                split_line.extend([''] * (9 - len(split_line)))  
            # Adds processed line to list of processed lines
            processed_lines.append(split_line)

        # Create DataFrame
        data = pd.DataFrame(processed_lines)
        return data

    # Finds exceptions that occur during file reading
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return pd.DataFrame()

# Function to calculate the mutation frequency ratio
def calculate_ratios(data, max_nucleotide_changes):
    mutation_data = {}

    # Iterates over each row in dataframe    
    for index, row in data.iterrows():
        if len(row) < 9:
            continue
        
        try:
            # Total number of reads (second column). converts to float
            total_reads = float(row[1])
            # Extracts nucleotide changes from 5th column
            nucleotide_changes = row[4] 
            # Extracts aa changes from 9th column
            amino_acid_changes = row[8] 
 
            # Checks if nucleotide_changes is a non-empty string         
            if isinstance(nucleotide_changes, str) and nucleotide_changes.strip():
                #Splits nucleotide changes into a list
                nucleotide_list = nucleotide_changes.split(',')
                # Continues if number of nucleotide changes does not exceed the maximum
                if len(nucleotide_list) <= max_nucleotide_changes:
                    # Checks if amino_acid_changes is a non-empty string
                    if isinstance(amino_acid_changes, str) and amino_acid_changes.strip():
                        # Slpits amino acid changes into a list and remove whitespace
                        mutation_list = [mutation.strip() for mutation in amino_acid_changes.split(';')]
                        for mutation in mutation_list:
                            # Checks if the mutation is already in the mutation_data dictionary
                            if mutation not in mutation_data:
                                # Created dictionary new entry in mutationdata disctionary for specific mutation, but only if an entry for that mutation does not already exist. 
                                mutation_data[mutation] = {'ratios': []}
                            try:
                                # Converts first column to a float, num. of reads for mutation
                                num_reads = float(row[0])
                                # calculates ratio
                                ratio = num_reads / total_reads
                                # Add ratio in list for mutation
                                mutation_data[mutation]['ratios'].append(ratio)
                            # Catches division by 0 errors 
                            except ZeroDivisionError:
                                continue
                            except Exception as e:
                                continue
                else:
                    continue
        
        except ValueError as e:
            continue
    
    # Sum the ratios for each mutation if there are multiple occurrences and stores in a dictionary
    ratios = {mutation: sum(values['ratios']) for mutation, values in mutation_data.items()}
    return ratios

# Function to extract numeric part from mutation for sorting
def extract_numeric_part(mutation):
    match = re.search(r'\d+', mutation)
    return int(match.group()) if match else float('inf')

# Define groups with specific file sets
groups = {
    "plasmid": [13, 14, 15, 16],
    "mutating_cells": [17, 18, 19, 20],
    "treated_cells": [21, 22, 23, 24]
}

# Initialize global dictionaries for enrichment data
global_enrichment_dict_single = {}
global_enrichment_dict_multi = {}

# Process each group independently
for group_name, file_set in groups.items():
    
    single_dict = {}
    multi_dict = {}
    
    # Initialize sets for each group
    group_mutations_single = set()
    group_mutations_multi = set()
    
    for file_num in file_set:
        file_path = f"{base_path}/ASM_{file_num}/{file_num}.variantCounts"
        
        # Read and process data
        data = read_and_process_data(file_path)
        
        # Skip if data couldn't be read
        if data.empty:
            print(f"Skipping file {file_num} due to empty data.")
            continue
        
        # Calculate ratios
        ratios_single = calculate_ratios(data, 3)  # Consider up to 3 nucleotide changes
        ratios_multi = calculate_ratios(data, 3)
        
        # Merge mutations into the group-specific mutation sets
        group_mutations_single.update(ratios_single.keys())
        group_mutations_multi.update(ratios_multi.keys())
        
        # Store frequency ratios for the current file
        single_dict[file_num] = {mutation: ratios_single.get(mutation, 0) for mutation in group_mutations_single}
        multi_dict[file_num] = {mutation: ratios_multi.get(mutation, 0) for mutation in group_mutations_multi}
    
    # Create and save frequency matrices for the group
    single_matrix = pd.DataFrame(single_dict).reindex(index=sorted(group_mutations_single, key=extract_numeric_part)).fillna(0)
    multi_matrix = pd.DataFrame(multi_dict).reindex(index=sorted(group_mutations_multi, key=extract_numeric_part)).fillna(0)
    
    single_matrix.to_csv(f'{group_name}_single_mutation_frequency_matrix_aa_nucleo_fixed_noumi.tsv', sep='\t')
    multi_matrix.to_csv(f'{group_name}_multi_mutation_frequency_matrix_aa_nucleo_fixed_noumi.tsv', sep='\t')

# Now for enrichment calculations, using comparisons across groups
comparisons = {
    "plasmid_vs_treated": [(13, 21), (14, 22), (15, 23), (16, 24)],
    "plasmid_vs_mutating": [(13, 17), (14, 18), (15, 19), (16, 20)],
    "mutating_vs_treated": [(17, 21), (18, 22), (19, 23), (20, 24)]
}

# Process each comparison independently for enrichment
for comp_name, comp_pairs in comparisons.items():
    for start_file, end_file in comp_pairs:
        start_path = f"{base_path}/ASM_{start_file}/{start_file}.variantCounts"
        end_path = f"{base_path}/ASM_{end_file}/{end_file}.variantCounts"
        
        # Read and process data
        start_data = read_and_process_data(start_path)
        end_data = read_and_process_data(end_path)
        
        # Skip if data couldn't be read
        if start_data.empty or end_data.empty:
            print(f"Skipping comparison {start_file} vs {end_file} due to empty data.")
            continue
        
        # Calculate ratios
        start_ratios_single = calculate_ratios(start_data, 3)  # Consider up to 3 nucleotide changes
        start_ratios_multi = calculate_ratios(start_data, 3)
        end_ratios_single = calculate_ratios(end_data, 3)
        end_ratios_multi = calculate_ratios(end_data, 3)
        
        # Calculate enrichment values with proper handling
        enrichment_values_single = {
            mutation: (
                end_ratios_single.get(mutation, np.nan) / start_ratios_single.get(mutation, np.nan)
                if start_ratios_single.get(mutation, np.nan) > 0 else (999999 if end_ratios_single.get(mutation, 0) > 0 else np.nan)
            )
            for mutation in set(start_ratios_single.keys()).union(set(end_ratios_single.keys()))
        }
        enrichment_values_multi = {
            mutation: (
                end_ratios_multi.get(mutation, np.nan) / start_ratios_multi.get(mutation, np.nan)
                if start_ratios_multi.get(mutation, np.nan) > 0 else (999999 if end_ratios_multi.get(mutation, 0) > 0 else np.nan)
            )
            for mutation in set(start_ratios_multi.keys()).union(set(end_ratios_multi.keys()))
        }
        
        # Add to global enrichment dictionary
        comparison_label = f"{start_file}vs{end_file}"
        global_enrichment_dict_single[comparison_label] = enrichment_values_single
        global_enrichment_dict_multi[comparison_label] = enrichment_values_multi

# Convert global enrichment dictionaries to DataFrames
enrichment_matrix_single = pd.DataFrame.from_dict(global_enrichment_dict_single, orient='index').T
enrichment_matrix_multi = pd.DataFrame.from_dict(global_enrichment_dict_multi, orient='index').T

# Define the desired column order
desired_order = [
    '13vs17', '14vs18', '15vs19', '16vs20',
    '13vs21', '14vs22', '15vs23', '16vs24',
    '17vs21', '18vs22', '19vs23', '20vs24'
]

# Reorder columns in the matrices according to the desired order
enrichment_matrix_single = enrichment_matrix_single[desired_order]
enrichment_matrix_multi = enrichment_matrix_multi[desired_order]

# Ensure the index is sorted by the numeric part in the mutation (if present)
enrichment_matrix_single.sort_index(key=lambda x: [extract_numeric_part(m) for m in x], inplace=True)
enrichment_matrix_multi.sort_index(key=lambda x: [extract_numeric_part(m) for m in x], inplace=True)

# Replace np.nan with 0 where applicable
enrichment_matrix_single.replace(np.nan, 0, inplace=True)
enrichment_matrix_multi.replace(np.nan, 0, inplace=True)

# Save the global enrichment matrices
enrichment_matrix_single.to_csv('all_single_mutation_enrichment_matrix_aa_nucleo_fixed_noumi.tsv', sep='\t')
enrichment_matrix_multi.to_csv('all_multi_mutation_enrichment_matrix_aa_nucleo_fixed_noumi.tsv', sep='\t')

print("Frequency and enrichment matrices have been saved.")

