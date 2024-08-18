import pandas as pd
import numpy as np

# Set path for ASMV outputs 
base_path = '/home/sp877/Independent_project/enchrichment_count/ASMV_output_noumi'

# Function to read and process the dataset, handling errors 
def read_and_process_data(file_path):
    try:
        data = pd.read_csv(file_path, sep='\t', header=None, on_bad_lines='skip')
        # Returns the loaded data as a DataFrame
        return data
    # Handles any exceptions that occur during the file reading process
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return pd.DataFrame()

# Function to calculate the mutation frequency ratio
def calculate_ratios(data, max_mutations):
    # Creates an empty dictionary to store the mutation info
    mutation_data = {}

    # Iterates over each row of the DataFrame    
    for index, row in data.iterrows():
        # Extracts total number of reads (second column)
        total_reads = row[1] 
        # Extracts mutations from 5th column (index 4) 
        mutations = row[4]  

        # Checks if mutation column contains a string        
        if isinstance(mutations, str):
            # Splits string into individual mutations and strips any surrounding whitespace
            mutation_list = [mutation.strip() for mutation in mutations.split(',')]
            # Checks if the number of mutations does not exceed the maximum allowed mutations
            if len(mutation_list) <= max_mutations:
                # Iterates over each mutation in the list
                for mutation in mutation_list:
                    # Checks if the mutation is not already in the mutation_data dictionary
                    if mutation not in mutation_data:
                        # create a dictionary with an empty list for storing ratios if the mutation is new
                        mutation_data[mutation] = {'ratios': []}
                    # Extracts the number of reads from the first column
                    num_reads = row[0]  
                    # Calculates the ratio of reads for the mutation to the total reads
                    ratio = num_reads / total_reads
                    # Adds calculated ratio to the list for the mutation
                    mutation_data[mutation]['ratios'].append(ratio)
    
    # Sum the ratios for each mutation if there are multiple ratios for same mutation
    ratios = {mutation: sum(values['ratios']) for mutation, values in mutation_data.items()}
    return ratios

# Define groups with specific file sets
groups = {
    "plasmid": [13, 14, 15, 16],
    "mutating_cells": [17, 18, 19, 20],
    "treated_cells": [21, 22, 23, 24]
}

# Create global dictionaries for enrichment data
# creates an empty dictionary to store single mutation enrichment data across all comparisons
global_enrichment_dict_single = {}
# creates empty dictionary to store multi-mutation enrichment data across all comparisons
global_enrichment_dict_multi = {}

# Process each group independently
# Iterates over each group defined in the groups dictionary
for group_name, file_set in groups.items():

    # creates empty dictionary to store single mutation ratios for the current group    
    single_dict = {}
    # creates empty dictionary to store multi-mutation ratios for the current group
    multi_dict = {}
    
    # create sets for each group
    group_mutations_single = set()
    group_mutations_multi = set()
    
    #  Iterates over each file number in the current group
    for file_num in file_set:
        file_path = f"{base_path}/ASM_{file_num}/{file_num}.variantCounts"
        
        # Read and process data
        data = read_and_process_data(file_path)
        
        # Skip if data couldn't be read
        if data.empty:
            print(f"Skipping file {file_num} due to empty data.")
            continue
        
        # Calculate ratios
        ratios_single = calculate_ratios(data, 1)
        ratios_multi = calculate_ratios(data, 3)
        
        # Merge mutations into the group-specific mutation sets
        group_mutations_single.update(ratios_single.keys())
        group_mutations_multi.update(ratios_multi.keys())
        
        # Store frequency ratios for the current file
        single_dict[file_num] = {mutation: ratios_single.get(mutation, 0) for mutation in group_mutations_single}
        multi_dict[file_num] = {mutation: ratios_multi.get(mutation, 0) for mutation in group_mutations_multi}
    
    # Create and save frequency matrices for the group
    single_matrix = pd.DataFrame(single_dict).reindex(index=sorted(group_mutations_single)).fillna(0)
    multi_matrix = pd.DataFrame(multi_dict).reindex(index=sorted(group_mutations_multi)).fillna(0)
    
    single_matrix.to_csv(f'{group_name}_single_mutation_frequency_matrix_fixed_noumi.tsv', sep='\t')
    multi_matrix.to_csv(f'{group_name}_multi_mutation_frequency_matrix_fixed_noumi.tsv', sep='\t')

# Set and define the comparison groups for enrichment calculation
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
        
        # Skip if data cant be read
        if start_data.empty or end_data.empty:
            print(f"Skipping comparison {start_file} vs {end_file} due to empty data.")
            continue
        
        # Calculate ratios
        start_ratios_single = calculate_ratios(start_data, 1)
        start_ratios_multi = calculate_ratios(start_data, 3)
        end_ratios_single = calculate_ratios(end_data, 1)
        end_ratios_multi = calculate_ratios(end_data, 3)
        
        # Calculate enrichment values with proper handling
        enrichment_values_single = {
            # Calculates the enrichment for each mutation by dividing the end ratio by the start ratio. Handles cases where the start ratio is zero or missing by assigning a large default value (999999) if the end ratio is non-zero.
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

# Ensure the index is sorted by position and mutation type
enrichment_matrix_single.sort_index(inplace=True)
enrichment_matrix_multi.sort_index(inplace=True)

# Replace np.nan with 0 
enrichment_matrix_single.replace(np.nan, 0, inplace=True)
enrichment_matrix_multi.replace(np.nan, 0, inplace=True)

# Save the global enrichment matrices
enrichment_matrix_single.to_csv('all_single_mutation_enrichment_matrix_fixed_noumi.tsv', sep='\t')
enrichment_matrix_multi.to_csv('all_multi_mutation_enrichment_matrix_fixed_noumi.tsv', sep='\t')

print("Frequency and enrichment matrices have been saved.")






