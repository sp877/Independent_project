import pandas as pd
from Bio import SeqIO
from collections import Counter, OrderedDict
import gzip

# Define the location of UMI in the sequence and the length of the forward index
FORWARD_INDEX_LENGTH = 5  
UMI_LENGTH = 15

# Function to extract UMI and sequence from gzipped Fastq file
def extract_umi_and_sequence(fastq_file):
    records = []  # Initialize an empty list to store UMI, sequence, and quality tuples
    # Open the gzipped FASTQ file
    with gzip.open(fastq_file, "rt") as handle:
        # Parse the FASTQ file
        for record in SeqIO.parse(handle, "fastq"):
            # Extract UMI after the forward index
            umi = str(record.seq[FORWARD_INDEX_LENGTH:FORWARD_INDEX_LENGTH + UMI_LENGTH])
            # Extract the remaining part of the sequence after the UMI
            seq = str(record.seq[FORWARD_INDEX_LENGTH + UMI_LENGTH:])
            # Extract the quality scores corresponding to the sequence (excluding forward index and UMI part)
            qual = record.letter_annotations["phred_quality"][FORWARD_INDEX_LENGTH + UMI_LENGTH:]
            # Append a tuple of (UMI, sequence, quality) to the records list
            records.append((umi, seq, qual))
            # Debug statement to print lengths and extracted values
            print(f"UMI: {umi}, Sequence Length: {len(seq)}, Quality Length: {len(qual)}")
    return records

# Function to group sequences by UMI and output them to a file
def group_by_umi(records, output_file, discarded_output_file, min_family_size=2):
    # Ordered dictionary to maintain the order of insertion. store grouped sequences by UMI
    grouped_records = OrderedDict()  
    # Ordered dictionary to store discarded records. Store discared sequences (that do not meet min family size)
    discarded_records = OrderedDict()  

    # Group the records by UMI. Iterates over each UMI, sequence & quality tuple in 'records' list
    for umi, sequence, quality in records:
        # Checks if UMI not already in grouped_records dictionary
        if umi not in grouped_records:
            # If UMI is new, puts empty list for that UMI in grouped_records
            grouped_records[umi] = []
        # Adds sequence and quality ;scored to the list for that UMI
        grouped_records[umi].append((sequence, quality))

    # Open the output files for writing the grouped and discared 
    with open(output_file, 'w') as f, open(discarded_output_file, 'w') as df:
        # Iterates over each UMI and its corresponding group of sequences
        for umi, group in grouped_records.items():
            # Checks if group size meets min family size or not
            if len(group) >= min_family_size:
                # Write the UMI as a header in the output file
                f.write(f">{umi}\n")
                #Iterates over each sequence and quality in group
                for sequence, quality in group:
                    # Convert quality scores to ASCII characters
                    quality_str = ''.join(map(chr, [q + 33 for q in quality]))
                    # Write the sequence and corresponding quality scores to the file
                    f.write(f"{sequence}\n{quality_str}\n")
                    # Debug statement to print sequence and quality
                    print(f"Group UMI: {umi}, Sequence: {sequence}, Quality: {quality_str}")
            # If group size smaller than the min fam size, sequences are discarded
            else:
                # Adds group to discarded_records dictionary
                discarded_records[umi] = group
                # Write the discarded UMI group to the discarded output file
                df.write(f">{umi}\n")
                for sequence, quality in group:
                    quality_str = ''.join(map(chr, [q + 33 for q in quality]))
                    df.write(f"{sequence}\n{quality_str}\n")
                    print(f"Discarded UMI: {umi}, Sequence: {sequence}, Quality: {quality_str}")

    return output_file, grouped_records, discarded_output_file

# Function to find the consensus sequence for a group of sequences
def find_consensus_sequence(sequences, qualities, threshold=1.0, am_threshold=0.4):
    if not sequences:
        return None, None

    # Determine the maximum length of the sequences in the group
    seq_len = max(len(seq) for seq in sequences)  
    # List to store the consensus sequence
    consensus_seq = []
    # List to store the consensus quality scores  
    consensus_qual = []  



    for i in range(seq_len):
        # Collect bases at position i from all sequences
        bases = [seq[i] for seq in sequences if i < len(seq)]
        # Count the frequency of each base at position i  
        base_count = Counter(bases) 
        # Find the most common base and its count 
        most_common_base, count = base_count.most_common(1)[0]  
         # Check if the most common base meets the consensus threshold
        if count / len(sequences) >= threshold:  
            consensus_seq.append(most_common_base)
            # Calculate the average quality score at position i
            avg_quality = int(sum(qual[i] for qual in qualities if i < len(qual)) / len(qualities))
            consensus_qual.append(avg_quality)
        else:
            # The most common base does not meet the threshold
            # Handle artifactual mutations
            am_fraction = (len(sequences) - count) / len(sequences)
            if am_fraction < am_threshold:
                # The proportion of sequences with a different base is less than am_threshold
                consensus_seq.append(most_common_base)
                avg_quality = int(sum(qual[i] for qual in qualities if i < len(qual)) / len(qualities))
                consensus_qual.append(avg_quality)
            else:
                # Too much disagreement at this position, abandon forming consensus for this UMI group
                return None, None

    return ''.join(consensus_seq), consensus_qual

# Function to process UMI groups and find consensus sequences
def process_umi_groups(grouped_records, consensus_output_file, min_family_size=2, consensus_threshold=1.0, am_threshold=0.4):
    # List to store consensus sequences
    consensus_sequences = []  
    # 
    for umi, sequences_quals in grouped_records.items():
        # Ensure UMI groups with at least min_family_size
        if len(sequences_quals) >= min_family_size:  
            # Extract sequences from the group
            sequences = [seq for seq, _ in sequences_quals] 
            # Extract quality scores from the group 
            qualities = [qual for _, qual in sequences_quals]  
            print(f"Processing UMI: {umi} with {len(sequences)} sequences")
            consensus_seq, consensus_qual = find_consensus_sequence(sequences, qualities, consensus_threshold, am_threshold)
            if consensus_seq:
                print(f"Consensus sequence for {umi}: {consensus_seq}")
                consensus_sequences.append((umi, len(sequences), consensus_seq, consensus_qual))
            else:
                print(f"No consensus sequence for {umi}")

    # Open the consensus output file for writing    
    with open(consensus_output_file, 'w') as output_handle:
        for umi, seq_count, seq, qual in consensus_sequences:
            output_handle.write(f"@{umi}_{seq_count}\n{seq}\n+\n{''.join(chr(q + 33) for q in qual)}\n")

    print(f"Consensus sequences written to {consensus_output_file}")

    return consensus_sequences

# Main function
def main(fastq_file, grouped_output_file, discarded_output_file, consensus_output_file, min_family_size=2, consensus_threshold=1.0, am_threshold=0.4):
    # Extract UMIs and sequences from the input FASTQ file
    records = extract_umi_and_sequence(fastq_file)
    # Group sequences by UMI and write them to an output file
    grouped_file, grouped_records, discarded_file = group_by_umi(records, grouped_output_file, discarded_output_file, min_family_size)
    # Process UMI groups to find consensus sequences
    process_umi_groups(grouped_records, consensus_output_file, min_family_size, consensus_threshold, am_threshold)

# Define input and output files
fastq_file = "Unknown_BT944-ZX01-017_merged.fq.gz"
grouped_output_file = "grouped_umis_017.txt"
discarded_output_file = "discarded_umis_017.txt"
consensus_output_file = "consensus_sequences_017.fastq"

# Run the main function
main(fastq_file, grouped_output_file, discarded_output_file, consensus_output_file, min_family_size=2, consensus_threshold=1.0, am_threshold=0.4)

