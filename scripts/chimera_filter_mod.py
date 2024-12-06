#############################################################################################################################
#                                     Candidate Chimera Filtering (Optimized with Reverse Complement for Negative Strand)
#############################################################################################################################

import os
import math
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
import time

start_time = time.time()

# Parse arguments
parser = argparse.ArgumentParser(description='A list of sequences from a CSV file.')
parser.add_argument('input', help='fusionCandidates')
parser.add_argument('JDis', help='JunctionDistance')
parser.add_argument('Sent', help='ShannonEntropy')
parser.add_argument('DNAseq', help='dnatopLevel')
parser.add_argument('bptSeq', help='Sequence around the bpt, value of JunctionSeq from config file')
parser.add_argument('bptDinuc')
parser.add_argument('count1', help='Top rows for the same transcript pair')
parser.add_argument('promis', help='PROMISCIOUS')
parser.add_argument('split', help='splitCount')
parser.add_argument('span', help='spanCount')
parser.add_argument('output', help='finalFusion')
args = parser.parse_args()

output_dir = os.path.dirname(args.output)

# Load data
filtered_table = pd.read_csv(args.input, delimiter='\t')
filtered_table.rename(columns={'0': 'Split_support_fbpt'}, inplace=True)

# Filter based on Junction Distance
filtered_table = filtered_table[abs(filtered_table['new_bpT1'] - filtered_table['new_bpT2']) >= int(args.JDis)]
dummy_fusion = filtered_table.copy()

# Define conditions and apply selections for strand-based coordinates
sent = int(args.Sent)
conditions = [
    (filtered_table['strand_x'] == '+') & (filtered_table['strand_y'] == '+'),
    (filtered_table['strand_x'] == '-') & (filtered_table['strand_y'] == '-'),
    (filtered_table['strand_x'] == '+') & (filtered_table['strand_y'] == '-'),
    (filtered_table['strand_x'] == '-') & (filtered_table['strand_y'] == '+')
]

dummy_fusion['start1'] = np.select(conditions, [
    (dummy_fusion['new_bpT1'] - sent)+1,
    dummy_fusion['new_bpT1']-1,
    (dummy_fusion['new_bpT1'] - sent)+1,
    dummy_fusion['new_bpT1']-1
])

dummy_fusion['end1'] = np.select(conditions, [
    dummy_fusion['new_bpT1']+1,
    (dummy_fusion['new_bpT1'] + sent)-1,
    dummy_fusion['new_bpT1']+1,
    (dummy_fusion['new_bpT1'] + sent)-1
])

dummy_fusion['start2'] = np.select(conditions, [
    dummy_fusion['new_bpT2']-1,
    (dummy_fusion['new_bpT2'] - sent)+1,
    (dummy_fusion['new_bpT2'] - sent)+1,
    dummy_fusion['new_bpT2']-1
])

dummy_fusion['end2'] = np.select(conditions, [
    (dummy_fusion['new_bpT2'] + sent)-1,
    dummy_fusion['new_bpT2']+1,
    dummy_fusion['new_bpT2']+1,
    (dummy_fusion['new_bpT2'] + sent)-1
])

# Load entire DNA sequence file into memory
sequence_dict = {record.id: str(record.seq) for record in SeqIO.parse(args.DNAseq, "fasta")}

# Calculate entropy with reverse complement handling for negative strands
def calculate_entropy(sequence):
    if sequence:
        count = Counter(sequence)
        total_bases = len(sequence)
        probabilities = [count[base] / total_bases for base in 'ATCG' if base in count]
        return -sum(p * math.log2(p) for p in probabilities if p > 0)
    return 0

# Fetch sequence for each breakpoint, reverse-complement if negative strand, and calculate entropy
left_entropies = []
right_entropies = []

for _, row in dummy_fusion.iterrows():
    left_seq = sequence_dict.get(row['chr_x'])[int(row['start1']):int(row['end1'])] if row['chr_x'] in sequence_dict else ""
    right_seq = sequence_dict.get(row['chr_y'])[int(row['start2']):int(row['end2'])] if row['chr_y'] in sequence_dict else ""

    # Reverse complement if strand is negative
    if row['strand_x'] == '-':
        left_seq = str(Seq(left_seq).reverse_complement())
    if row['strand_y'] == '-':
        right_seq = str(Seq(right_seq).reverse_complement())

    # Calculate entropy for each sequence
    left_entropies.append(calculate_entropy(left_seq))
    right_entropies.append(calculate_entropy(right_seq))

dummy_fusion['LeftBreakEntropy'] = np.round(left_entropies, 2)
dummy_fusion['RightBreakEntropy'] = np.round(right_entropies, 2)

# Drop unnecessary columns and save the result
dummy_fusion.drop(['end2', 'start2', 'end1', 'start1'], axis=1, inplace=True)
#dummy_fusion.to_csv(args.output, sep='\t', index=False)

#print(dummy_fusion)
print("Entropy calculation is completed.")

print("Script completed in", time.time() - start_time, "seconds.")


## SPLICE PATTERN  
import pandas as pd
import numpy as np
import os
from Bio.Seq import Seq
import subprocess

# Define conditions for each combination of strand_x and strand_y
conditions = [
    (filtered_table['strand_x'] == '+') & (filtered_table['strand_y'] == '+'),
    (filtered_table['strand_x'] == '-') & (filtered_table['strand_y'] == '-'),
    (filtered_table['strand_x'] == '+') & (filtered_table['strand_y'] == '-'),
    (filtered_table['strand_x'] == '-') & (filtered_table['strand_y'] == '+')
]

# Define corresponding start1, end1, start2, end2 values for each condition
start1_values = [
    filtered_table['new_bpT1'],
    filtered_table['new_bpT1'] - int(args.bptDinuc),
    filtered_table['new_bpT1'],
    filtered_table['new_bpT1'] - int(args.bptDinuc) 
]

end1_values = [
    filtered_table['new_bpT1'] + int(args.bptDinuc),
    filtered_table['new_bpT1'],
    filtered_table['new_bpT1'] + int(args.bptDinuc),
    filtered_table['new_bpT1']
]

start2_values = [
    filtered_table['new_bpT2'] - int(args.bptDinuc),
    filtered_table['new_bpT2'],
    filtered_table['new_bpT2'],
    filtered_table['new_bpT2'] - int(args.bptDinuc)
]

end2_values = [
    filtered_table['new_bpT2'],
    filtered_table['new_bpT2'] + int(args.bptDinuc),
    filtered_table['new_bpT2'] + int(args.bptDinuc),
    filtered_table['new_bpT2']
]

# Apply the conditions using np.select
dummy_fusion['start1'] = np.select(conditions, start1_values, default=np.nan)
dummy_fusion['end1'] = np.select(conditions, end1_values, default=np.nan)
dummy_fusion['start2'] = np.select(conditions, start2_values, default=np.nan)
dummy_fusion['end2'] = np.select(conditions, end2_values, default=np.nan)


from Bio import SeqIO
from Bio.Seq import Seq

def extract_and_reverse_sequences(dummy_fusion, sequence_dict):
    """
    This function extracts left and right sequences based on the start and end positions, 
    reverse complements them if the strand is negative, and stores them in new columns.
    
    Parameters:
    - dummy_fusion (pandas DataFrame): The DataFrame containing the data.
    - sequence_dict (dict): A dictionary containing sequences with record ids as keys.

    Returns:
    - dummy_fusion (pandas DataFrame): The updated DataFrame with left and right sequences.
    """
    
    # Iterate over each row in the DataFrame
    for idx, row in dummy_fusion.iterrows():
        # Extract the left and right sequences based on the start and end positions
        left_seq = sequence_dict.get(row['chr_x'])[int(row['start1']):int(row['end1'])] if row['chr_x'] in sequence_dict else ""
        right_seq = sequence_dict.get(row['chr_y'])[int(row['start2']):int(row['end2'])] if row['chr_y'] in sequence_dict else ""

        # Reverse complement if strand is negative
        if row['strand_x'] == '-':
            left_seq = str(Seq(left_seq).reverse_complement())
        if row['strand_y'] == '-':
            right_seq = str(Seq(right_seq).reverse_complement())

        # Store the sequences in new columns in the DataFrame
        dummy_fusion.at[idx, 'left_seq'] = left_seq
        dummy_fusion.at[idx, 'right_seq'] = right_seq
    
    return dummy_fusion

# Example usage:
# Parse the DNA sequences from the fasta file
sequence_dict = {record.id: str(record.seq) for record in SeqIO.parse(args.DNAseq, "fasta")}

# Call the function
dummy_fusion = extract_and_reverse_sequences(dummy_fusion, sequence_dict)

# Construct SplicePattern
dummy_fusion['SplicePattern'] = dummy_fusion['left_seq'] + '_' + dummy_fusion['right_seq']
dummy_fusion.drop(['end2', 'start2', 'end1', 'start1', 'left_seq', 'right_seq'], axis=1, inplace=True)


##  JUNCTION SEQUENCE 
# Define conditions for each combination of strand_x and strand_y
conditions = [
    (filtered_table['strand_x'] == '+') & (filtered_table['strand_y'] == '+'),
    (filtered_table['strand_x'] == '-') & (filtered_table['strand_y'] == '-'),
    (filtered_table['strand_x'] == '+') & (filtered_table['strand_y'] == '-'),
    (filtered_table['strand_x'] == '-') & (filtered_table['strand_y'] == '+')
]

# Define corresponding start1, end1, start2, end2 values for each condition
start1_values = [
    (filtered_table['new_bpT1'] - int(args.bptSeq))+1,
    filtered_table['new_bpT1']-1,
    (filtered_table['new_bpT1'] - int(args.bptSeq))+1,
    filtered_table['new_bpT1']-1
]

end1_values = [
    filtered_table['new_bpT1']+1,
    (filtered_table['new_bpT1'] + int(args.bptSeq))-1,
    filtered_table['new_bpT1']+1,
    (filtered_table['new_bpT1'] + int(args.bptSeq))-1
]

start2_values = [
    filtered_table['new_bpT2']-1,
    (filtered_table['new_bpT2'] - int(args.bptSeq))+1,
    (filtered_table['new_bpT2'] - int(args.bptSeq))+1,
    filtered_table['new_bpT2']-1
]

end2_values = [
    (filtered_table['new_bpT2'] + int(args.bptSeq))-1,
    filtered_table['new_bpT2']+1,
    filtered_table['new_bpT2']+1,
    (filtered_table['new_bpT2'] + int(args.bptSeq))-1
]

# Apply the conditions using np.select
dummy_fusion['start1'] = np.select(conditions, start1_values, default=np.nan)
dummy_fusion['end1'] = np.select(conditions, end1_values, default=np.nan)
dummy_fusion['start2'] = np.select(conditions, start2_values, default=np.nan)
dummy_fusion['end2'] = np.select(conditions, end2_values, default=np.nan)

# Function to calculate homology between two sequences
def calculate_homology(seq1, seq2):
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of the same length.")
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
    return (matches / len(seq1)) * 100

# Parse the DNA sequences from the fasta file
sequence_dict = {record.id: str(record.seq) for record in SeqIO.parse(args.DNAseq, "fasta")}

# Call the function
dummy_fusion = extract_and_reverse_sequences(dummy_fusion, sequence_dict)

# Calculate homology and store it in a new column in the DataFrame
for idx, row in dummy_fusion.iterrows():
    left_seq = row['left_seq']
    right_seq = row['right_seq']
    
    # Calculate homology between the two sequences if they are of the same length
    if len(left_seq) == len(right_seq) and left_seq and right_seq:
        homology = calculate_homology(left_seq, right_seq)
    else:
        homology = 0.0  # If sequences are not of the same length or empty, set homology to 0
    
    # Store the calculated homology in the DataFrame
    dummy_fusion.at[idx, '%Homology'] = homology


dummy_fusion['FusionJunctionSequence'] = dummy_fusion['left_seq'] + '_' + dummy_fusion['right_seq']


#print(dummy_fusion)


# I think this should come at the end
dummy_fusion['chr_x'] = dummy_fusion['chr_x'].apply(lambda x: 'chr' + str(x))
dummy_fusion['chr_y'] = dummy_fusion['chr_y'].apply(lambda x: 'chr' + str(x))
dummy_fusion['Fusion_Name'] = dummy_fusion[['gene_id_x','gene_id_y']].apply(lambda x: '::'.join(x), axis=1)
dummy_fusion["5'Breakpoint"] = dummy_fusion[['chr_x', 'new_bpT1', 'strand_x']].apply(lambda x: ':'.join(map(str, x)), axis=1)
dummy_fusion["3'Breakpoint"] = dummy_fusion[['chr_y', 'new_bpT2', 'strand_y']].apply(lambda x: ':'.join(map(str, x)), axis=1)
dummy_fusion["5'Gene ID"] = dummy_fusion["gene_id_x"]
dummy_fusion["3'Gene ID"] = dummy_fusion["gene_id_y"]
dummy_fusion["5'Transcript_ID"]=dummy_fusion['Transcript_ID_T1']
dummy_fusion["3'Transcript_ID"]=dummy_fusion['Transcript_ID_T2']
dummy_fusion['Sites'] = np.where(dummy_fusion['chr_x'] == dummy_fusion['chr_y'], 'INTRACHROMOSOMAL', 'INTERCHROMOSOMAL')


    # Apply filtering based on the given thresholds
dummy_fusion = dummy_fusion.loc[~((dummy_fusion['Split_support_fbpt'] <= int(args.split)) & 
                                      (dummy_fusion['span_read_count'] <= int(args.span)))]
    
    # Create a new column 'SRC' by summing 'Split_support_fbpt' and 'span_read_count'
dummy_fusion['SRC'] = dummy_fusion['Split_support_fbpt'] + dummy_fusion['span_read_count']
    
    # Drop the 'Split_support_fbpt' and 'span_read_count' columns
#dummy_fusion.drop(['Split_support_fbpt', 'span_read_count'], axis=1, inplace=True)
    
    # Create a new column 'SpliceSiteClass' based on 'SplicePattern'
dummy_fusion['SpliceSiteClass'] = np.where(dummy_fusion['SplicePattern'] == 'GT-AG', 'CanonicalPattern', 'NonCanonicalPattern')
    
#dummy_fusion.to_csv("final_out.csv",sep=',')

dummy_fusion2 = dummy_fusion[['Fusion_Name',"5'Gene ID","5'Breakpoint","5'Transcript_ID","3'Gene ID","3'Transcript_ID","3'Breakpoint",'LeftBreakEntropy','RightBreakEntropy', 'FusionJunctionSequence', 'SplicePattern','Split_support_fbpt','span_read_count', 'span_read_count_uniq_queryid_count','SRC','Sites', 't1_region' ,'t2_region' ]]


##filtering::::::::::::::::::::::::::::::::::::::::::::::
# Assuming df is your DataFrame
# Group by '5'Gene ID' and '3'Gene ID'

df_top3 = (
    dummy_fusion2.groupby(["5'Gene ID", "3'Gene ID"], group_keys=False)  # Ensure group_keys=False to avoid adding a multi-index
    .apply(lambda x: x.sort_values(by='SRC', ascending=False).head(int(args.count1)))  # Sort by 'SRC' descending and get top N rows
    .reset_index(drop=True)  # Reset the index
)


# Count unique '3'Gene ID' for each '5'Gene ID'
count_df = df_top3.groupby("5'Gene ID")["3'Gene ID"].nunique().reset_index()


# Filter out the rows where the count of unique strings in col2 is more than 3
to_discard = count_df[count_df["3'Gene ID"] > int(args.promis)]["5'Gene ID"]
#to_discard
# Remove the rows from the original DataFrame where col1 has more than 3 unique col2 values
df_filtered = df_top3[~df_top3["5'Gene ID"].isin(to_discard)]

df_filtered.drop_duplicates(inplace = True)

#df_filtered.to_csv('s272_final_chimera_out.csv', sep ='\t', index=False)
#:::::::::::::::::::::::::::::::::::::::::::
#further filtering to get the unique chimera based on entropy and SRC count 


# Calculate total entropy
df_filtered['total_entropy'] = df_filtered['LeftBreakEntropy'] + df_filtered['RightBreakEntropy']

# Group by specified columns
#group_columns = ['Fusion_Name', "5'Gene ID", "5'Transcript_ID", "3'Gene ID", "3'Transcript_ID"]
group_columns = ['Fusion_Name']
# Custom aggregation logic
def custom_filter(group):
    # Keep the rows with maximum SRC
    max_src_rows = group[group['SRC'] == group['SRC'].max()]
    # If tie, keep the row with maximum total entropy
    if len(max_src_rows) > 1:
        return max_src_rows.loc[max_src_rows['total_entropy'].idxmax()]
    return max_src_rows.iloc[0]

# Apply the custom filter to each group
final_df = df_filtered.groupby(group_columns, group_keys=False).apply(custom_filter).reset_index(drop=True)

# Drop the total_entropy column if not needed
final_df.drop(columns=['total_entropy'], inplace=True)


#:::::::::::::::::::::::::::::::::::::::::::


if final_df.empty:
    print("No chimera detected for this sample.\n")
else:
# df_filtered.to_csv('FinalFuZen254.csv', sep ='\t', index=False)
    final_df.to_csv(args.output, sep ='\t', index=False)
## OLD

output_dir = Path(output_dir)

# Iterate over all .fasta and .bed files and delete them
for file in output_dir.glob('*.fasta'):
    try:
        file.unlink()  # Delete the file
        #print(f"Deleted: {file}")
    except Exception as e:
        print(f"Error deleting {file}: {e}")


for file in output_dir.glob('*.bed'):
    try:
        file.unlink()  # Delete the file
        #print(f"Deleted: {file}")
    except Exception as e:
        print(f"Error deleting {file}: {e}")

end_time = time.time()


print(f"Execution Time: {end_time - start_time} seconds")























