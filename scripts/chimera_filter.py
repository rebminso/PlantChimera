#############################################################################################################################
#                                     Candidate Chimera Filtering
#############################################################################################################################

import re
import os
import math
import glob
import argparse
import subprocess
import pandas as pd
import numpy as np
from Bio import SeqIO
#from Bio.Blast.Applications import NcbiblastnCommandline
#from Bio.Blast import NCBIXML
from collections import Counter
import warnings
warnings.filterwarnings("ignore", category=pd.errors.SettingWithCopyWarning)
import resource
import time

start_time = time.time()

parser = argparse.ArgumentParser(description='a list of seq from csv file')
parser.add_argument('input',help='fusionCandidates')
parser.add_argument('JDis', help='JunctionDistance')
parser.add_argument('Sent', help='ShannonEntropy')
parser.add_argument('DNAseq',help='dnatopLevel')
parser.add_argument('bptSeq', help='Sequence around the bpt')
parser.add_argument('bptDinuc')
parser.add_argument('count1', help='top rows for the same transcript pair')
parser.add_argument('promis', help='PROMISCIOUS')
parser.add_argument('split', help='splitCount')
parser.add_argument('span', help='spanCount')
parser.add_argument('output', help='finalFusion')
args = parser.parse_args()

output_dir = os.path.dirname(args.output)

filtered_table = pd.read_csv(args.input,delimiter='\t')
filtered_table.rename(columns={'0': 'Split_support_fbpt'}, inplace=True)

# JUNCTION DISTANCE
# # Filter based on the condition abs(new_bp1 - new_bp2) >= 100
filtered_table =  filtered_table[abs((filtered_table['new_bpT1']) - (filtered_table['new_bpT2'])) >= int(args.JDis)]

dummy_fusion = filtered_table

# filtered_data = Mod4_1
#### SHANNONG ENTORPY

# Define conditions for each combination of strand_x and strand_y
conditions = [
    (filtered_table['strand_x'] == '+') & (filtered_table['strand_y'] == '+'),
    (filtered_table['strand_x'] == '-') & (filtered_table['strand_y'] == '-'),
    (filtered_table['strand_x'] == '+') & (filtered_table['strand_y'] == '-'),
    (filtered_table['strand_x'] == '-') & (filtered_table['strand_y'] == '+')
]

# Define corresponding start1, end1, start2, end2 values for each condition
start1_values = [
    filtered_table['new_bpT1'] - int(args.Sent),
    filtered_table['new_bpT1'],
    filtered_table['new_bpT1'] - int(args.Sent),
    filtered_table['new_bpT1']
]

end1_values = [
    filtered_table['new_bpT1'],
    filtered_table['new_bpT1'] + int(args.Sent),
    filtered_table['new_bpT1'],
    filtered_table['new_bpT1'] + int(args.Sent)
]

start2_values = [
    filtered_table['new_bpT2'],
    filtered_table['new_bpT2'] - int(args.Sent),
    filtered_table['new_bpT2'] - int(args.Sent),
    filtered_table['new_bpT2']
]

end2_values = [
    filtered_table['new_bpT2'] + int(args.Sent),
    filtered_table['new_bpT2'],
    filtered_table['new_bpT2'],
    filtered_table['new_bpT2'] + int(args.Sent)
]

# Apply the conditions using np.select
dummy_fusion['start1'] = np.select(conditions, start1_values, default=np.nan)
dummy_fusion['end1'] = np.select(conditions, end1_values, default=np.nan)
dummy_fusion['start2'] = np.select(conditions, start2_values, default=np.nan)
dummy_fusion['end2'] = np.select(conditions, end2_values, default=np.nan)

# Select the required columns
left_fusion_cord = dummy_fusion[['chr_x', 'Transcript_ID_T1', 'start1', 'end1', 'gene_name_x', 'strand_x', 't1_region']]

# Function to convert DataFrame to BED file
def df1_to_bed(df, bed_file_path):
    with open(bed_file_path, 'w') as f:
        for index, row in df.iterrows():
            f.write(f"{row['chr_x']}\t{int(row['start1'])}\t{int(row['end1'])}\t{row['gene_name_x']}\t0\t{row['strand_x']}\t0\t{row['t1_region']}\n")

# Convert DataFrame to BED file
bed_file = 'left_fusion_cord.bed'
left_bed_file_path = os.path.join(output_dir, bed_file)

df1_to_bed(left_fusion_cord, left_bed_file_path)

# left_fusion_cord.to_csv("left_fusion_cord.bed", sep='\t', header=False, index=False)
right_fusion_cord =  dummy_fusion[['chr_y', 'Transcript_ID_T2', 'start2','end2','gene_name_y', 'strand_y', 't2_region']]

# Function to convert DataFrame to BED file
def df2_to_bed(df, bed_file_path):
    with open(bed_file_path, 'w') as f:
        for index, row in df.iterrows():
            f.write(f"{row['chr_y']}\t{int(row['start2'])}\t{int(row['end2'])}\t{row['gene_name_y']}\t0\t{row['strand_y']}\t0\t{row['t2_region']}\n")

# Convert DataFrame to BED file
bed_file = 'right_fusion_cord.bed'
right_bed_file_path = os.path.join(output_dir, bed_file)
df2_to_bed(right_fusion_cord, right_bed_file_path)

def run_bedtools_getfasta(fasta_file, bed_file, output_file):
    # Construct the bedtools command
    command = ['bedtools', 'getfasta', '-fi', fasta_file, '-bed', bed_file, '-name', '-s', '-fo', output_file]

    # Execute the command
    try:
        subprocess.run(command, check=True)
        #print(f"Output saved to {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running bedtools: {e}")

# Example usage
fasta_file = args.DNAseq
#left_bed_file = "left_fusion_cord.bed"
left_output_file = "left_fusion.fasta"
left_output_file_path = os.path.join(output_dir,left_output_file)
#right_bed_file = "right_fusion_cord.bed"
right_output_file = "right_fusion.fasta"
right_output_file_path = os.path.join(output_dir,right_output_file)

# Run bedtools getfasta for left fusion cord
run_bedtools_getfasta(fasta_file, left_bed_file_path, left_output_file_path)

# Run bedtools getfasta for right fusion cord
run_bedtools_getfasta(fasta_file, right_bed_file_path, right_output_file_path)

def calculate_entropy(sequence):
    # Calculate the frequency of each nucleotide
    count = Counter(sequence)
    total_bases = len(sequence)

    # Calculate the probability of each nucleotide
    probabilities = [count[base] / total_bases for base in 'ATCG']

    # Calculate entropy
    entropy = -sum(p * math.log2(p) for p in probabilities if p > 0)

    return entropy

def read_sequences(file_path):
    sequences = []
    for record in SeqIO.parse(file_path, "fasta"):
        sequences.append(str(record.seq))
    return sequences

# # Paths to the files containing the sequences
# left_file_path = "left_sequences.fasta"
# right_file_path = "right_sequences.fasta"
# Read sequences from files
left_sequences = read_sequences(left_output_file_path)
right_sequences = read_sequences(right_output_file_path)

# Calculate entropy for each sequence
left_entropies = [calculate_entropy(seq) for seq in left_sequences]
right_entropies = [calculate_entropy(seq) for seq in right_sequences]

# Add the calculated entropy values to the original DataFrame
dummy_fusion['LeftBreakEntropy'] = left_entropies
dummy_fusion['RightBreakEntropy'] = right_entropies

#dummy_fusion.shape
dum_fusion = dummy_fusion.drop(['end2','start2','end1','start1'],axis=1)
dum_fusion['LeftBreakEntropy']= dum_fusion['LeftBreakEntropy'].round(2)
dum_fusion['RightBreakEntropy']= dum_fusion['RightBreakEntropy'].round(2)



# dum_fusion


# I think this should come at the end
dum_fusion['chr_x'] = dum_fusion['chr_x'].apply(lambda x: 'chr' + str(x))
dum_fusion['chr_y'] = dum_fusion['chr_y'].apply(lambda x: 'chr' + str(x))
dum_fusion['Fusion_Name'] = dum_fusion[['Transcript_ID_T1','Transcript_ID_T2']].apply(lambda x: '_'.join(x), axis=1)
# dum_fusion['LeftBpt'] = dum_fusion[[ 'chr_x', 'new_bp1','final_transcript_strand_1']].apply(lambda x: ':'.join(x), axis=1)
dum_fusion["5'Breakpoint"] = dum_fusion[['chr_x', 'new_bpT1', 'strand_x']].apply(lambda x: ':'.join(map(str, x)), axis=1)
dum_fusion["3'Breakpoint"] = dum_fusion[['chr_y', 'new_bpT2', 'strand_y']].apply(lambda x: ':'.join(map(str, x)), axis=1)
# dum_fusion['LeftGene'] = dum_fusion[['gene_name_x','Transcript_ID_T1']].apply(lambda x: '_'.join(x), axis=1)
# dum_fusion['RightGene'] = dum_fusion[['gene_name_y','Transcript_ID_T2']].apply(lambda x: '_'.join(x), axis=1)
dum_fusion["5'Gene ID"] = dum_fusion["Transcript_ID_T1"]
dum_fusion["3'Gene ID"] = dum_fusion["Transcript_ID_T2"]
dum_fusion['Sites'] = np.where(dum_fusion['chr_x'] == dum_fusion['chr_y'], 'INTRACHROMOSOMAL', 'INTERCHROMOSOMAL')
dum_fusion2 = dum_fusion[['Fusion_Name',"5'Gene ID","5'Breakpoint","3'Gene ID","3'Breakpoint",'LeftBreakEntropy','RightBreakEntropy', 'Split_support_fbpt','span_read_count', 'Sites', 't1_region' ,'t2_region' ]]
dum_fusion2.to_csv('dum_fusion2.csv',sep='\t',index=False)
dum_fusion2



# !apt-get install bedtools
#    totransfer


dummy_fusion = filtered_table
# dummy_fusion['start1'], dummy_fusion['end1'], dummy_fusion['start2'], dummy_fusion['end2'] = (dummy_fusion['new_bpT1'] - int(args.bptSeq) , dummy_fusion['new_bpT1'], dummy_fusion['new_bpT2'] , dummy_fusion['new_bpT2'] + int(args.bptSeq))
# left_fusion_cord=  dummy_fusion[['chr_x', 'Transcript_ID_T1', 'start1','end1','gene_name_x', 'strand_x', 't1_region']]

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

# Select the required columns
left_fusion_cord = dummy_fusion[['chr_x', 'Transcript_ID_T1', 'start1', 'end1', 'gene_name_x', 'strand_x', 't1_region']]


# Function to convert DataFrame to BED file
def df1_to_bed(df, bed_file):
    with open(bed_file, 'w') as f:
        for index, row in df.iterrows():
            f.write(f"{row['chr_x']}\t{int(row['start1'])}\t{int(row['end1'])}\t{row['gene_name_x']}\t0\t{row['strand_x']}\t0\t{row['t1_region']}\n")

# Convert DataFrame to BED file
bed_file = 'left_fusion_cord.bed'
left_bed_file_path = os.path.join(output_dir, bed_file)
df1_to_bed(left_fusion_cord, left_bed_file_path)

# left_fusion_cord.to_csv("left_fusion_cord.bed", sep='\t', header=False, index=False)
right_fusion_cord =  dummy_fusion[['chr_y', 'Transcript_ID_T2', 'start2','end2','gene_name_y', 'strand_y', 't2_region']]

# Function to convert DataFrame to BED file
def df2_to_bed(df, bed_file):
    with open(bed_file, 'w') as f:
        for index, row in df.iterrows():
            f.write(f"{row['chr_y']}\t{int(row['start2'])}\t{int(row['end2'])}\t{row['gene_name_y']}\t0\t{row['strand_y']}\t0\t{row['t2_region']}\n")

# Convert DataFrame to BED file
bed_file = 'right_fusion_cord.bed'
right_bed_file_path = os.path.join(output_dir, bed_file)
df2_to_bed(right_fusion_cord, right_bed_file_path)

def run_bedtools_getfasta(fasta_file, bed_file, output_file):
    # Construct the bedtools command
    command = ['bedtools', 'getfasta', '-fi', fasta_file, '-bed', bed_file, '-name', '-s', '-fo', output_file]

    # Execute the command
    try:
        subprocess.run(command, check=True)
        #print(f"Output saved to {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running bedtools: {e}")

# Example usage
fasta_file = args.DNAseq
#left_bed_file = "left_fusion_cord.bed"
left_output_file = "left_Sequencefusion.fasta"
left_output_file_path = os.path.join(output_dir,left_output_file)

#right_bed_file = "right_fusion_cord.bed"
right_output_file = "right_Sequencefusion.fasta"
right_output_file_path = os.path.join(output_dir,right_output_file)

# Run bedtools getfasta for left fusion cord
run_bedtools_getfasta(fasta_file, left_bed_file_path, left_output_file_path)
# Run bedtools getfasta for right fusion cord
run_bedtools_getfasta(fasta_file, right_bed_file_path, right_output_file_path)

AlSEQ = pd.read_csv(left_output_file_path,delimiter='\t',header=None)
AlSEQ = AlSEQ[~AlSEQ[0].str.contains(">", na=False)]
AlSEQ.rename(columns = {0: 'SequenceLeft'},inplace=True)
ArSEQ = pd.read_csv(right_output_file_path,delimiter='\t',header=None)
ArSEQ = ArSEQ[~ArSEQ[0].str.contains(">", na=False)]
ArSEQ.rename(columns = {0: 'SequenceRight'},inplace=True)
dum_fusion2.join([AlSEQ, ArSEQ])


# If they are not of the same length, handle the mismatch
if len(dum_fusion2) == len(AlSEQ) == len(ArSEQ):
    # Concatenate the DataFrames along columns (axis=1)
    result1 = pd.concat([dum_fusion2.reset_index(drop=True), 
                         AlSEQ.reset_index(drop=True), 
                         ArSEQ.reset_index(drop=True)], axis=1)
    # Replace NaN values with '-' in SequenceLeft and SequenceRight, then convert to string
    result1['SequenceLeft'] = result1['SequenceLeft'].fillna('-').astype(str)
    result1['SequenceRight'] = result1['SequenceRight'].fillna('-').astype(str)
    
    # Create the SplicePattern by joining SequenceLeft and SequenceRight with an underscore
    result1['SplicePattern'] = result1[['SequenceLeft', 'SequenceRight']].apply(lambda x: '_'.join(x), axis=1)
    
    # Drop the unnecessary columns if needed
    result1.drop(['SequenceLeft', 'SequenceRight'], axis=1, inplace=True)
else:
    print("Error: DataFrames are of different lengths.")

## extract dinucleotide seqeunce around btpt

dummy_fusion = filtered_table
# Define conditions for each combination of strand_x and strand_y
conditions = [
    (filtered_table['strand_x'] == '+') & (filtered_table['strand_y'] == '+'),
    (filtered_table['strand_x'] == '-') & (filtered_table['strand_y'] == '-'),
    (filtered_table['strand_x'] == '+') & (filtered_table['strand_y'] == '-'),
    (filtered_table['strand_x'] == '-') & (filtered_table['strand_y'] == '+')
]

# Define corresponding start1, end1, start2, end2 values for each condition
start1_values = [
    filtered_table['new_bpT1'] - int(args.bptSeq),
    filtered_table['new_bpT1'],
    filtered_table['new_bpT1'] - int(args.bptSeq),
    filtered_table['new_bpT1']
]

end1_values = [
    filtered_table['new_bpT1'],
    filtered_table['new_bpT1'] + int(args.bptSeq),
    filtered_table['new_bpT1'],
    filtered_table['new_bpT1'] + int(args.bptSeq)
]

start2_values = [
    filtered_table['new_bpT2'],
    filtered_table['new_bpT2'] - int(args.bptSeq),
    filtered_table['new_bpT2'] - int(args.bptSeq),
    filtered_table['new_bpT2']
]

end2_values = [
    filtered_table['new_bpT2'] + int(args.bptSeq),
    filtered_table['new_bpT2'],
    filtered_table['new_bpT2'],
    filtered_table['new_bpT2'] + int(args.bptSeq)
]

# Apply the conditions using np.select
dummy_fusion['start1'] = np.select(conditions, start1_values, default=np.nan)
dummy_fusion['end1'] = np.select(conditions, end1_values, default=np.nan)
dummy_fusion['start2'] = np.select(conditions, start2_values, default=np.nan)
dummy_fusion['end2'] = np.select(conditions, end2_values, default=np.nan)

left_fusion_cord=  dummy_fusion[['chr_x', 'Transcript_ID_T1', 'start1','end1','gene_name_x', 'strand_x', 't1_region']]

# Function to convert DataFrame to BED file
def df1_to_bed(df, bed_file):
    with open(bed_file, 'w') as f:
        for index, row in df.iterrows():
            f.write(f"{row['chr_x']}\t{int(row['start1'])}\t{int(row['end1'])}\t{row['gene_name_x']}\t0\t{row['strand_x']}\t0\t{row['t1_region']}\n")

# Convert DataFrame to BED file
bed_file = 'left_fusion_cord.bed'
left_bed_file_path =  os.path.join(output_dir, bed_file)
df1_to_bed(left_fusion_cord, left_bed_file_path)

# left_fusion_cord.to_csv("left_fusion_cord.bed", sep='\t', header=False, index=False)
right_fusion_cord =  dummy_fusion[['chr_y', 'Transcript_ID_T2', 'start2','end2','gene_name_y', 'strand_y', 't2_region']]

# Function to convert DataFrame to BED file
def df2_to_bed(df, bed_file):
    with open(bed_file, 'w') as f:
        for index, row in df.iterrows():
            f.write(f"{row['chr_y']}\t{int(row['start2'])}\t{int(row['end2'])}\t{row['gene_name_y']}\t0\t{row['strand_y']}\t0\t{row['t2_region']}\n")

# Convert DataFrame to BED file
bed_file = 'right_fusion_cord.bed'
right_bed_file_path =  os.path.join(output_dir, bed_file)
df2_to_bed(right_fusion_cord, right_bed_file_path)

def run_bedtools_getfasta(fasta_file, bed_file, output_file):
    # Construct the bedtools command
    command = ['bedtools', 'getfasta', '-fi', fasta_file, '-bed', bed_file, '-name', '-s', '-fo', output_file]

    # Execute the command
    try:
        subprocess.run(command, check=True)
        #print(f"Output saved to {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running bedtools: {e}")


fasta_file = args.DNAseq
#left_bed_file = "left_fusion_cord.bed"
left_output_file = "left_anchorfusion.fasta"
left_output_file_path = os.path.join(output_dir, left_output_file)

#right_bed_file = "right_fusion_cord.bed"
right_output_file = "right_anchorfusion.fasta"
right_output_file_path = os.path.join(output_dir, right_output_file)

# Run bedtools getfasta for left fusion cord
run_bedtools_getfasta(fasta_file, left_bed_file_path, left_output_file_path)

# Run bedtools getfasta for right fusion cord
run_bedtools_getfasta(fasta_file, right_bed_file_path, right_output_file_path)


Al = pd.read_csv(left_output_file_path,delimiter='\t',header=None)
Al = Al[~Al[0].str.contains(">", na=False)]
Al.rename(columns = {0: 'AnchorL'},inplace=True)
Ar = pd.read_csv(right_output_file_path,delimiter='\t',header=None)
Ar = Ar[~Ar[0].str.contains(">", na=False)]
Ar.rename(columns = {0: 'AnchorRight'},inplace=True)
# Example DataFrames with sequences
# Assuming both DataFrames Al and Ar have a column 'sequence' containing DNA sequences

# # Function to calculate homology between two sequences
def calculate_homology(seq1, seq2):
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of the same length.")
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
    return (matches / len(seq1)) * 100

# Calculate homology for each row
# Al['%Homology'] = Al['AnchorLeft'].combine(Ar['AnchorRight'], calculate_homology) round value to 1 decimal point
Al['%Homology'] = Al['AnchorL'].combine(Ar['AnchorRight'], lambda x, y: round(calculate_homology(x, y), 1))
Al['AnchorLeft'] = Al[['%Homology','AnchorL']].apply(lambda x: '&'.join(map(str, x)), axis=1)
dum_fusion2.join([Al, Ar])
result2 = pd.concat([dum_fusion2.reset_index(drop=True), Al.reset_index(drop=True), Ar.reset_index(drop=True)], axis=1)
result2["FusionJunctionSequence"] = result2[['AnchorLeft', 'AnchorRight']].apply(lambda x: '-'.join(map(str, x)), axis=1)
result2.drop(['AnchorLeft','AnchorRight','%Homology', 'AnchorL' ],axis=1,inplace=True)

result2[['%Homology', 'FusionJunctionSequence']] = result2['FusionJunctionSequence'].str.split('&', expand=True)
# result2.drop(columns=['span_read_count', 'Split_support_fbpt', 'RightBreakEntropy', 'LeftBreakEntropy'],inplace=True

#merge result1 and result2 dataframe and add columns which are uniq , concatenate df and add columns if are not present in result1
merged = pd.merge(result1, result2, how='outer')
# Filtering
# Removing fusion with one split read and 0 spanning read
## ADD PARAMETER
merged = merged.loc[~((merged['Split_support_fbpt']<=int(args.split)) & (merged['span_read_count']<=int(args.span)))]
merged['SRC'] = merged['Split_support_fbpt'] + merged['span_read_count']
merged.drop(['Split_support_fbpt','span_read_count'],axis=1,inplace=True)
merged['SpliceSiteClass'] = np.where(merged['SplicePattern'] == 'GT-AG', 'CanonicalPattern', 'NonCanonicalPattern')
# merged.to_csv('FinalFuZen272.csv', sep ='\t', index=False)
## REMOVE PROMISCIOUS GENES

# dum_fusion[["5'Breakpoint", "3'Breakpoint"]].value_counts()

# Assuming df is your DataFrame
df_top3 = (
    merged.groupby(["5'Gene ID", "3'Gene ID"])
    .apply(lambda x: x.nlargest(int(args.count1), 'SRC'))  # Get top rows by 'SRC' value
    .reset_index(drop=True)  # Reset the index
)

# Count unique '3'Gene ID' for each '5'Gene ID'
count_df = df_top3.groupby("5'Gene ID")["3'Gene ID"].nunique().reset_index()


# Filter out the rows where the count of unique strings in col2 is more than 3
to_discard = count_df[count_df["3'Gene ID"] > int(args.promis)]["5'Gene ID"]
to_discard
# Remove the rows from the original DataFrame where col1 has more than 3 unique col2 values
df_filtered = df_top3[~df_top3["5'Gene ID"].isin(to_discard)]
if df_filtered.empty:
    print("No chimera detected for this sample.\n")
else:
# df_filtered.to_csv('FinalFuZen254.csv', sep ='\t', index=False)
    df_filtered.to_csv(args.output, sep ='\t', index=False)
## OLD


end_time = time.time()
usage = resource.getrusage(resource.RUSAGE_SELF)
max_memory = usage.ru_maxrss  # Memory in kilobytes

#print(f"Execution Time: {end_time - start_time} seconds")
print(f"Max Memory Usage: {max_memory / 1024} MB")






















