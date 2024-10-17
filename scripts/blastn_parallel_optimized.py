#############################################################################################################################
#                                                    BLAST parallel (Optimized)
#############################################################################################################################

import pandas as pd
import argparse
import subprocess
import tempfile
import os
import glob
from Bio import SeqIO
from Bio.Blast import NCBIXML
from tqdm import tqdm  # import tqdm for progress bar
import resource
import time

start_time = time.time()

# Function to generate a single BLASTn command
def generate_blastn_command(query_file, subject_file, output_file, num_threads):
    blastn_command = [
        'blastn',
        '-query', query_file,
        '-subject', subject_file,
        '-out', output_file,
        '-outfmt', '5',
        '-strand', 'both',
        '-word_size', args.wordsize,
        '-perc_identity', args.perIdentity,
        #'-num_threads', str(num_threads) ## '-num_threads' is removed because it's ignored with -subjec
    ]
    return " ".join(blastn_command)  # Join the list into a single string

# Parse arguments
parser = argparse.ArgumentParser(description="Script to run BLASTn")
parser.add_argument('--input', required=True, help='Comma-separated list of input FASTQ files')
parser.add_argument('--transcript', required=True, help='Transcript FASTA file')
parser.add_argument('--safile', required=True, help='SA file in CSV format')
parser.add_argument('--temp_path', required=True, help='Temporary path for output files')
parser.add_argument('--num_threads', type=int, default=1, help='Number of threads to use for BLASTn')
parser.add_argument('--output', required=True, help='Output file for BLAST results')
parser.add_argument('--wordsize')
parser.add_argument('--perIdentity')
parser.add_argument('--max_parallel', type=int, default=4, help='Maximum number of parallel jobs to run')
args = parser.parse_args()

# Function to process sequences one by one using a generator
def filter_sequences(input_files, transcript_file_path, id_list_file, temp_path, file_suffix=".fa"):
    srr_id_list = set()
    transcript_id_list = set()
    with open(id_list_file, 'r') as f:
        for line in f:
            columns = line.strip().split('\t')
            if len(columns) >= 3:
                srr_id_list.add(columns[0])
                transcript_id_list.add(columns[1])
                transcript_id_list.add(columns[2])
    # Log the IDs for debugging purposes
    #print(f"Extracted IDs: {srr_id_list}")

    os.makedirs(temp_path, exist_ok=True)

    # Process each input file
    input_file_list = input_files.split(",")
    for ind, input_file in enumerate(input_file_list):
        if not (input_file.endswith('.fastq') or input_file.endswith('.fq')):
            print(f"Error: {input_file} is not a valid FASTQ file. Please provide files with .fastq or .fq extension.")
            continue  # Skip this file if it's not valid

        with open(input_file, 'r') as in_file:
        # Set the file format to fastq
            file_format = 'fastq'
            # Process records from the input file
            for record in tqdm(SeqIO.parse(in_file, file_format), desc=f"Processing {os.path.basename(input_file)}"):
                # Check if the record ID is in the srr_id_list
                if record.id in srr_id_list:
                    # Construct the path for the output file
                    temp_file_path = os.path.join(temp_path, f"{record.id}_{ind + 1}{file_suffix}")

                    # Only write if the file doesn't already exist
                    if not os.path.exists(temp_file_path):
                        SeqIO.write(record, temp_file_path, "fasta")
                        #print(f"Wrote record {record.id} to {temp_file_path}")  # Log successful write
                    else:
                        print(f"File {temp_file_path} already exists. Skipping...")


    # Process the transcript file lazily (using generator)
    with open(transcript_file_path, 'r') as transcript_file:
        for x in tqdm(SeqIO.parse(transcript_file, 'fasta'), desc="Processing Transcripts"):
            if x.id in transcript_id_list:
                temp_file_path = os.path.join(temp_path, f"{x.id}.fa")
                if not os.path.exists(temp_file_path):
                    SeqIO.write(x, temp_file_path, "fasta")

print("Extracting matching transcripts and sample read sequences...")
filter_sequences(args.input, args.transcript, args.safile, args.temp_path)

# Read and process SA file in chunks
print("Preparing BLASTn commands...")

# List to store BLASTn commands (processed in batches)
blastn_commands = []

# Process SA file in chunks to reduce memory consumption
chunk_size = 100  # Process in batches of 100 rows
for chunk in pd.read_csv(args.safile, sep="\t", chunksize=chunk_size):
    chunk['read_id'] = chunk.groupby('read_id').cumcount().astype(str) + '#' + chunk['read_id']
    for index, row in chunk.iterrows():
        for x in range(len(args.input.split(","))):
            # Generate BLASTn command for transcript1
            query_file = f"{args.temp_path}/{row['read_id'][2:]}_{str(x+1)}.fa"
            subject_file1 = f"{args.temp_path}/{row['transcript1_id']}.fa"
            output_file1 = f"{args.temp_path}/{row['read_id']}_{str(x+1)}_{row['transcript1_id']}.xml"
            blastn_commands.append(generate_blastn_command(query_file, subject_file1, output_file1, args.num_threads))

            # Generate BLASTn command for transcript2
            subject_file2 = f"{args.temp_path}/{row['transcript2_id']}.fa"
            output_file2 = f"{args.temp_path}/{row['read_id']}_{str(x+1)}_{row['transcript2_id']}.xml"
            blastn_commands.append(generate_blastn_command(query_file, subject_file2, output_file2, args.num_threads))

# Save BLASTn commands to a temporary file
commands_file = os.path.join(args.temp_path, 'blastn_commands.txt')
with open(commands_file, 'w') as f:
    for command in blastn_commands:
        f.write(command + '\n')

# Run the commands in parallel using GNU Parallel
print(f"Running BLASTn using GNU Parallel with {args.max_parallel} parallel jobs...")
subprocess.run(['parallel', '-j', str(args.max_parallel), '-a', commands_file])

print("BLASTn runs completed.")

# Parse the output BLASTn XML files and save to CSV
df_blast_rows = []  # Store results as list of rows
for file in glob.glob(args.temp_path + "/*.xml"):
    filename = file.split("/")[-1]  # Get the filename from the path
    #q_ID = file.split("/")[-1][:file.split("/")[-1].rfind("_")]
    q_ID = (filename.split('_1_')[0] + '_1') if '_1_' in filename else (filename.split('_2_')[0] + '_2') if '_2_' in filename else filename
    result_handle = open(file)
    blast_records = NCBIXML.parse(result_handle)
    for blast_record in blast_records:
        query_ID = blast_record.query.split(' ')[0]
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                df_blast_rows.append({
                    'query_ID': q_ID,
                    'Transcript_ID': alignment.title.split(' ')[1],
                    'query_start': hsp.query_start,
                    'query_end': hsp.query_end,
                    'subject_start': hsp.sbjct_start,
                    'subject_end': hsp.sbjct_end,
                    'strand': "/".join(hsp.strand),
                    'expect': hsp.expect,
                    'identity': hsp.identities * 100 / hsp.align_length
                })

# Convert to DataFrame and save
df_blast = pd.DataFrame(df_blast_rows)
df_blast.to_csv(args.output, sep='\t', index=False)

#clean up memory

end_time = time.time()
usage = resource.getrusage(resource.RUSAGE_SELF)
max_memory = usage.ru_maxrss  # Memory in kilobytes

print(f"Execution Time: {end_time - start_time} seconds")
print(f"Max Memory Usage: {max_memory / 1024} MB")

