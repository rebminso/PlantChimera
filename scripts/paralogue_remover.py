# Import required libraries
import pandas as pd
import argparse
import os

# Create the argument parser
parser = argparse.ArgumentParser(description='Process SAfile.txt to remove duplication and paralogue gene pairs')

# Add required arguments
parser.add_argument('inputSA', help='Path to the previous output file (input1)')
parser.add_argument('gtf', help='Path to the GTF file (input2)')
parser.add_argument('transcript_info', help='Path to the output transcript_info.txt file from GTF')
parser.add_argument('outputSA', help='Path to the output new_SAfile')
parser.add_argument('paralogue_file', help='Path to the paralogue_gene.txt file (input3)', nargs='?', default='')

# Parse the arguments
args = parser.parse_args()

# Read the input file
fusion_df = pd.read_csv(args.inputSA, delimiter='\t')

# Read the GTF file
gtf_df = pd.read_csv(args.gtf, sep='\t', header=None, comment='#', low_memory=False)

def id_version_column(base_id, version):
    if version == 'NA':
        return base_id
    else:
        return f"{base_id}.{version}"
    
def transform_gtf_to_df(input_df):
    # Define the keys that will be used as headers in the CSV
    keys = ["gene_id", "gene_version", "transcript_id", "transcript_version", "gene_name", "gene_source",
            "gene_biotype", "transcript_name", "transcript_source",
            "transcript_biotype", "exon_number", "exon_version", "exon_id", "tag"]

    # Initialize an empty list to store the rows of the output DataFrame
    output_rows = []

    for _, row in input_df.iterrows():
        # Split the 9th column by semicolons to get the key-value pairs
        kv_pairs = row[8].split('; ')

        # Create a dictionary to store the key-value pairs, initialize with 'NA'
        kv_dict = {key: 'NA' for key in keys}

        for pair in kv_pairs:
            # Split each pair by space to separate key and value
            key_value = pair.split(' ', 1)
            # Only process if we have both key and value (some pairs may not be well-formed)
            if len(key_value) == 2:
                key, value = key_value
                value = value.strip('"').strip('";')
                # Add to dictionary if key is in the predefined list
                if key in keys:
                    kv_dict[key] = value

        # Combine the values from columns 0 to 7 with the kv_dict
        combined_row = {f"col_{i}": row[i] for i in range(8)}
        combined_row.update(kv_dict)
        
        # Create new columns using the create_combined_column function
        combined_row['gene_id_new'] = id_version_column(combined_row.get('gene_id', 'NA'), combined_row.get('gene_version', 'NA'))
        combined_row['transcript_id_new'] = id_version_column(combined_row.get('transcript_id', 'NA'), combined_row.get('transcript_version', 'NA'))
        combined_row['exon_id_new'] = id_version_column(combined_row.get('exon_id', 'NA'), combined_row.get('exon_version', 'NA'))
        
        # Append the dictionary values as a new row to the output list
        output_rows.append(combined_row)

    # Create a new DataFrame from the output rows
    output_df = pd.DataFrame(output_rows)

    # Drop unnecessary columns
    columns_to_drop = ['gene_id', 'gene_version', 'transcript_id', 'transcript_version', 'exon_id', 'exon_version', 'exon_number']
    output_df = output_df.drop(columns=columns_to_drop, errors='ignore')
    new_column_names = {
        "col_0": "chr",
        "col_1": "source",
        "col_2": "feature",
        "col_3": "start",
        "col_4": "end",
        "col_5": "score",
        "col_6": "strand",
        "col_7": "frame",
        "transcript_id_new": "transcript_id",
        "gene_id_new": "gene_id"
    }
    output_df = output_df.rename(columns=new_column_names)

    return output_df

# Filter rows where the 'feature' column value is "transcript"
transcript_df = gtf_df[gtf_df.iloc[:, 2] == "transcript"]
print("Extracting transcripts information from the GTF file for the next step...")

# Transform the DataFrame
new_transcript_df = transform_gtf_to_df(transcript_df)

# Write the transcript information to the output file
new_transcript_df.to_csv(args.transcript_info, sep=',')

# Merge two dataframes with transcript_id as primary key
merged_df = pd.merge(fusion_df, new_transcript_df[['chr', 'feature', 'start', 'end', 'strand', 'gene_id', 'transcript_id', 'gene_name']], left_on='transcript1_id', right_on='transcript_id', how='inner')
merged_df2 = pd.merge(merged_df, new_transcript_df[['chr', 'feature', 'start', 'end', 'strand', 'gene_id', 'transcript_id', 'gene_name']], left_on='transcript2_id', right_on='transcript_id', how='inner')
merged_df3 = merged_df2[['read_id', 'transcript1_id', 'transcript2_id', 'gene_id_x', 'transcript_id_x', 'gene_name_x', 'gene_id_y', 'transcript_id_y', 'gene_name_y']]
merged_df4 = merged_df3.loc[~(merged_df3['gene_id_x'] == merged_df3['gene_id_y'])]

# Remove the paralogue gene fusion if paralogue_file is provided
if args.paralogue_file:
    if os.path.exists(args.paralogue_file):
        # Read the paralogue file
        para = pd.read_csv(args.paralogue_file, sep='\t')

        # Create tuples of (gene_id, para_gene_id) from para DataFrame
        para_tuples = set(para.apply(lambda row: (row['GeneID'], row['ParalogGeneID']), axis=1))

        # Create tuples of (gene1, gene2) and (gene2, gene1) from merged_df4 DataFrame
        merged_df4_tuples = merged_df4.apply(lambda row: {(row['gene_id_x'], row['gene_id_y']), (row['gene_id_y'], row['gene_id_x'])}, axis=1)

        # Function to check if any tuple in merged_df4_tuples is in para_tuples
        def is_in_para_tuples(tuples_set):
            return any(t in para_tuples for t in tuples_set)

        # Filter rows where any (gene1, gene2) or (gene2, gene1) is in para_tuples
        filtered_merged_df5 = merged_df4[~merged_df4_tuples.apply(is_in_para_tuples)]

        # Select the necessary columns and remove duplicates
        filtered_merged_df6 = filtered_merged_df5[['read_id', 'transcript1_id', 'transcript2_id']]
        final_df6 = filtered_merged_df6.drop_duplicates()
    else:
        print(f"Warning: Paralogue file {args.paralogue_file} does not exist. No filtering will be done.")
        final_df6 = merged_df4[['read_id', 'transcript1_id', 'transcript2_id']].drop_duplicates()
else:
    # If the paralogue file is not provided, use the original merged_df4 without filtering
    final_df6 = merged_df4[['read_id', 'transcript1_id', 'transcript2_id']].drop_duplicates()

# Write the DataFrame (filtered or original) to the output file
final_df6.to_csv(args.outputSA, sep='\t', index=False)
