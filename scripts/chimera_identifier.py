#############################################################################################################################
#                                     Chimeric Transcript Table Generation
# ############################################################################################################################


import re
import pandas as pd
import argparse
import numpy as np
from itertools import product
import glob
import resource
import time
import os
from joblib import Parallel, delayed
import time
from memory_profiler import memory_usage
import warnings
warnings.filterwarnings("ignore", category=pd.errors.SettingWithCopyWarning)

start_time = time.time()

parser = argparse.ArgumentParser(description='a list of seq from csv file')
parser.add_argument('gtf',help='arb gtf')
parser.add_argument('transcript_info')
parser.add_argument('exons_position_adjusted')
parser.add_argument('allsplit',help='main input file')
parser.add_argument('ftable')
parser.add_argument('allsplit_fusiondf')
args = parser.parse_args()


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
            key, value = pair.split(' ', 1)
            value = value.strip('"')
            value = value.strip('";')
            # Add to dictionary if key is in the predefined list
            if key in keys:
                kv_dict[key] = value

        # Combine the values from columns 0 to 7 with the kv_dict
        combined_row = {f"col_{i}": row[i] for i in range(8)}
        combined_row.update(kv_dict)

        # Create the new columns using the create_combined_column function
        # join gene_id, transcript_id, exon_id with their version
        combined_row['gene_id_new'] = id_version_column(combined_row.get('gene_id', 'NA'), combined_row.get('gene_version', 'NA'))
        combined_row['transcript_id_new'] = id_version_column(combined_row.get('transcript_id', 'NA'), combined_row.get('transcript_version', 'NA'))
        combined_row['exon_id_new'] = id_version_column(combined_row.get('exon_id', 'NA'), combined_row.get('exon_version', 'NA'))

        # Append the dictionary values as a new row to the output list
        output_rows.append(combined_row)

    # Create a new DataFrame from the output rows
    output_df = pd.DataFrame(output_rows)

    # Drop 'gene_id', 'gene_version', 'transcript_id', 'transcript_version', 'exon_id', and 'exon_version' columns
    columns_to_drop = ['gene_id', 'gene_version', 'transcript_id', 'transcript_version', 'exon_id', 'exon_version','exon_number']
    output_df = output_df.drop(columns=columns_to_drop, errors='ignore')
    # Rename col_0 to col_7
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


def adjust_exon_positions_in_gtfdf(df):
    # Reset index for easier row referencing
    df = df.reset_index(drop=True)

    prev_end = 0
    new_exon_start = []
    new_exon_end = []

    for index, row in df.iterrows():
        if index == 0 or row['transcript_id'] != df.iloc[index - 1]['transcript_id']:
            # For the first exon of each transcript
            start = 1
            end = (row['end'] - row['start']) + 1
        else:
            # For subsequent exons within the same transcript
            start = prev_end + 1
            end = (row['end'] - row['start']) + start

        # Update previous end for the next iteration
        prev_end = end

        # Append new values to lists
        new_exon_start.append(start)
        new_exon_end.append(end)

    # Add the new columns to the DataFrame
    df['new_exon_start'] = new_exon_start
    df['new_exon_end'] = new_exon_end

    # Efficiently swap exon start and end positions for transcripts with negative strand
    df.loc[df['strand'] == '-', ['start', 'end']] = df.loc[df['strand'] == '-', ['end', 'start']].values

    return df


def adjust_transcript_position_with_profiling(df1, df2):
    start_time = time.time()  # Start tracking time

    # Function to track memory usage
    def process_with_memory_tracking():
        return adjust_transcript_position(df1, df2)

    # Measure memory usage and get the result
    mem_usage, result = memory_usage(process_with_memory_tracking, retval=True, max_usage=True)

    end_time = time.time()  # End tracking time

    total_time = end_time - start_time  # Calculate total time taken
    max_memory = mem_usage  # Maximum memory used

    print(f"Time taken: {total_time:.2f} seconds")
    print(f"Maximum memory usage: {max_memory:.2f} MiB")

    return result

# Main function (without profiling)
def adjust_transcript_position(df1, df2):
    # Ensure the columns to be modified exist
    for col in ['new_start_T1', 'new_start_T2', 'new_end_T1', 'new_end_T2']:
        if col not in df2.columns:
            df2[col] = np.nan  # Initialize with NaN

    # Precompute conditions for matching based on transcript_id and exon ranges
    def process_single_row(index2, row2, df1):
        # Create a copy to modify
        row_modified = row2.copy()

        # Match rows for Transcript T1
        matching_rows_t1 = df1[
            (df1['transcript_id'] == row2['Transcript_ID_T1']) &
            ((row2['start_T1'] >= df1['new_exon_start']) & (row2['end_T1'] <= df1['new_exon_end']))
        ]

        # Match rows for Transcript T2
        matching_rows_t2 = df1[
            (df1['transcript_id'] == row2['Transcript_ID_T2']) &
            ((row2['start_T2'] >= df1['new_exon_start']) & (row2['end_T2'] <= df1['new_exon_end']))
        ]

        # Handle overlapping regions for T1
        matching_rows_t3 = df1[
            (df1['transcript_id'] == row2['Transcript_ID_T1']) &
            (
                ((row2['start_T1'] >= df1['new_exon_start']) & (row2['start_T1'] <= df1['new_exon_end'])) |
                ((row2['end_T1'] >= df1['new_exon_start']) & (row2['end_T1'] <= df1['new_exon_end']))
            )
        ]

        # Handle overlapping regions for T2
        matching_rows_t4 = df1[
            (df1['transcript_id'] == row2['Transcript_ID_T2']) &
            (
                ((row2['start_T2'] >= df1['new_exon_start']) & (row2['start_T2'] <= df1['new_exon_end'])) |
                ((row2['end_T2'] >= df1['new_exon_start']) & (row2['end_T2'] <= df1['new_exon_end']))
            )
        ]

        # Logic for adjusting start and end positions for T1
        for _, row1 in matching_rows_t3.iterrows():
            if row2['start_T1'] >= row1['new_exon_start'] and row2['start_T1'] <= row1['new_exon_end']:
                if row1['strand'] == '+':
                    row_modified['new_start_T1'] = int((row2['start_T1'] - row1['new_exon_start']) + row1['exon_start'])
                    if abs(row2['end_T1'] - row1['new_exon_end']) < 3:
                        row_modified['new_end_T1'] = row1['exon_end']
                elif row1['strand'] == '-':
                    row_modified['new_start_T1'] = int(row1['exon_start'] - (row2['start_T1'] - row1['new_exon_start']))
                    if abs(row2['end_T1'] - row1['new_exon_end']) < 3:
                        row_modified['new_end_T1'] = row1['exon_end']

            elif row2['end_T1'] >= row1['new_exon_start'] and row2['start_T1'] < row1['new_exon_start'] and row2['end_T1'] <= row1['new_exon_end']:
                if row1['strand'] == '+':
                    if abs(row2['end_T1'] - row1['new_exon_start']) < 3:
                        row_modified['new_end_T1'] = row1['exon_start'] + (row2['end_T1'] - row1['new_exon_start'])
                    else:
                        row_modified['new_end_T1'] = int(row1['exon_start'] + (row2['end_T1'] - row1['new_exon_start']))
                    if abs(row1['new_exon_start'] - row2['start_T1']) < 3:
                        row_modified['new_start_T1'] = row1['exon_start']
                elif row1['strand'] == '-':
                    if abs(row2['end_T1'] - row1['new_exon_start']) < 3:
                        row_modified['new_end_T1'] = row1['exon_start'] - (row2['end_T1'] - row1['new_exon_start'])
                    else:
                        row_modified['new_end_T1'] = int(row1['exon_start'] - (row2['end_T1'] - row1['new_exon_start']))
                    if abs(row1['new_exon_start'] - row2['start_T1']) < 3:
                        row_modified['new_start_T1'] = row1['exon_start']

        # Logic for adjusting start and end positions for T2
        for _, row1 in matching_rows_t4.iterrows():
            if row2['start_T2'] >= row1['new_exon_start'] and row2['start_T2'] <= row1['new_exon_end']:
                if row1['strand'] == '+':
                    row_modified['new_start_T2'] = int(row1['exon_start'] - (row2['start_T2'] - row1['new_exon_start']))
                    if abs(row2['end_T2'] - row1['new_exon_end']) < 3:
                        row_modified['new_end_T2'] = row1['exon_end']
                elif row1['strand'] == '-':
                    row_modified['new_start_T2'] = int(row1['exon_start'] - (row2['start_T2'] - row1['new_exon_start']))
                    if abs(row2['end_T2'] - row1['new_exon_end']) < 3:
                        row_modified['new_end_T2'] = row1['exon_end']

            elif row2['end_T2'] >= row1['new_exon_start'] and row2['end_T2'] <= row1['new_exon_end'] and row2['start_T2'] < row1['new_exon_start']:
                if row1['strand'] == '+':
                    if abs(row2['end_T2'] - row1['new_exon_start']) < 3:
                        row_modified['new_end_T2'] = row1['exon_start'] + (row2['end_T2'] - row1['new_exon_start'])
                    else:
                        row_modified['new_end_T2'] = int((row2['end_T2'] - row1['new_exon_start']) + row1['exon_start'])
                    if abs(row1['new_exon_start'] - row2['start_T2']) < 3:
                        row_modified['new_start_T2'] = row1['exon_start']
                elif row1['strand'] == '-':
                    if abs(row2['end_T2'] - row1['new_exon_start']) < 3:
                        row_modified['new_end_T2'] = row1['exon_start'] - (row2['end_T2'] - row1['new_exon_start'])
                    else:
                        row_modified['new_end_T2'] = int(row1['exon_start'] - (row2['end_T2'] - row1['new_exon_start']))
                    if abs(row1['new_exon_start'] - row2['start_T2']) < 3:
                        row_modified['new_start_T2'] = row1['exon_start']

        # Logic for when transcripts cover one exon (T1)
        for _, row1 in matching_rows_t1.iterrows():
            if row1['strand'] == '+':
                row_modified['new_start_T1'] = (row2['start_T1'] - row1['new_exon_start']) + row1['exon_start']
                row_modified['new_end_T1'] = row_modified['new_start_T1'] + (row2['end_T1'] - row2['start_T1'])
            elif row1['strand'] == '-':
                row_modified['new_start_T1'] = row1['exon_start'] - (row2['start_T1'] - row1['new_exon_start'])
                row_modified['new_end_T1'] = row_modified['new_start_T1'] - (row2['end_T1'] - row2['start_T1'])

        # Logic for when transcripts cover one exon (T2)
        for _, row1 in matching_rows_t2.iterrows():
            if row1['strand'] == '+':
                row_modified['new_start_T2'] = (row2['start_T2'] - row1['new_exon_start']) + row1['exon_start']
                row_modified['new_end_T2'] = row_modified['new_start_T2'] + (row2['end_T2'] - row2['start_T2'])
            elif row1['strand'] == '-':
                row_modified['new_start_T2'] = row1['exon_start'] - (row2['start_T2'] - row1['new_exon_start'])
                row_modified['new_end_T2'] = row_modified['new_start_T2'] - (row2['end_T2'] - row2['start_T2'])

        # Ensure integer type and handle NaN
        for col in ['new_start_T1', 'new_start_T2', 'new_end_T1', 'new_end_T2']:
            if pd.notnull(row_modified[col]):
                row_modified[col] = int(row_modified[col])

        return row_modified

    # Split df2 into 8 chunks for parallel processing
    df2_split = np.array_split(df2, 8)

    # Function to process a single chunk
    def process_chunk(chunk, df1):
        return chunk.apply(lambda row: process_single_row(row.name, row, df1), axis=1)

    # Process each chunk in parallel
    df2_processed = Parallel(n_jobs=8)(
        delayed(process_chunk)(chunk, df1) for chunk in df2_split
    )

    # Combine results back into a single DataFrame
    df2_result = pd.concat(df2_processed).sort_index()

    return df2_result


# WHEN T!Strands are created
# REPRESENTING THE LOGCALLTY TRUE STRANDS wITH THE HEP OG FINAL TRANSRIPT AND GTF STRAND TO FIND BPTs
# Create ultimate strand for the transcripts to decide breakpoints
def transcript_strand(row):
    # Update T1strand based on conditions
    if row['final_transcript_strand_1'] == '-' and row['strand_x'] == '-':
        row['T1strand'] = '+'
    elif row['final_transcript_strand_1'] == '-' and row['strand_x'] == '+':
        row['T1strand'] = '-'
    elif row['final_transcript_strand_1'] == '+' and row['strand_x'] == '-':
        row['T1strand'] = '-'
    elif row['final_transcript_strand_1'] == '+' and row['strand_x'] == '+':
        row['T1strand'] = '+'

    # Update T2strand based on conditions
    if row['final_transcript_strand_2'] == '-' and row['strand_y'] == '-':
        row['T2strand'] = '+'
    elif row['final_transcript_strand_2'] == '-' and row['strand_y'] == '+':
        row['T2strand'] = '-'
    elif row['final_transcript_strand_2'] == '+' and row['strand_y'] == '-':
        row['T2strand'] = '-'
    elif row['final_transcript_strand_2'] == '+' and row['strand_y'] == '+':
        row['T2strand'] = '+'

    return row

# Function to determine bpT1 and bpT2 based on the strand information
def determine_bp(row):

    row['bpT1'] = row['end_T1']
    row['bpT2'] = row['start_T2']

    return row

def split_and_sort(range_str):
    start, end = map(int, range_str.split('-'))
    return min(start, end), max(start, end)

def switch_columns(row):
    if row['start_T1'] > row['end_T1']:
        row['start_T1'], row['end_T1'] = row['end_T1'], row['start_T1']
    if row['start_T2'] > row['end_T2']:
        row['start_T2'], row['end_T2'] = row['end_T2'], row['start_T2']
    return row

def bp_position_annotation(df1, df2, threshold=3):
    # Create a dictionary to map transcript IDs to their matching rows in df1
    transcript_dict = {
        'T1': df1[df1['transcript_id'].isin(df2['Transcript_ID_T1'])],
        'T2': df1[df1['transcript_id'].isin(df2['Transcript_ID_T2'])]
    }

    # Iterate over df2 rows
    for index2, row2 in df2.iterrows():
        # Get matching rows for T1 and T2 from the dictionary
        matching_rows_t1 = transcript_dict['T1'][transcript_dict['T1']['transcript_id'] == row2['Transcript_ID_T1']]
        matching_rows_t2 = transcript_dict['T2'][transcript_dict['T2']['transcript_id'] == row2['Transcript_ID_T2']]

        # Process T1 matches
        for index1, row1 in matching_rows_t1.iterrows():
            if abs(row2['bpT1'] - row1['new_exon_end']) <= threshold:
                if row2['start_T1'] >= row1['new_exon_start']:
                    df2.at[index2, 't1_region'] = "exon_ter"
                elif row2['start_T1'] < row1['new_exon_start'] < row2['bpT1']:
                    df2.at[index2, 't1_region'] = "exon_ter"  # "exon_ter/exon"
            elif row2['bpT1'] < row1['new_exon_end']:
                if row2['start_T1'] >= row1['new_exon_start']:
                    df2.at[index2, 't1_region'] = "exon_mid"
                elif row2['start_T1'] < row1['new_exon_start'] < row2['bpT1']:
                    df2.at[index2, 't1_region'] = "exon_mid"  # "exon_mid/exon"
                elif row2['bpT1'] == row1['new_exon_start']:
                    df2.at[index2, 't1_region'] = "exon_ter"

        # Process T2 matches
        for index1, row1 in matching_rows_t2.iterrows():
            if abs(row2['bpT2'] - row1['new_exon_start']) <= threshold:
                if row2['end_T2'] <= row1['new_exon_end']:
                    df2.at[index2, 't2_region'] = "exon_ter"
                else:
                    df2.at[index2, 't2_region'] = "exon_ter"  # "exon_ter/exon"
            elif row2['bpT2'] > row1['new_exon_start']:
                if row2['end_T2'] <= row1['new_exon_end'] or (row2['end_T2'] - row1['new_exon_end']) <= threshold:
                    df2.at[index2, 't2_region'] = "exon_mid"
                elif row2['end_T2'] > row1['new_exon_end'] and row2['start_T2'] < row1['new_exon_end']:
                    df2.at[index2, 't2_region'] = "exon_mid"  # "exon_mid/exon"
                elif row2['bpT2'] == row1['new_exon_end'] and row2['start_T2'] > row1['new_exon_start']:
                    df2.at[index2, 't2_region'] = "exon_ter"

    return df2

def bp_correction(df, output_path=None):

   # Apply conditions based on 'retain' column
    df.loc[df['retain'].str.contains('HEAD', case=False, na=False), 'cor_bpT1'] = df['new_bpT1'].astype(int)
    df.loc[df['retain'].str.contains('HEAD', case=False, na=False) & (df['overlap_query'] >= 0), 'cor_bpT2'] = (df['new_bpT2'] + df['overlap_query']).astype(int)

    df.loc[df['retain'].str.contains('TAIL', case=False, na=False), 'cor_bpT1'] = (df['new_bpT1'] - df['overlap_query']).astype(int)
    df.loc[df['retain'].str.contains('TAIL', case=False, na=False) & (df['overlap_query'] >= 0), 'cor_bpT2'] = df['new_bpT2'].astype(int)

    df.loc[df['retain'].str.contains('HEAD', case=False, na=False) & (df['gap_query'] >= 0), 'cor_bpT1'] = df['new_bpT1'].astype(int)
    df.loc[df['retain'].str.contains('HEAD', case=False, na=False) & (df['gap_query'] >= 0), 'cor_bpT2'] = (df['new_bpT2'] - df['gap_query']).astype(int)

    df.loc[df['retain'].str.contains('TAIL', case=False, na=False) & (df['gap_query'] >= 0), 'cor_bpT1'] = (df['new_bpT1'] + df['gap_query']).astype(int)
    df.loc[df['retain'].str.contains('TAIL', case=False, na=False) & (df['gap_query'] >= 0), 'cor_bpT2'] = df['new_bpT2'].astype(int)

        # Ensure the columns can hold string data
    df['cor_bpT1'] = df['cor_bpT1'].astype('object')
    df['cor_bpT2'] = df['cor_bpT2'].astype('object')
    #Apply conditions based on 't1_region' and 't2_region' columns and 'retain' column
    condition_exon_ter = df['t1_region'].str.contains('exon_ter', case=False, na=False) | df['t2_region'].str.contains('exon_ter', case=False, na=False)
    condition_retain_NONE_or_BOTH = df['retain'].str.contains('NONE|BOTH', case=False, na=False)

    #Apply the complex condition
    complex_condition = condition_exon_ter & condition_retain_NONE_or_BOTH

    def calculate_cor_bpT2(row, bp):
        if row['gap_query'] == 0 and row['overlap_query'] == 0:
            return int(row[bp]) # or any other logic you want for this case
        elif row['gap_query'] > 0 and row['overlap_query'] == 0:
            return int(row[bp] - row['gap_query'])
        elif row['overlap_query'] > 0 and row['gap_query'] == 0:
            return int(row[bp] + row['overlap_query'])
        else:
            return int(row[bp])

    def calculate_cor_bpT1(row, bp):
        if row['gap_query'] == 0 and row['overlap_query'] == 0:
            return int(row[bp])  # or any other logic you want for this case
        elif row['gap_query'] > 0 and row['overlap_query'] == 0:
            return int(row[bp] + row['gap_query'])
        elif row['overlap_query'] > 0 and row['gap_query'] == 0:
            return int(row[bp] - row['overlap_query'])
        else:
            return int(row[bp])

    def calculate_cor_bp(row, bp):
        if row['gap_query'] == 0 and row['overlap_query'] == 0:
            return int(row[bp])
        elif row['gap_query'] > 0 and row['overlap_query'] == 0:
            return str(str(row[bp]) + ' +/- ' + str(int(row['gap_query'])))
        elif row['overlap_query'] > 0 and row['gap_query'] == 0:
            return str(str(row[bp]) + ' +/- ' + str(int(row['overlap_query'])))
        else:
            return int(row[bp])

    #Apply conditions to modify 'cor_bpT1' and 'cor_bpT2' accordingly
    df.loc[complex_condition & df['t1_region'].str.contains('exon_ter', case=False, na=False), 'cor_bpT1'] = df['new_bpT1'].astype(int)
    df.loc[complex_condition & df['t1_region'].str.contains('exon_ter', case=False, na=False), 'cor_bpT2'] = df.apply(calculate_cor_bpT2, axis=1, bp='new_bpT2')

    df.loc[complex_condition & df['t2_region'].str.contains('exon_ter', case=False, na=False), 'cor_bpT2'] = df['new_bpT2'].astype(int)
    df.loc[complex_condition & df['t2_region'].str.contains('exon_ter', case=False, na=False), 'cor_bpT1'] = df.apply(calculate_cor_bpT1, axis=1, bp='new_bpT1')

    df.loc[condition_retain_NONE_or_BOTH & df['t1_region'].str.contains('exon_mid', case=False, na=False) & df['t2_region'].str.contains('exon_mid', case=False, na=False), 'cor_bpT1'] = df.apply(calculate_cor_bp, axis=1, bp='new_bpT1')
    df.loc[condition_retain_NONE_or_BOTH & df['t1_region'].str.contains('exon_mid', case=False, na=False) & df['t2_region'].str.contains('exon_mid', case=False, na=False), 'cor_bpT2'] = df.apply(calculate_cor_bp, axis=1, bp='new_bpT2')

    #Identify rows where both t1_region and t2_region are exon_ter  and exon_ter
    both_conditions = df['t1_region'].str.contains('exon_ter', case=False, na=False) & df['t2_region'].str.contains('exon_ter', case=False, na=False) & df['retain'].str.contains('NONE|BOTH', case=False, na=False)
    rows_to_duplicate = df[both_conditions]


     #   Apply the modifications directly to the existing DataFrame
    if not rows_to_duplicate.empty:

        df.loc[both_conditions, 'cor_bpT1'] = rows_to_duplicate['new_bpT1'].astype(int)
        for index, row in rows_to_duplicate.iterrows():
            if row['overlap_query'] > 0:
                df.at[index, 'cor_bpT2'] = int(row['new_bpT2'] + row['overlap_query'])
            elif row['gap_query'] > 0:
                df.at[index, 'cor_bpT2'] = int(row['new_bpT2'] + row['gap_query'])

        reversed_rows = rows_to_duplicate.copy()
        reversed_rows['cor_bpT2'] ==  reversed_rows['cor_bpT2'].astype(int)
        for index, row in reversed_rows.iterrows():
            if row['overlap_query'] > 0:
                df.at[index, 'cor_bpT1'] = int(row['new_bpT1'] + row['overlap_query'])
            elif row['gap_query'] > 0:
                df.at[index, 'cor_bpT1'] = int(row['new_bpT1'] + row['gap_query'])

     #    Append the reversed rows to the DataFrame
        df = pd.concat([df, reversed_rows], ignore_index=True)



    #If an output path is provided, save the modified DataFrame to a new CSV file
    if output_path:
        df.to_csv(output_path, index=False)

    return df

# Read files    
input_file = args.gtf
output1_file = args.transcript_info
output2_file = args.exons_position_adjusted

# Check if the output file already exists
if not os.path.exists(output1_file) or not os.path.exists(output2_file):
    # Read the GTF file into a DataFrame
    df1= pd.read_csv(input_file, sep='\t', header=None, skiprows=5, low_memory=False)

    # Filter rows where the 'feature' column value is "exon"
    transcript_df = df1[df1.iloc[:, 2] == "transcript"]

    print("Extracting tanscripts information from the gtf and modifying the eoxn position for next step...")
    # Transform the DataFrame
    transcript_df2 = transform_gtf_to_df(transcript_df)

    # Save the trancript_info.csv
    transcript_df2.to_csv(output1_file, index=False)

    exon_df = df1[df1.iloc[:, 2] == "exon"]

   # print(exon_df.columns)
    exon_df2 = transform_gtf_to_df(exon_df)
    #exon_df2.to_csv("exon_df2.csv")
    df3 = exon_df2[['chr', 'feature', 'strand','start', 'end', 'gene_id', 'transcript_id']]
    #df3.to_csv("df3.csv")
    #df3.loc[:, ["start", "end"]] = df3.loc[:, ["start", "end"]].rename(columns={"start": "exon_start", "end": "exon_end"})
    #df3.rename(columns= {"start": "exon_start", "end": "exon_end"},inplace=True)
    df4 = adjust_exon_positions_in_gtfdf(df3)

    df4.rename(columns= {"start": "exon_start", "end": "exon_end"},inplace=True)
    # Save the  exon_position_adjusted.csv
    df4.to_csv(output2_file, index=False)

else:
    # File already exists, so import the file and move to the next step
    df4 = pd.read_csv(output2_file, low_memory= False)
    transcript_df2 = pd.read_csv(output1_file, low_memory= False)

# Proceed with the next step : to get the genomic cordiantes , bp annoattion, bp correction for split 
# quantificaiton of aplit and span supporting read for chimera
# file import  allplit and fusion_df

allsplit = pd.read_csv(args.allsplit,delimiter='\t')
allsplit_anno =  pd.merge(allsplit,transcript_df2[['chr','feature','start','end','strand','gene_id', 'transcript_id', 'gene_name','gene_source', 'gene_biotype', 'transcript_name','transcript_source', 'transcript_biotype']], left_on='Transcript_ID_T1', right_on = 'transcript_id', how='inner')

##merge two dataframe with transcript_id as primary key
allsplit_anno2 = pd.merge(allsplit_anno,transcript_df2[['chr','feature','start','end','strand','gene_id', 'transcript_id', 'gene_name','gene_source', 'gene_biotype', 'transcript_name','transcript_source', 'transcript_biotype']], left_on='Transcript_ID_T2', right_on = 'transcript_id', how='inner')

# Apply the function to the DataFrame
allsplit_anno3 = allsplit_anno2.apply(transcript_strand, axis=1)
# allsplit_anno3.to_csv('garima.csv',sep='\t')

allsplit_anno3[['start_T1', 'end_T1']] = allsplit_anno3['subject_start_end_T1'].apply(lambda x: pd.Series(split_and_sort(x)))
allsplit_anno3[['start_T2', 'end_T2']] = allsplit_anno3['subject_start_end_T2'].apply(lambda x: pd.Series(split_and_sort(x)))
# Function to split and sort the start and end positions


# Apply the function to each row
allsplit_anno4 = allsplit_anno3.apply(determine_bp, axis=1)

# Display the desired columns
df5 = allsplit_anno4[['query_ID', 'chr_x', 'Transcript_ID_T1', 'gene_name_x', 'start_T1', 'final_transcript_strand_1','final_transcript_strand_2',
       'end_T1', 'gene_id_x','T1strand', 'strand_x',
       'chr_y', 'Transcript_ID_T2', 'gene_name_y', 'gene_id_y', 'start_T2',
       'end_T2', 'T2strand', 'strand_y',
       'overlap_query', 'gap_query', 'retain', 'bpT1', 'bpT2']]



# Switch columns efficiently
# Use a boolean mask to switch start and end columns
mask = df5['start_T1'] > df5['end_T1']
df5.loc[mask, ['start_T1', 'end_T1']] = df5.loc[mask, ['end_T1', 'start_T1']].values

mask = df5['start_T2'] > df5['end_T2']
df5.loc[mask, ['start_T2', 'end_T2']] = df5.loc[mask, ['end_T2', 'start_T2']].values

#print(df5)
#print(df5.columns)
# Ensure columns are of correct type
#df5['new_exon_start'] = df5['new_exon_start'].astype(int)                                       #..................not defined and not in df5 
#df5['new_exon_end'] = df5['new_exon_end'].astype(int)

# Use pd.to_numeric with errors='coerce' for numeric columns
numeric_columns = ['start_T1', 'end_T1', 'start_T2', 'end_T2', 'bpT1', 'bpT2']
df5[numeric_columns] = df5[numeric_columns].apply(pd.to_numeric, errors='coerce')

#df4.to_csv("df4.csv",sep='\t')
#df5.to_csv("df5.csv",sep='\t')
print("Extracting Transcript's Genomic coordinates for split reads...")
# Function to adjust transcript positions (assuming it's defined)
df6 = adjust_transcript_position_with_profiling(df4, df5)
print("Extraction of Transcript's Genomic coordinates is completed!\n")

print("Annotation of Breakpoints is running...")
df7 = bp_position_annotation(df4, df6)

#df7.to_csv("df7_adjusted_position.csv")
print("Annotaion of Breakpoints is completed!\n")

#define new breakpoint  based on genomic coordinates
df7['new_bpT1'] = df7['new_end_T1']
df7['new_bpT2'] = df7['new_start_T2']

#df7.to_csv("df7.csv",sep='\t',index=False)

print("Breakpoint correciton is running...")
# Calling bp_correction function
filtered_df = bp_correction(df7, 'SRR11279601_bp_corrected3.csv')
filtered_df = df7
print("Breakpoint correction is completed!\n")

# filtered_df.to_csv("bp_corrcetion.csv" ,sep=',', index=False)

# ft = pd.read_csv('/content/SRR11279601_bp_corrected3.csv', delimiter=',')
splitTab = filtered_df[['query_ID', 'chr_x', 'Transcript_ID_T1', 'gene_id_x','gene_name_x', 'strand_x', 'chr_y', 'Transcript_ID_T2', 'gene_id_y','gene_name_y',
       'strand_y', 'overlap_query', 'gap_query', 'retain', 'bpT1', 'bpT2', 'cor_bpT1', 'cor_bpT2', 'new_bpT1', 	'new_bpT2', 't1_region','t2_region']].sort_values(by=['query_ID', 'bpT1', 'bpT2'])

splitTab = splitTab.drop_duplicates()

#splitTab.to_csv("os_splitTab.csv",sep='\t')
del filtered_df
# splitTab.to_csv('collpasing.csv')
# Function to conditionally split the column and retain original value
def conditional_split(row):
    value = str(row['cor_bpT1'])  # Convert to string
    if '+' in value or '-' in value:
        nbp1, nobp1 = value.split('+', 1) if '+' in value else value.split('-', 1)
        return pd.Series([nbp1, nobp1, value])
    else:
        return pd.Series([value, None, value])

# Apply the function and create new columns
splitTab[['nbp1', 'nobp1', 'final_bpT1']] = splitTab.apply(conditional_split, axis=1)

# Function to conditionally split the column and retain original value
def conditional_split(row):
    value = str(row['cor_bpT2'])  # Convert to string
    if '+' in value or '-' in value:
        nbp1, nobp1 = value.split('+', 1) if '+' in value else value.split('-', 1)
        return pd.Series([nbp1, nobp1, value])
    else:
        return pd.Series([value, None, value])

# Apply the function and create new columns
splitTab[['nbp2', 'nobp2', 'final_bpT2']] = splitTab.apply(conditional_split, axis=1)
splitTab_1 = splitTab[['chr_x', 'Transcript_ID_T1', 'gene_id_x','gene_name_x', 'chr_y', 'Transcript_ID_T2','gene_id_y', 
       'gene_name_y', 'nbp1', 'nbp2', 'strand_x', 'strand_y','query_ID',
        't1_region', 't2_region']]

splitTab_1.rename(columns={"nbp1": "new_bpT1", "nbp2": "new_bpT2"}, inplace=True)
# new
#splitTab_1[['Transcript_ID_T1','no1']] = splitTab_1['Transcript_ID_T1'].str.split('.', expand=True)                     #..........used gene id instead split transcript id
#splitTab_1[['Transcript_ID_T2','no2']] = splitTab_1['Transcript_ID_T2'].str.split('.', expand=True)


pair_counts = splitTab_1.groupby(['query_ID','gene_id_x', 'gene_id_y']).size()

pairs_exactly_once = pair_counts[pair_counts == 1].index
splitTab_1_once= splitTab_1[splitTab_1.set_index(['query_ID','gene_id_x', 'gene_id_y']).index.isin(pairs_exactly_once)]

pairs_more_than_once = pair_counts[pair_counts > 1].index

# Filter the DataFrame for only the pairs that occur more than once
splitTab_1_filtered = splitTab_1[splitTab_1.set_index(['query_ID','gene_id_x', 'gene_id_y']).index.isin(pairs_more_than_once)]



# Step 2: Group by Transcript_ID_T1 and Transcript_ID_T2
grouped = splitTab_1_filtered.groupby(['Transcript_ID_T1', 'Transcript_ID_T2'])
grouped.head()

results = []

for (transcript1, transcript2), group in grouped:
    # Most frequent (mode) values of new_bpT1 and new_bpT2
    new_bpT1_modes = group['new_bpT1'].mode()
    new_bpT2_modes = group['new_bpT2'].mode()

    # Default values in case of empty modes
    if new_bpT1_modes.empty:
        new_bpT1_modes = pd.Series([0])  # or any other default value you'd like
    if new_bpT2_modes.empty:
        new_bpT2_modes = pd.Series([0])  # or any other default value you'd like

    # All combinations of the most frequent values
    for bpT1, bpT2 in product(new_bpT1_modes, new_bpT2_modes):
        results.append({
            'query_ID': group['query_ID'].iloc[0],
            'chr_x': group['chr_x'].iloc[0],
            'Transcript_ID_T1': transcript1,
            'gene_id_x': group['gene_id_x'].iloc[0],  
            'gene_name_x': group['gene_name_x'].iloc[0],
            'chr_y': group['chr_y'].iloc[0],
            'Transcript_ID_T2': transcript2,
            'gene_id_y': group['gene_id_y'].iloc[0],  
            'gene_name_y': group['gene_name_y'].iloc[0],
            'new_bpT1': bpT1,
            'new_bpT2': bpT2,
            'strand_x': group['strand_x'].iloc[0],
            'strand_y': group['strand_y'].iloc[0],
            't1_region': group['t1_region'].iloc[0],  
            't2_region': group['t2_region'].iloc[0],
        })

# Create the final DataFrame
output_df = pd.DataFrame(results)

RedundantSplit = pd.concat([splitTab_1_once, output_df], ignore_index=True)

RedundantSplit = RedundantSplit.drop_duplicates()
# RedundantSplit.to_csv('mod4spit', sep='\t', index=False)
# Convert to numeric values, coercing non-numeric values to NaN
RedundantSplit['new_bpT1'] = pd.to_numeric(RedundantSplit['new_bpT1'], errors='coerce')
# Fill NaN values with 0
RedundantSplit['new_bpT1'] = RedundantSplit['new_bpT1'].fillna(0)
# Round the values to the nearest integer and convert to integers
RedundantSplit['new_bpT1'] = RedundantSplit['new_bpT1'].apply(lambda x: int(round(x)))
RedundantSplit['new_bpT2'] = pd.to_numeric(RedundantSplit['new_bpT2'], errors='coerce')
# Fill NaN values with 0
RedundantSplit['new_bpT2'] = RedundantSplit['new_bpT2'].fillna(0)
# Round the values to the nearest integer and convert to integers
RedundantSplit['new_bpT2'] = RedundantSplit['new_bpT2'].apply(lambda x: int(round(x)))


# filtered_df = RedundantSplit
# sp1 = filtered_df[[ 'chr_x', 'Transcript_ID_T1',  'gene_name_x', 'chr_y', 'Transcript_ID_T2', 'gene_name_y','new_bpT1',  'new_bpT2','strand_x', 'strand_y' ]].value_counts().reset_index()
df_with_splitcount = RedundantSplit[[ 'chr_x', 'Transcript_ID_T1','gene_id_x','gene_name_x', 'chr_y', 'Transcript_ID_T2','gene_id_y' ,'gene_name_y','new_bpT1',  'new_bpT2','strand_x', 'strand_y','t1_region','t2_region' ]].value_counts(dropna=False).reset_index()

# sp1.to_csv('mod4_sp1_0.5.csv',sep='\t', index=False)
df_with_splitcount = df_with_splitcount.rename(columns={'count': 'Split_support_fbpt'})
# THIS NEDS TO BE CORRECTED BECAUSE SO THAT IT WILL TAKE VALUES TAKEN AFTER CONSIDERING OVERLAPS AND GAPS======================lOOKOUT!!!!!                    # DONE
# sp1 =  nft1[[ 'chr_x', 'Transcript_ID_T1', 'gene_name_x', 'chr_y', 'Transcript_ID_T2', 'gene_name_y']].value_counts().reset_index()
df_with_splitcount = df_with_splitcount.rename(columns={'count': 'Split_support_fbpt'})
#df_with_splitcount.to_csv("df_with_splitcount.csv",sep='\t')

# READ spanning_fsuion_df.csv
fusion_df = pd.read_csv(args.allsplit_fusiondf,delimiter='\t')
fusion_df.rename(columns={"Start1": "start_T1", "End1": "end_T1", "Start2": "start_T2", "End2": "end_T2"}, inplace=True)
# filterin only the rows of span_fusiondf whose t1 and t2 id is present in the df_with_splitcount dataframe
# Create a set of tuples for the pairs in df_with_splitcount
mask = fusion_df[['Transcript_ID_T1', 'Transcript_ID_T2']].apply(tuple, axis=1).isin(df_with_splitcount[['Transcript_ID_T1', 'Transcript_ID_T2']].apply(tuple, axis=1))
# Filter fusion_df based on the mask
filtered_fusion_df = fusion_df[mask]

#transcript_df2 is transcript_info.txt file 
#annotation of transripts of filtered_fusion_df
#merge two dataframe with trnacript_id as primary key
span_anno =  pd.merge(filtered_fusion_df,transcript_df2[['chr','feature','start','end','strand','gene_id', 'transcript_id', 'gene_name','gene_source', 'gene_biotype', 'transcript_name','transcript_source', 'transcript_biotype']], left_on='Transcript_ID_T1', right_on = 'transcript_id', how='inner')
del filtered_fusion_df
del fusion_df
#merge two dataframe with transcript_id as primary key
df5 = pd.merge(span_anno,transcript_df2[['chr','feature','start','end','strand','gene_id', 'transcript_id', 'gene_name','gene_source', 'gene_biotype', 'transcript_name','transcript_source', 'transcript_biotype']], left_on='Transcript_ID_T2', right_on = 'transcript_id', how='inner')

df5 = df5.apply(switch_columns, axis=1)  # switch column when start>end

# Ensure columns are of correct type
df4['new_exon_start'] = df4['new_exon_start'].astype(int)
df4['new_exon_end'] = df4['new_exon_end'].astype(int)
numeric_columns = ['start_T1', 'end_T1', 'start_T2', 'end_T2']
df5[numeric_columns] = df5[numeric_columns].apply(pd.to_numeric, errors='coerce')


#df5.to_csv("s272_df5.csv",sep='\t')
#print(df5)

print("Extracting Transcript's Genomic corodinates for spanning reads...")
# tion to change transcriptomic ordinates to genomic coordinates
df1 = adjust_transcript_position_with_profiling(df4, df5)


print("Extraction of Transcript's Genomic corodinates is completed!\n")


table1 = df_with_splitcount
#df1.to_csv("s272_df1.csv",sep='\t')
#table1.to_csv("s272_table1.csv",sep='\t')
#this step need modification to get the span read counts
# Initialize new columns in table1 to hold the counts

    # Process tables
if df1.empty:
    print("No chimera detected for this sample.")
else:
    # Initialize new columns in table1 to hold the counts
    table1['span_read_count'] = 0
    table1['new_bpT1'] = table1['new_bpT1'].astype(float).astype(int)
    table1['new_bpT2'] = table1['new_bpT2'].astype(float).astype(int)

    # Iterate over table1 and apply the counting logic
    for idx, row in table1.iterrows():
        count1 = 0
        count2 = 0

        # Strand information
        strand_x = row['strand_x']
        strand_y = row['strand_y']

        # Define condition1
        condition1 = (
            ((df1['Transcript_ID_T1'] == row['Transcript_ID_T1']) & 
            (df1['strand_x'] == row['strand_x']) &  # Match strand_x for Transcript_ID_T1
            (df1['Transcript_ID_T2'] == row['Transcript_ID_T2']) & 
            (df1['strand_y'] == row['strand_y'])) &  # Match strand_y for Transcript_ID_T2
            (
                (strand_x == '-' and strand_y == '-' and
                (df1['new_start_T1'] >= row['new_bpT1']) & 
                (df1['new_end_T1'] >= row['new_bpT1']) &
                (df1['new_start_T2'] <= row['new_bpT2']) & 
                (df1['new_end_T2'] <= row['new_bpT2'])) |
                
                (strand_x == '-' and strand_y == '+' and
                (df1['new_start_T1'] >= row['new_bpT1']) & 
                (df1['new_end_T1'] >= row['new_bpT1']) &
                (df1['new_start_T2'] >= row['new_bpT2']) & 
                (df1['new_end_T2'] >= row['new_bpT2'])) |

                (strand_x == '+' and strand_y == '+' and
                (df1['new_start_T1'] <= row['new_bpT1']) & 
                (df1['new_end_T1'] <= row['new_bpT1']) &
                (df1['new_start_T2'] >= row['new_bpT2']) & 
                (df1['new_end_T2'] >= row['new_bpT2'])) |

                (strand_x == '+' and strand_y == '-' and
                (df1['new_start_T1'] <= row['new_bpT1']) & 
                (df1['new_end_T1'] <= row['new_bpT1']) &
                (df1['new_start_T2'] <= row['new_bpT2']) & 
                (df1['new_end_T2'] <= row['new_bpT2']))
            )
        )

        # Define condition2
        condition2 = (
            ((df1['Transcript_ID_T1'] == row['Transcript_ID_T2']) & 
            (df1['strand_x'] == row['strand_y']) &  # Match strand_y for Transcript_ID_T1
            (df1['Transcript_ID_T2'] == row['Transcript_ID_T1']) & 
            (df1['strand_y'] == row['strand_x'])) &  # Match strand_x for Transcript_ID_T2
            (
                (strand_x == '-' and strand_y == '-' and
                (df1['new_start_T2'] >= row['new_bpT1']) & 
                (df1['new_end_T2'] >= row['new_bpT1']) &
                (df1['new_start_T1'] <= row['new_bpT2']) & 
                (df1['new_end_T1'] <= row['new_bpT2'])) |

                (strand_x == '-' and strand_y == '+' and
                (df1['new_start_T2'] >= row['new_bpT1']) & 
                (df1['new_end_T2'] >= row['new_bpT1']) &
                (df1['new_start_T1'] >= row['new_bpT2']) & 
                (df1['new_end_T1'] >= row['new_bpT2'])) |

                (strand_x == '+' and strand_y == '+' and
                (df1['new_start_T2'] <= row['new_bpT1']) & 
                (df1['new_end_T2'] <= row['new_bpT1']) &
                (df1['new_start_T1'] >= row['new_bpT2']) & 
                (df1['new_end_T1'] >= row['new_bpT2'])) |

                (strand_x == '+' and strand_y == '-' and
                (df1['new_start_T2'] <= row['new_bpT1']) & 
                (df1['new_end_T2'] <= row['new_bpT1']) &
                (df1['new_start_T1'] <= row['new_bpT2']) & 
                (df1['new_end_T1'] <= row['new_bpT2']))
            )
        )

        #Extract Query_id values for rows satisfying the conditions
        satisfied_rows_condition1 = df1[condition1]
        satisfied_rows_condition2 = df1[condition2]

        # Combine satisfied rows and get unique Query_ids
        unique_query_ids = pd.concat([satisfied_rows_condition1, satisfied_rows_condition2])['query_ID'].unique()

        # Count the unique Query_ids
        unique_query_count = len(unique_query_ids)

        # Update span_read_count with the unique Query_id count for the fusion pair in table1
        table1.at[idx, 'span_read_count_uniq_queryid_count'] = unique_query_count


        # Sum counts for the main conditions
        count1 += condition1.sum() if isinstance(condition1, pd.Series) else int(condition1)
        count2 += condition2.sum() if isinstance(condition2, pd.Series) else int(condition2)

        #table1.at[idx, 'count1'] = count1
        #table1.at[idx, 'count2'] = count2

        # Update counts in table1
        table1.at[idx, 'span_read_count'] = count1 + count2


    # Save the updated table
    table1.to_csv(args.ftable, sep='\t', index=False)



'''
#clean up memory
del allsplit_anno4
del allsplit_anno3
del allsplit_anno2
del allsplit_anno
del df7 
del df6
del df5 
del df4
del table1 
del df1
del fusion_df
del new_bpT1_modes
del new_bpT2_modes
del RedundantSplit
del exon_df2
del pairs_exactly_once
del transcript1
del transcript2
del transcript_df
del transcript_df2
del RedundantSplit
del splitTab
del splitTab_1
del splitTab_1_filtered
del splitTab_1_once
'''

end_time = time.time()
usage = resource.getrusage(resource.RUSAGE_SELF)
max_memory = usage.ru_maxrss  # Memory in kilobytes
print(f"Execution Time: {end_time - start_time} seconds")
print(f"Max Memory Usage: {max_memory / 1024} MB")
