
#############################################################################################################################
#                                      Discordant Read Analysis: Split and Spanning Discovery
#############################################################################################################################

import re
import glob, os
import argparse
import math
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings("ignore", category=pd.errors.SettingWithCopyWarning)
import resource
import time

start_time = time.time()

parser = argparse.ArgumentParser(description='blast output analysis')
parser.add_argument( 'df_blast')
parser.add_argument('df_SAfile')
parser.add_argument('fusion_df')
parser.add_argument('df_blast_output')
args= parser.parse_args()


# Read the CSV files
df_blast = pd.read_csv(args.df_blast, delimiter='\t')
df_SAfile = pd.read_csv(args.df_SAfile, sep="\t")
df_SAfile['read_id'] = df_SAfile.groupby('read_id').cumcount().astype(str) + '#' + df_SAfile['read_id']
# bringing the coordinates in a single column
df_blast['query_start_end'] = df_blast['query_start'].astype(str) + "-" + df_blast['query_end'].astype(str)
df_blast['subject_start_end'] = df_blast['subject_start'].astype(str) + "-" + df_blast['subject_end'].astype(str)
df_blast.drop(columns=['query_start', 'query_end', 'subject_start', 'subject_end'], inplace=True)
df_blast[['query_ID', 'frame_id1']] = df_blast['query_ID'].str.split('_', expand=True)
# Create merge keys
df_SAfile['merge_key_T1'] = df_SAfile['read_id'] + '_' + df_SAfile['transcript1_id']
df_SAfile['merge_key_T2'] = df_SAfile['read_id'] + '_' + df_SAfile['transcript2_id']
df_blast['merge_key_T1'] = df_blast['query_ID'] + '_' + df_blast['Transcript_ID']
df_blast['merge_key_T2'] = df_blast['query_ID'] + '_' + df_blast['Transcript_ID']

df_SAfile['SAfile_index'] = df_SAfile.index

df_SAfile['SAfile_index'] = df_SAfile['SAfile_index'].apply(lambda x: x if isinstance(x, list) else [x])

# Format the SAfile_index with prefixes
df_SAfile['SAfile_index_formatted_T1'] = df_SAfile['SAfile_index'].apply(lambda x: [f'T1_{i}' for i in x])
df_SAfile['SAfile_index_formatted_T2'] = df_SAfile['SAfile_index'].apply(lambda x: [f'T2_{i}' for i in x])

# Create dictionaries for fast lookup
index_T1_dict = df_SAfile.set_index('merge_key_T1')['SAfile_index_formatted_T1'].to_dict()
index_T2_dict = df_SAfile.set_index('merge_key_T2')['SAfile_index_formatted_T2'].to_dict()

 #Map the formatted indices to df_blast and handle missing values
df_blast['SAfile_index_T1'] = df_blast['merge_key_T1'].map(index_T1_dict)
df_blast['SAfile_index_T2'] = df_blast['merge_key_T2'].map(index_T2_dict)

# Replace NaN values with empty lists
df_blast['SAfile_index_T1'] = df_blast['SAfile_index_T1'].apply(lambda x: x if isinstance(x, list) else [])
df_blast['SAfile_index_T2'] = df_blast['SAfile_index_T2'].apply(lambda x: x if isinstance(x, list) else [])
df_blast['SA_row_num']= df_blast['SAfile_index_T1'] + df_blast['SAfile_index_T2']
df_blast = df_blast.explode('SA_row_num')
df_blast['query_ID'] = df_blast['query_ID'] + "_" + df_blast['frame_id1']
ab = df_blast[['query_ID','Transcript_ID', 'strand', 'expect', 'identity', 'query_start_end', 'subject_start_end', 'SA_row_num']]
df_blast=ab
df_blast = df_blast.groupby(['query_ID', 'SA_row_num']).agg({'Transcript_ID': lambda x: list(map(str, x)), \
                                                'query_start_end': lambda x: list(map(str, x)), \
                                                'subject_start_end': lambda x: list(map(str, x)), \
                                                'strand': lambda x: list(map(str,x)), \
                                                'expect': lambda x: list(map(str, x)), \
                                                'identity': lambda x: list(map(str, x))}).reset_index()


df_blast[['SA_row_num', 'row_id']] = df_blast['SA_row_num'].str.split('_', expand=True)
del ab
# df_blast['query_start_end'] = df_blast['query_start'].astype(str) + "-" + df_blast['query_end'].astype(str)
# df_blast['subject_start_end'] = df_blast['subject_start'].astype(str) + "-" + df_blast['subject_end'].astype(str)
# df_blast.drop(columns=['query_start', 'query_end', 'subject_start', 'subject_end'], inplace=True)
# df_blast[['query_ID', 'frame_id1']] = df_blast['query_ID'].str.split('_', expand=True)
# # To trace back to which index pos in SAfile.txt transcript in df_blast belong to
# df_blast["SA_row_num1"] = df_blast.apply(lambda x: df_SAfile.index[(df_SAfile["read_id"] == x["query_ID"]) & (df_SAfile["transcript1_id"] == x["Transcript_ID"])].to_list(), axis=1)
# df_blast["SA_row_num2"] = df_blast.apply(lambda x: df_SAfile.index[(df_SAfile["read_id"] == x["query_ID"]) & (df_SAfile["transcript2_id"] == x["Transcript_ID"])].to_list(), axis=1)
# df_blast["SA_row_num1"] = df_blast["SA_row_num1"].apply(lambda x: ["T1_" + str(i) for i in x])
# df_blast["SA_row_num2"] = df_blast["SA_row_num2"].apply(lambda x: ["T2_" + str(i) for i in x])
# df_blast["SA_row_num"] = df_blast["SA_row_num1"]+df_blast["SA_row_num2"]
# df_blast = df_blast.explode("SA_row_num")
# df_blast['query_ID'] = df_blast['query_ID'] + "_" + df_blast['frame_id1']
# df_blast.drop(columns=['SA_row_num1', 'SA_row_num2', 'frame_id1'], inplace=True)
# df_blast = df_blast.groupby(['query_ID', 'SA_row_num']).agg({'Transcript_ID': lambda x: list(map(str, x)), \
#                                                 'query_start_end': lambda x: list(map(str, x)), \
#                                                 'subject_start_end': lambda x: list(map(str, x)), \
#                                                 'strand': lambda x: list(map(str,x)), \
#                                                 'expect': lambda x: list(map(str, x)), \
#                                                 'identity': lambda x: list(map(str, x))}).reset_index()
# df_blast[['SA_row_num', 'row_id']] = df_blast['SA_row_num'].str.split('_', expand=True)

df_blast.to_csv("cicer_df_blast.csv",sep='\t')
df_blast['row_id'] = pd.to_numeric(df_blast['row_id'])
df_blast.sort_values(['row_id'], inplace=True)
df_blast[["query_ID","SA_row_num", "row_id", "Transcript_ID"]]
df_blast_T1 = df_blast[df_blast['SA_row_num'] == 'T1']
df_blast_T2 = df_blast[df_blast['SA_row_num'] == 'T2']
df_blast_T1.drop(columns=['row_id', 'SA_row_num'], inplace=True)
df_blast_T2.drop(columns=['row_id', 'SA_row_num'], inplace=True)
df_blast = pd.merge(df_blast_T1, df_blast_T2, on='query_ID', how='outer', suffixes=('_T1', '_T2'))
df_blast[['query_ID', 'frame_id1']] = df_blast['query_ID'].str.split('_', expand=True)
df_blast.loc[df_blast['Transcript_ID_T1']==df_blast['Transcript_ID_T2']]
df_blast_uniq = df_blast[~df_blast.duplicated(['query_ID'], keep=False)]
df_blast = df_blast[df_blast.duplicated(['query_ID'], keep=False)]
df_blast['query_ID'] = df_blast['query_ID'] + "_" + df_blast['frame_id1']
df_blast_uniq['query_ID'] = df_blast_uniq['query_ID'] + "_" + df_blast_uniq['frame_id1']
df_blast.drop(columns=['frame_id1'], inplace=True)
df_blast_uniq.drop(columns=['frame_id1'], inplace=True)
df_blast.sort_values(['query_ID'], inplace=True)
df_blast_uniq.sort_values(['query_ID'], inplace=True)
df_blast = df_blast.explode(["Transcript_ID_T1","query_start_end_T1","subject_start_end_T1","strand_T1","expect_T1","identity_T1"])
df_blast = df_blast.explode(["Transcript_ID_T2","query_start_end_T2","subject_start_end_T2","strand_T2","expect_T2","identity_T2"])
df_blast_uniq = df_blast_uniq.explode(["Transcript_ID_T1","query_start_end_T1","subject_start_end_T1","strand_T1","expect_T1","identity_T1"])
df_blast_uniq = df_blast_uniq.explode(["Transcript_ID_T2","query_start_end_T2","subject_start_end_T2","strand_T2","expect_T2","identity_T2"])
df_blast_NaNs = df_blast[df_blast.isna().any(axis=1)]
df_blast.dropna(inplace=True)
df_blast[['query_ID', 'frame_id1']] = df_blast['query_ID'].str.split('_', expand=True)
df_blast_2 = df_blast[~df_blast.duplicated(['query_ID'], keep=False)]
df_blast = df_blast[df_blast.duplicated(['query_ID'], keep=False)]
# This line creates a DataFrame df_blast_f1 containing rows from df_blast where the 'frame_id1' column is '1'.
df_blast_f1 = df_blast[df_blast['frame_id1'] == '1']
df_blast_f1
# Similarly, this line creates a DataFrame df_blast_f2 containing rows from df_blast where the 'frame_id1' column is '2'.
df_blast_f2 = df_blast[df_blast['frame_id1'] == '2']
df_blast_f2
# This line creates a list of 'query_ID' values from df_blast_f1 that are not present in df_blast_f2.
df_blast_f1_list = df_blast_f1[~df_blast_f1['query_ID'].isin(df_blast_f2['query_ID'])]['query_ID'].tolist()
df_blast_f2_list = df_blast_f2[~df_blast_f2['query_ID'].isin(df_blast_f1['query_ID'])]['query_ID'].tolist()
df_blast_3 = df_blast[df_blast['query_ID'].isin(df_blast_f1_list + df_blast_f2_list)]
df_blast['query_ID'] = df_blast['query_ID'] + "_" + df_blast['frame_id1']
df_blast_2['query_ID'] = df_blast_2['query_ID'] + "_" + df_blast_2['frame_id1']
df_blast_3['query_ID'] = df_blast_3['query_ID'] + "_" + df_blast_3['frame_id1']
df_blast.drop(columns=['frame_id1'], inplace=True)
df_blast_2.drop(columns=['frame_id1'], inplace=True)
df_blast_3.drop(columns=['frame_id1'], inplace=True)
# Removes rows from df_blast where the 'query_ID' is present in the 'query_ID' column of df_blast_3.
df_blast = df_blast[~df_blast['query_ID'].isin(df_blast_3['query_ID'])]
# Concatenates df_blast_uniq, df_blast_NaNs, df_blast_2, and df_blast_3 along the row axis.
df_blast_uniq = pd.concat([df_blast_uniq, df_blast_NaNs, df_blast_2, df_blast_3], axis=0)
# Sorts df_blast_uniq by the 'query_ID' column.
df_blast_uniq.sort_values(['query_ID'], inplace=True)
# Splits the 'query_ID' column of df_blast_uniq by underscore and assigns the result to new columns 'query_ID' and 'frame_id1'.
df_blast_uniq[['query_ID', 'frame_id1']] = df_blast_uniq['query_ID'].str.split('_', expand=True)
# Groups df_blast_uniq by 'query_ID' and keeps only groups with a length of 2.
df_blast_uniq_gp1 = df_blast_uniq.groupby(['query_ID']).filter(lambda x: len(x) == 2)
# Creates a list of 'query_ID' values from df_blast_uniq_gp1 where no NaN values are present.
df_blast_uniq_gp1_list = df_blast_uniq_gp1[~df_blast_uniq_gp1.isna().any(axis=1)]['query_ID'].tolist()
df_blast_uniq_gp1 = df_blast_uniq_gp1[df_blast_uniq_gp1['query_ID'].isin(df_blast_uniq_gp1_list)]
df_blast_uniq_gp1['query_ID'] = df_blast_uniq_gp1['query_ID'] + "_" + df_blast_uniq_gp1['frame_id1']
df_blast_uniq['query_ID'] = df_blast_uniq['query_ID'] + "_" + df_blast_uniq['frame_id1']
df_blast_uniq_gp1.drop(columns=['frame_id1'], inplace=True)
df_blast_uniq.drop(columns=['frame_id1'], inplace=True)
# OKAY, HERE IVE MY spnaning READS
df_blast_uniq = df_blast_uniq[~df_blast_uniq['query_ID'].isin(df_blast_uniq_gp1['query_ID'])]
df_blast_uniq[['query_ID', 'frameID']] = df_blast_uniq['query_ID'].str.split('_',expand=True)
# keeping only those that are in a pair
df_blast_uniq_double = df_blast_uniq[df_blast_uniq.duplicated(['query_ID'], keep=False)]
df_blast_uniq_double_1 = df_blast_uniq_double.loc[df_blast_uniq_double['frameID']=='1']
df_blast_uniq_double_2 = df_blast_uniq_double.loc[df_blast_uniq_double['frameID']=='2']
df_spanning = pd.merge(df_blast_uniq_double_1,df_blast_uniq_double_2, on = "query_ID", suffixes=('_fd','_rv'))
# df_spanning
# Shifting towards left

def shift_left_and_fill_nan(row):
    non_nans = row.dropna().tolist()
    return pd.Series(non_nans + [np.nan] * (len(row) - len(non_nans)))
# Apply the function to each row
shifted_df_spanning = df_spanning.apply(shift_left_and_fill_nan, axis=1)
# Keepin only those reads that are mapping to two different transcripts
trim_shifted_df_spanning=shifted_df_spanning.loc[shifted_df_spanning[1]!=shifted_df_spanning[8]]
ab4 = shifted_df_spanning.loc[(shifted_df_spanning[7]=='1')& (shifted_df_spanning[14]=='2')]
ab5 = ab4.loc[ab4[1]!=ab4[8]]
gk_span = ab5[[1,2,3,4,5,6,8,9,10,11,12,13]]
# Function to split and sort dash-separated values
def split_and_sort(range_str):
    start, end = map(int, range_str.split('-'))
    return min(start, end), max(start, end)
# Split and sort the columns
gk_span[['Start1', 'End1']] = gk_span[3].apply(lambda x: pd.Series(split_and_sort(x)))
gk_span[['Start2', 'End2']] = gk_span[10].apply(lambda x: pd.Series(split_and_sort(x)))
gk_span.rename(columns={1:"Transcript_ID_T1", 8:"Transcript_ID_T2"},inplace=True)
gk_span_save = gk_span[['Transcript_ID_T1','Start1','End1','Transcript_ID_T2', 'Start2', 'End2']]
# gk_span_save = gk_span_save.to_csv('fusion_df.csv',sep='\t',index=False)
gk_span_save = gk_span_save.to_csv(args.fusion_df,sep='\t',index=False)


ab6 = shifted_df_spanning.loc[~((shifted_df_spanning[7]=='1')& (shifted_df_spanning[14]=='2'))]
df_blast = pd.concat([df_blast, df_blast_uniq_gp1], axis=0)
id =  ab6[0].to_list()
extra_sp = df_blast_uniq_double.loc[df_blast_uniq_double['query_ID'].isin(id)]
extra_sp['query_ID'] = extra_sp['query_ID']+ '_' + extra_sp['frameID']
extra_sp = extra_sp.drop('frameID',axis=1)
df_blast = pd.concat([df_blast,extra_sp])
df_blast['strand_symbol_1'] = df_blast['strand_T1'].str.split('/').str[1].replace({'Plus': '+', 'Minus': '-', 'NaN': ' '})
df_blast['strand_symbol_2'] = df_blast['strand_T2'].str.split('/').str[1].replace({'Plus': '+', 'Minus': '-', 'NaN': ' '})
df_blast_gp1 = df_blast.groupby(['query_ID']).first().reset_index()
df_blast_gp1 = df_blast.groupby(['query_ID']).first().reset_index()
df_blast_gp1_fd = df_blast_gp1[df_blast_gp1['query_ID'].str.endswith('_1')]
df_blast_gp1_rv = df_blast_gp1[df_blast_gp1['query_ID'].str.endswith('_2')]
df_blast_gp1_fd[['query_ID', 'frame_id1']] = df_blast_gp1_fd['query_ID'].str.split('_', expand=True)
df_blast_gp1_rv[['query_ID', 'frame_id1']] = df_blast_gp1_rv['query_ID'].str.split('_', expand=True)
df_blast_gp1 = pd.merge(df_blast_gp1_fd, df_blast_gp1_rv, on='query_ID', how='outer', suffixes=('_fd', '_rv'))

#### -------------------------------------------Final transcript strand addition------------
# drop row if Transcript_Id_fd and rv, both are NaN
condition = ((df_blast_gp1['Transcript_ID_T1_fd'].isna()) & (df_blast_gp1['Transcript_ID_T2_fd'].isna()) | (df_blast_gp1['Transcript_ID_T1_rv'].isna()) & (df_blast_gp1['Transcript_ID_T2_rv'].isna()))
df_blast_gp1 = df_blast_gp1[~condition]
def assignStrand_fd(x):
    if (x['strand_symbol_1_fd'] == '+' and x['strand_symbol_1_rv'] == '-'):
        return '+'
    elif (x['strand_symbol_1_fd'] == '-' and x['strand_symbol_1_rv'] == '+'):
        return '-'
    elif (x['strand_symbol_1_fd'] == '+' and x['strand_symbol_1_rv'] is None):
        return '+'
    elif (x['strand_symbol_1_fd'] == '-' and x['strand_symbol_1_rv'] is None):
        return '-'
    elif (x['strand_symbol_1_fd'] is None and x['strand_symbol_1_rv'] == '-'):
        return '+'
    elif (x['strand_symbol_1_fd'] is None and x['strand_symbol_1_rv'] == '+'):
        return '-'
    else:
        return 'No value'

def assignStrand_rv(x):
    if (x['strand_symbol_2_fd'] == '+' and x['strand_symbol_2_rv'] == '-'):
        return '+'
    elif (x['strand_symbol_2_fd'] == '-' and x['strand_symbol_2_rv'] == '+'):
        return '-'
    elif (x['strand_symbol_2_fd'] == '+' and x['strand_symbol_2_rv'] is None):
        return '+'
    elif (x['strand_symbol_2_fd'] == '-' and x['strand_symbol_2_rv'] is None):
        return '-'
    elif (x['strand_symbol_2_fd'] is None and x['strand_symbol_2_rv'] == '+'):     # CHECK THIS # Done checking changed to -
        return '-'
    elif (x['strand_symbol_2_fd'] is None and x['strand_symbol_2_rv'] == '-'):
        return '+'
    else:
        return 'No value'

df_blast_gp1['final_transcript_strand_1'] = df_blast_gp1.apply(assignStrand_fd, axis=1)
df_blast_gp1['final_transcript_strand_2'] = df_blast_gp1.apply(assignStrand_rv, axis=1)
df_blast_gp1.loc[(df_blast_gp1['strand_T1_fd']=='Plus/Plus')&(df_blast_gp1['strand_T2_fd']=='Plus/Minus')]
df_blast_gp1['final_transcript_strand'] = df_blast_gp1[['final_transcript_strand_1', 'final_transcript_strand_2']].values.tolist()
condition = ((df_blast_gp1['Transcript_ID_T1_fd'].isna()) & (df_blast_gp1['Transcript_ID_T2_fd'].isna()) | (df_blast_gp1['Transcript_ID_T1_rv'].isna()) & (df_blast_gp1['Transcript_ID_T2_rv'].isna()))
df_blast_gp1 = df_blast_gp1[~condition]
df_blast_gp1['final_frame_id'] = df_blast_gp1[['frame_id1_fd', 'frame_id1_rv']].values.tolist()
df_blast_gp1 = df_blast_gp1.explode(["final_transcript_strand", "final_frame_id"])
df_blast_gp1['query_ID'] = df_blast_gp1['query_ID'] + "_" + df_blast_gp1['final_frame_id']
df_blast_gp1['Transcript_ID_T1_fd'] = df_blast_gp1.apply(lambda x: None if x['Transcript_ID_T1_fd'] is None else (x['Transcript_ID_T1_fd'] + "_" + x['final_transcript_strand'] if x['query_ID'].endswith('_1') else ''), axis=1)
df_blast_gp1['Transcript_ID_T2_fd'] = df_blast_gp1.apply(lambda x: None if x['Transcript_ID_T2_fd'] is None else (x['Transcript_ID_T2_fd'] + "_" + x['final_transcript_strand'] if x['query_ID'].endswith('_2') else ''), axis=1)
df_blast_gp1['Transcript_ID_T1_rv'] = df_blast_gp1.apply(lambda x: None if x['Transcript_ID_T1_rv'] is None else (x['Transcript_ID_T1_rv'] + "_" + x['final_transcript_strand'] if x['query_ID'].endswith('_1') else ''), axis=1)
df_blast_gp1['Transcript_ID_T2_rv'] = df_blast_gp1.apply(lambda x: None if x['Transcript_ID_T2_rv'] is None else (x['Transcript_ID_T2_rv'] + "_" + x['final_transcript_strand'] if x['query_ID'].endswith('_2') else ''), axis=1)
df_blast_gp1['Transcript_ID_T1_final'] = df_blast_gp1[['Transcript_ID_T1_fd', 'Transcript_ID_T1_rv']].values.tolist()
df_blast_gp1['Transcript_ID_T2_final'] = df_blast_gp1[['Transcript_ID_T2_fd', 'Transcript_ID_T2_rv']].values.tolist()
df_blast_gp1['Transcript_ID_T1_final'] = df_blast_gp1['Transcript_ID_T1_final'].apply(lambda x: set(list(filter(None, x))))
df_blast_gp1['Transcript_ID_T2_final'] = df_blast_gp1['Transcript_ID_T2_final'].apply(lambda x: set(list(filter(None, x))))
df_blast_gp1 = df_blast_gp1.explode(["Transcript_ID_T1_final"])
df_blast_gp1 = df_blast_gp1.explode(["Transcript_ID_T2_final"])
df_blast_gp1.drop(columns=['frame_id1_fd', 'frame_id1_rv', 'final_frame_id'], inplace=True)
df_blast_gp1[['query_ID', 'final_frame_id']] = df_blast_gp1['query_ID'].str.split('_', expand=True)
df_blast_gp1_short = df_blast_gp1.groupby(['query_ID']).agg({'Transcript_ID_T1_final': lambda x: list(map(str, x)),'Transcript_ID_T2_final': lambda x: list(map(str, x))}).reset_index()
df_blast_gp1_short['Transcript_ID_T1_final'] = df_blast_gp1_short['Transcript_ID_T1_final'].apply(lambda x: list(filter(lambda y: len(y) > 5, x)))                 # WHAT's this
df_blast_gp1_short['Transcript_ID_T2_final'] = df_blast_gp1_short['Transcript_ID_T2_final'].apply(lambda x: list(filter(lambda y: len(y) > 5, x)))
df_blast_gp1_short['Transcript_ID_T1_final'] = df_blast_gp1_short['Transcript_ID_T1_final'].apply(lambda x: ','.join(x))
df_blast_gp1_short['Transcript_ID_T2_final'] = df_blast_gp1_short['Transcript_ID_T2_final'].apply(lambda x: ','.join(x))
df_blast_gp1['query_ID'] = df_blast_gp1['query_ID'] + "_" + df_blast_gp1['final_frame_id']
df_blast['final_transcript_strand'] = df_blast['query_ID'].map(df_blast_gp1.set_index('query_ID')['final_transcript_strand'])
df_blast.to_csv(args.df_blast_output,sep='\t',index=False)

#clean up memory
del ab4
del ab5
# del shifted_df_spanning
del df_blast_gp1_short
del df_blast_gp1
del condition
# del df_blast_gp1
del df_blast_gp1_fd
del df_blast_gp1_rv
del gk_span_save
del id
del ab6
del extra_sp
del df_SAfile
del index_T1_dict
del index_T2_dict
del df_blast_T1
del df_blast_T2
del df_blast_uniq
del df_blast_f1
del df_blast_f2 
del df_blast_f1_list
del df_blast_f2_list
del df_blast_uniq_gp1_list

end_time = time.time()
usage = resource.getrusage(resource.RUSAGE_SELF)
max_memory = usage.ru_maxrss  # Memory in kilobytes

print(f"Execution Time: {end_time - start_time} seconds")
print(f"Max Memory Usage: {max_memory / 1024} MB") 