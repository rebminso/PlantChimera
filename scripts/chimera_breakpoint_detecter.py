#############################################################################################################################
#                                     Breakpoint Identification using Split Reads
#############################################################################################################################

import re
import glob
import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO
import warnings
warnings.filterwarnings("ignore", category=pd.errors.SettingWithCopyWarning)
import resource
import time

start_time = time.time()

parser = argparse.ArgumentParser(description='a list of seq from csv file')
parser.add_argument('input',help='previous_input')
parser.add_argument('overlap_query', help='set it to 50')
parser.add_argument('gap', help='set it to 5')
parser.add_argument('output',help='intfile')
parser.add_argument('AnchorLength')
parser.add_argument('Rlength')
parser.add_argument('subject_query', help='set it to 20')
parser.add_argument('gap_subject', help='set it 0')
args = parser.parse_args()


df_blast = pd.read_csv(args.input,delimiter='\t')
df_blast_1 = df_blast[df_blast['query_ID'].str.endswith('_1')]
df_blast_2 = df_blast[df_blast['query_ID'].str.endswith('_2')]
df_blast_1[['query_ID', 'frame_id1']] = df_blast_1['query_ID'].str.split('_', expand=True)
df_blast_2[['query_ID', 'frame_id1']] = df_blast_2['query_ID'].str.split('_', expand=True)

# Bring all the readpairs in a single line
full_df = df_blast_1.merge(df_blast_2, on='query_ID', how='outer', suffixes=('_fd', '_rv'))

# Define the functions
def assignStrand_fd(x):
    if (x['strand_symbol_1_fd'] == '+' and x['strand_symbol_1_rv'] == '-'):
        return '+'
    elif (x['strand_symbol_1_fd'] == '-' and x['strand_symbol_1_rv'] == '+'):
        return '-'
    elif (x['strand_symbol_1_fd'] == '+' and pd.isna(x['strand_symbol_1_rv'])):
        return '+'
    elif (x['strand_symbol_1_fd'] == '-' and pd.isna(x['strand_symbol_1_rv'])):
        return '-'
    elif (pd.isna(x['strand_symbol_1_fd']) and x['strand_symbol_1_rv'] == '-'):
        return '+'
    elif (pd.isna(x['strand_symbol_1_fd']) and x['strand_symbol_1_rv'] == '+'):
        return '-'
    else:
        return 'No value'

def assignStrand_rv(x):
    if (x['strand_symbol_2_fd'] == '+' and x['strand_symbol_2_rv'] == '-'):
        return '+'
    elif (x['strand_symbol_2_fd'] == '-' and x['strand_symbol_2_rv'] == '+'):
        return '-'
    elif (x['strand_symbol_2_fd'] == '+' and pd.isna(x['strand_symbol_2_rv'])):
        return '+'
    elif (x['strand_symbol_2_fd'] == '-' and pd.isna(x['strand_symbol_2_rv'])):
        return '-'
    elif (pd.isna(x['strand_symbol_2_fd']) and x['strand_symbol_2_rv'] == '-'):
        return '+'
    elif (pd.isna(x['strand_symbol_2_fd']) and x['strand_symbol_2_rv'] == '+'):
        return '-'
    else:
        return 'No value'

# Applying the functions to the DataFrame
full_df['final_transcript_strand_1'] = full_df.apply(assignStrand_fd, axis=1)
full_df['final_transcript_strand_2'] = full_df.apply(assignStrand_rv, axis=1)

st = full_df[['query_ID', 'Transcript_ID_T1_fd', 'query_start_end_T1_fd','subject_start_end_T1_fd', 'strand_T1_fd', 'expect_T1_fd',
       'identity_T1_fd', 'Transcript_ID_T2_fd', 'query_start_end_T2_fd',
       'subject_start_end_T2_fd', 'strand_T2_fd', 'expect_T2_fd',
       'identity_T2_fd', 'strand_symbol_1_fd', 'strand_symbol_2_fd',
       'final_transcript_strand_fd', 'frame_id1_fd', 'Transcript_ID_T1_rv',
       'query_start_end_T1_rv', 'subject_start_end_T1_rv', 'strand_T1_rv',
       'expect_T1_rv', 'identity_T1_rv', 'Transcript_ID_T2_rv',
       'query_start_end_T2_rv', 'subject_start_end_T2_rv', 'strand_T2_rv',
       'expect_T2_rv', 'identity_T2_rv', 'strand_symbol_1_rv',
       'strand_symbol_2_rv', 'final_transcript_strand_rv', 'frame_id1_rv',
       'final_transcript_strand_1', 'final_transcript_strand_2']]
#normal_splitRead
samstrand = st.loc[(full_df['final_transcript_strand_1'] =='-') & (full_df['final_transcript_strand_2'] =='-') |(full_df['final_transcript_strand_1'] =='+') & (full_df['final_transcript_strand_2'] =='+') ]
# normal_complex_Read
samstrand_box2 = st.loc[(full_df['final_transcript_strand_1'] == '+') & (full_df['final_transcript_strand_2'] == '-') |(full_df['final_transcript_strand_1'] == '-') & (full_df['final_transcript_strand_2'] == '+') ]
# final head and tail transcript strands of complex reads
cmpl = samstrand_box2[['query_ID','Transcript_ID_T1_fd','Transcript_ID_T2_fd','Transcript_ID_T1_rv','Transcript_ID_T2_rv', 'final_transcript_strand_1', 	'final_transcript_strand_2']]
# Create the new columns
cmpl['head_strand'] = cmpl['Transcript_ID_T1_fd'].fillna('') + ';' + cmpl['Transcript_ID_T1_rv'].fillna('')
cmpl['tail_strand'] = cmpl['Transcript_ID_T2_fd'].fillna('') + ';' + cmpl['Transcript_ID_T2_rv'].fillna('')
def process_row(row):
    if pd.isna(row):  # Check if the value is NaN
        return ''
    if not isinstance(row, str):  # Ensure the value is a string
        return row
    # Split the values, convert to set to keep only unique values, and join back with ';'
    unique_values = set(row.split(';'))
    return ';'.join(unique_values)

# Apply the function to the specific column
cmpl['headT'] = cmpl['head_strand'].apply(process_row)
cmpl['tailT'] = cmpl['tail_strand'].apply(process_row)

cmpl = cmpl[['query_ID','headT', 'tailT',  'final_transcript_strand_1', 'final_transcript_strand_2']]
cmpl['headT'] = cmpl['headT'].str.strip(';')
cmpl['tailT'] = cmpl['tailT'].str.strip(';')
cmpl_inv = cmpl.rename(columns={
    'headT': 'tailT',
    'tailT': 'headT',
    'final_transcript_strand_1': 'final_transcript_strand_2',
    'final_transcript_strand_2': 'final_transcript_strand_1'
})
cmpl_inv
cmpl_all = pd.concat([cmpl,cmpl_inv],axis=0)
cmpl_all
# final head and tail transcript strands of Normal split reads
Nrm = samstrand[['query_ID','Transcript_ID_T1_fd','Transcript_ID_T2_fd','Transcript_ID_T1_rv','Transcript_ID_T2_rv', 'final_transcript_strand_1', 	'final_transcript_strand_2']]
# Create the new columns
Nrm['head_strand'] = Nrm['Transcript_ID_T1_fd'].fillna('') + ';' + Nrm['Transcript_ID_T1_rv'].fillna('')
Nrm['tail_strand'] = Nrm['Transcript_ID_T2_fd'].fillna('') + ';' + Nrm['Transcript_ID_T2_rv'].fillna('')
def process_row(row):
    if pd.isna(row):  # Check if the value is NaN
        return ''
    if not isinstance(row, str):  # Ensure the value is a string
        return row
    # Split the values, convert to set to keep only unique values, and join back with ';'
    unique_values = set(row.split(';'))
    return ';'.join(unique_values)
# Apply the function to the specific column
Nrm['headT'] = Nrm['head_strand'].apply(process_row)
Nrm['tailT'] = Nrm['tail_strand'].apply(process_row)
Nrm = Nrm[['query_ID','headT', 'tailT',  'final_transcript_strand_1', 'final_transcript_strand_2']]
Nrm['headT'] = Nrm['headT'].str.strip(';')
Nrm['tailT'] = Nrm['tailT'].str.strip(';')
Nrm_inv = Nrm.rename(columns={
    'headT': 'tailT',
    'tailT': 'headT',
    'final_transcript_strand_1': 'final_transcript_strand_2',
    'final_transcript_strand_2': 'final_transcript_strand_1'
})
Nrm_all = pd.concat([Nrm,Nrm_inv],axis=0)
# new variable just to keep the number of columns same
normalsplitRead =  full_df.loc[(full_df['final_transcript_strand_1'] =='-') & (full_df['final_transcript_strand_2'] =='-') |(full_df['final_transcript_strand_1'] =='+') & (full_df['final_transcript_strand_2'] =='+') ]
normalRead_strands = normalsplitRead[['query_ID','final_transcript_strand_1','final_transcript_strand_2']]
# these are the complex reads, reads with -+,+-
samstrand_box2.to_csv('Complex_reads.csv',index=False,sep='\t')
samstrand_box2.drop(['final_transcript_strand_1', 'final_transcript_strand_2'],axis=1,inplace=True)
# Now I am trying to mimic what been done for normal reads
samstrand_box2[['query_start_T1_fd', 'query_end_T1_fd']] = samstrand_box2['query_start_end_T1_fd'].str.split('-', expand=True)
samstrand_box2[['query_start_T2_fd', 'query_end_T2_fd']] = samstrand_box2['query_start_end_T2_fd'].str.split('-', expand=True)
samstrand_box2[['query_start_T1_rv', 'query_end_T1_rv']] = samstrand_box2['query_start_end_T1_rv'].str.split('-', expand=True)
samstrand_box2[['query_start_T2_rv', 'query_end_T2_rv']] = samstrand_box2['query_start_end_T2_rv'].str.split('-', expand=True)
samstrand_box2[['subject_start_T1_fd', 'subject_end_T1_fd']] = samstrand_box2['subject_start_end_T1_fd'].str.split('-', expand=True)
samstrand_box2[['subject_start_T2_fd', 'subject_end_T2_fd']] = samstrand_box2['subject_start_end_T2_fd'].str.split('-', expand=True)
samstrand_box2[['subject_start_T1_rv', 'subject_end_T1_rv']] = samstrand_box2['subject_start_end_T1_rv'].str.split('-', expand=True)
samstrand_box2[['subject_start_T2_rv', 'subject_end_T2_rv']] = samstrand_box2['subject_start_end_T2_rv'].str.split('-', expand=True)
for i in ['subject_start_T1_fd', 'subject_end_T1_fd', 'subject_start_T2_fd', 'subject_end_T2_fd', 'subject_start_T1_rv', 'subject_end_T1_rv', 'subject_start_T2_rv', 'subject_end_T2_rv']:
    samstrand_box2[i] = pd.to_numeric(samstrand_box2[i])
for i in ['query_start_T1_fd', 'query_end_T1_fd', 'query_start_T2_fd', 'query_end_T2_fd', 'query_start_T1_rv', 'query_end_T1_rv', 'query_start_T2_rv', 'query_end_T2_rv']:
    samstrand_box2[i] = pd.to_numeric(samstrand_box2[i])
samstrand_box2['subject_start_T1_fd'], samstrand_box2['subject_end_T1_fd'] = np.where(samstrand_box2['subject_start_T1_fd'] > samstrand_box2['subject_end_T1_fd'], (samstrand_box2['subject_end_T1_fd'], samstrand_box2['subject_start_T1_fd']), (samstrand_box2['subject_start_T1_fd'], samstrand_box2['subject_end_T1_fd']))
samstrand_box2['subject_start_T2_fd'], samstrand_box2['subject_end_T2_fd'] = np.where(samstrand_box2['subject_start_T2_fd'] > samstrand_box2['subject_end_T2_fd'], (samstrand_box2['subject_end_T2_fd'], samstrand_box2['subject_start_T2_fd']), (samstrand_box2['subject_start_T2_fd'], samstrand_box2['subject_end_T2_fd']))
samstrand_box2['subject_start_T1_rv'], samstrand_box2['subject_end_T1_rv'] = np.where(samstrand_box2['subject_start_T1_rv'] > samstrand_box2['subject_end_T1_rv'], (samstrand_box2['subject_end_T1_rv'], samstrand_box2['subject_start_T1_rv']), (samstrand_box2['subject_start_T1_rv'], samstrand_box2['subject_end_T1_rv']))
samstrand_box2['subject_start_T2_rv'], samstrand_box2['subject_end_T2_rv'] = np.where(samstrand_box2['subject_start_T2_rv'] > samstrand_box2['subject_end_T2_rv'], (samstrand_box2['subject_end_T2_rv'], samstrand_box2['subject_start_T2_rv']), (samstrand_box2['subject_start_T2_rv'], samstrand_box2['subject_end_T2_rv']))
def find_overlap_T1(row):
    range1_start, range1_end, range2_start, range2_end = row['subject_start_T1_fd'], row['subject_end_T1_fd'], row['subject_start_T1_rv'], row['subject_end_T1_rv']
    overlap_start = max(range1_start, range2_start)
    overlap_end = min(range1_end, range2_end)
    if np.isnan(row['subject_start_T1_fd']) or np.isnan(row['subject_start_T1_rv']) or np.isnan(row['subject_end_T1_fd']) or np.isnan(row['subject_end_T1_rv']):
        return np.nan
    return overlap_end - overlap_start +  1 if overlap_start <= overlap_end else  0
samstrand_box2['overlap_subject_T1'] = samstrand_box2.apply(find_overlap_T1, axis=1)
samstrand_box2.to_csv('df_complexReads.txt',sep='\t',index=False)
samstrand_box2 = pd.read_csv('df_complexReads.txt',delimiter='\t')
samstrand_box2[['Transcript_ID_T1_fd', 'Transcript_ID_T2_fd', 'Transcript_ID_T1_rv', 'Transcript_ID_T2_rv']].tail()
df = samstrand_box2

# Function to compare values, return the formatted string, and the difference
def compare_values(row, col1, col2):
    try:
        # Check if cells are empty
        if pd.isna(row[col1]) or pd.isna(row[col2]):
            return "NAN", 0, 0

        a, b = map(int, row[col1].split('-'))
        c, d = map(int, row[col2].split('-'))

        if b >= c:
            if a >= c and b >= d:
                if a < d:
                    formatted_value = "{}-{}".format(a, d)
                    difference = abs(a - d) + 1
                    gap = 0
                else:
                    #formatted_value = "no overlap"
                    if a == d:
                          formatted_value= "{}-{}".format(a,d)
                          difference = abs(a-d)+1
                          gap = 0
                    else:
                          formatted_value = "no overlap"
                          difference = 0
                          gap = abs(a - d) - 1

            elif a >= c and b <= d:
                formatted_value = "{}-{}".format(a, b)
                difference = abs(a - b) + 1
                gap = 0
            elif a <= c and b >= d:
                formatted_value = "{}-{}".format(c, d)
                difference = abs(c - d) + 1
                gap = 0
            elif a <= c and b <= d:
                formatted_value = "{}-{}".format(c, b)
                difference = abs(c - b) + 1
                gap = 0
        elif b < c:
            formatted_value = "no overlap"
            difference = 0
            gap = abs(b - c) - 1

        return formatted_value, difference, gap
    except KeyError as e:
        print("Column not found: {}".format(e))
        return None, None, None

# Apply the function to each row of the DataFrame for fd_overlap and rv_overlap
df[['fd_overlap', 'overlap_query_fd','gap_query_fd']] = df.apply(lambda row: pd.Series(compare_values(row, 'query_start_end_T1_fd', 'query_start_end_T2_fd')), axis=1)
df[['rv_overlap', 'overlap_query_rv','gap_query_rv']] = df.apply(lambda row: pd.Series(compare_values(row, 'query_start_end_T1_rv', 'query_start_end_T2_rv')), axis=1)
df.drop([ 'fd_overlap', 'rv_overlap', ],axis = 1, inplace=True)
samstrand_box2 = df
del df
def find_overlap_T2(row):
    range1_start, range1_end, range2_start, range2_end = row['subject_start_T2_fd'], row['subject_end_T2_fd'], row['subject_start_T2_rv'], row['subject_end_T2_rv']
    overlap_start = max(range1_start, range2_start)
    overlap_end = min(range1_end, range2_end)
    if np.isnan(row['subject_start_T2_fd']) or np.isnan(row['subject_start_T2_rv']) or np.isnan(row['subject_end_T2_fd']) or np.isnan(row['subject_end_T2_rv']):
        return np.nan
    return overlap_end - overlap_start +  1 if overlap_start <= overlap_end else  0
# def find_overlap_query_fd(row):
#     range1_start, range1_end, range2_start, range2_end = row['query_start_T1_fd'], row['query_end_T1_fd'], row['query_start_T2_fd'], row['query_end_T2_fd']
#     overlap_start = max(range1_start, range2_start)
#     overlap_end = min(range1_end, range2_end)
#     if np.isnan(row['query_start_T1_fd']) or np.isnan(row['query_start_T2_fd']) or np.isnan(row['query_end_T1_fd']) or np.isnan(row['query_end_T2_fd']):
#         return np.nan
#     return overlap_end - overlap_start +  1 if overlap_start <= overlap_end else  0
# def find_overlap_query_rv(row):
#     range1_start, range1_end, range2_start, range2_end = row['query_start_T1_fd'], row['query_end_T1_rv'], row['query_start_T2_fd'], row['query_end_T2_rv']
#     overlap_start = max(range1_start, range2_start)
#     overlap_end = min(range1_end, range2_end)
#     if np.isnan(row['query_start_T1_fd']) or np.isnan(row['query_start_T2_fd']) or np.isnan(row['query_end_T1_rv']) or np.isnan(row['query_end_T2_rv']):
#         return np.nan
#     return overlap_end - overlap_start +  1 if overlap_start <= overlap_end else  0
samstrand_box2['overlap_subject_T2'] = samstrand_box2.apply(find_overlap_T2, axis=1)
# samstrand_box2['overlap_query_fd_md'] = samstrand_box2.apply(find_overlap_query_fd, axis=1)
# samstrand_box2['overlap_query_rv_md'] = samstrand_box2.apply(find_overlap_query_rv, axis=1)
# here namrat's code will be added that will include overlap_query_fd and overlap_query_rv
# samstrand_box2.head(4)

# samstrand_box2.iloc[:20,-11:]



# samstrand_box2.head(20)

# samstrand_box2.iloc[:20,-11:]

# samstrand_box2['gap_query_fd_md'] = samstrand_box2.apply(lambda x: abs(min(x['query_end_T1_fd'], x['query_end_T2_fd']) - max(x['query_start_T1_fd'], x['query_start_T2_fd'])) if x['overlap_query_fd'] == 0 else 0, axis=1)
# samstrand_box2['gap_query_rv_md'] = samstrand_box2.apply(lambda x: abs(min(x['query_end_T1_rv'], x['query_end_T2_rv']) - max(x['query_start_T1_rv'], x['query_start_T2_rv'])) if x['overlap_query_rv'] == 0 else 0, axis=1)
# samstrand_box2.iloc[:20,-11:]

# samstrand_box2.iloc[:20,-11:]

# replace NAN with 0
# samstrand_box2['gap_query_fd_md'] = samstrand_box2['gap_query_fd_md'].fillna(0)
# samstrand_box2['gap_query_rv_md'] = samstrand_box2['gap_query_rv_md'].fillna(0)

# # subtracting 1 from query if the vlaue if more than 0
# # samstrand_box2['overlap_query_fd'] =
# samstrand_box2['gap_query_fd_md'] = samstrand_box2['gap_query_fd_md'].apply(lambda x: x - 1 if x > 0 else 0)
# samstrand_box2['gap_query_rv_md'] = samstrand_box2['gap_query_rv_md'].apply(lambda x: x - 1 if x > 0 else 0)
# samstrand_box2.iloc[:20,-11:]

# po = samstrand_box2.loc[samstrand_box2['query_ID'].str.contains('SRR11404272.10007456') & (samstrand_box2['gap_query_fd'] == 18 )]
# po.to_csv('gapquery.csv',index=False,sep='\t')

# samstrand_box2['gap_query_fd'] = samstrand_box2.apply(lambda x: abs(min(x['query_end_T1_fd'], x['query_end_T2_fd']) - max(x['query_start_T1_fd'], x['query_start_T2_fd'])) if x['overlap_query_fd'] == 0 else 0, axis=1)
# samstrand_box2['gap_query_rv'] = samstrand_box2.apply(lambda x: abs(min(x['query_end_T1_rv'], x['query_end_T2_rv']) - max(x['query_start_T1_rv'], x['query_start_T2_rv'])) if x['overlap_query_rv'] == 0 else 0, axis=1)
samstrand_box2['gap_subject_T1'] = samstrand_box2.apply(lambda x: abs(min(x['subject_end_T1_fd'], x['subject_end_T1_rv']) - max(x['subject_start_T1_fd'], x['subject_start_T1_rv'])) if x['overlap_subject_T1'] == 0 else 0, axis=1)
samstrand_box2['gap_subject_T2'] = samstrand_box2.apply(lambda x: abs(min(x['subject_end_T2_fd'], x['subject_end_T2_rv']) - max(x['subject_start_T2_fd'], x['subject_start_T2_rv'])) if x['overlap_subject_T2'] == 0 else 0, axis=1)
samstrand_box2.iloc[:, -9:] = samstrand_box2.iloc[:, -9:].fillna(0)
##_____________________________________________________________________PARAMETER__________________OVERLAPGAP__________________________________________________________ (4% of 150)
samstrand_box2['gap_condition'] = samstrand_box2.apply(lambda x: 'NO' if (x['gap_query_fd'] <= int(args.gap))  and x['gap_query_rv'] <= int(args.gap)  else 'YES', axis=1)
# samstrand_box2.iloc[:,-9:]
##_____________________________________________________________________PARAMETER___________________OVERLAPQUERY_________________________________________________________ (30% of 150)
samstrand_box2['overlap_query_condition'] = samstrand_box2.apply(lambda x: 'NO' if (x['overlap_query_fd'] <=  int(args.overlap_query) and x['overlap_query_rv'] <= int(args.overlap_query)) else 'YES', axis=1)
# samstrand_box2.iloc[:,-11:]
# samstrand_box2.apply(lambda x: 'NO' if (x['overlap_subject_T1'] == 0 and x['overlap_subject_T2'] == 0) else \
                                                        #    'NO' if ((x['overlap_subject_T1'] > 0 or x['overlap_subject_T2'] > 0) and (x['gap_subject_T1'] <= 5 and x['gap_subject_T2'] <= 5)) else 'YES', axis=1)
# If both overlap_subject_T1 and overlap_subject_T2 are 0, and both gap_subject_T1 and gap_subject_T2 are less than or equal to 5, return 'NO'.
# If either overlap_subject_T1 or overlap_subject_T2 is greater than 0, return 'NO'.
##_____________________________________________________________________PARAMETER________________OVERLAPSUBJECT/ ANCHOR LENGTH____________________________________________________________ (20% of 150)
samstrand_box2['overlap_subject_condition'] = samstrand_box2.apply(lambda row: 'NO' if row['overlap_subject_T1'] == int(args.subject_query) and row['overlap_subject_T2'] == int(args.subject_query) and row['gap_subject_T1'] <=  int(args.gap_subject) and row['gap_subject_T2'] <= int(args.gap_subject) else \
    'NO' if row['overlap_subject_T1'] > int(args.subject_query) or row['overlap_subject_T2'] > int(args.subject_query) else 'YES', axis=1)
# samstrand_box2.iloc[:,-9:]
# samstrand_box2['overlap_subject_condition'] = samstrand_box2.apply(lambda x: 'NO' if (x['overlap_subject_T1'] == 0 and x['overlap_subject_T2'] == 0) else \
#                                                            'NO' if ((x['overlap_subject_T1'] > 0 or x['overlap_subject_T2'] > 0) and (x['gap_subject_T1'] <= 5 and x['gap_subject_T2'] <= 5)) else 'YES', axis=1)
# samstrand_box2.head(4)
# ##CHECKKK!!!!!!
# samstrand_box2.columns

# # last_10_columns = samstrand_box2.iloc[:, -10:]
# checkup_tab = samstrand_box2[[ 'query_start_end_T1_fd', 'subject_start_end_T1_fd',  'query_start_end_T2_fd', 'subject_start_end_T2_fd', 'query_start_end_T1_rv', 'query_start_end_T2_rv', 'subject_start_end_T2_rv', 'subject_start_end_T1_rv',  'overlap_subject_T1', 'overlap_subject_T2', 'overlap_query_fd',
#        'overlap_query_rv', 'gap_query_fd', 'gap_query_rv',
#        'gap_subject_T1', 'gap_subject_T2', 'gap_condition','overlap_query_condition',
#        'overlap_subject_condition']]
# checkup_tab
# # checkup_tab.to_csv('complex_overlaps.txt',sep='\t',index=False)
# checkup_tab.to_csv('complex_overlaps_2.txt',sep='\t',index=False)

samstrand_box2['reads_skipped'] = samstrand_box2.apply(lambda x: 'YES' if (x['overlap_query_condition'] == 'YES' or x['overlap_subject_condition'] == 'YES' or x['gap_condition'] == 'YES') else 'NO', axis=1)

df_blast_skipped_3 = samstrand_box2[samstrand_box2['reads_skipped'] == 'YES']

samstrand_box2.drop(columns=[col for col in samstrand_box2.columns if 'skipped' in col], inplace=True)
df_blast_skipped_3.drop(columns=[col for col in df_blast_skipped_3.columns if 'skipped' in col], inplace=True)

samstrand_box2 = samstrand_box2[~samstrand_box2.isin(df_blast_skipped_3.to_dict(orient='list')).all(1)]
if samstrand_box2.shape[0] > 0:
    samstrand_box2['query_start_T1_fd'] = pd.to_numeric(samstrand_box2['query_start_T1_fd'])
    samstrand_box2['query_start_T2_fd'] = pd.to_numeric(samstrand_box2['query_start_T2_fd'])

    samstrand_box2['query_start_T1_rv'] = pd.to_numeric(samstrand_box2['query_start_T1_rv'])
    samstrand_box2['query_start_T2_rv'] = pd.to_numeric(samstrand_box2['query_start_T2_rv'])
    samstrand_box2.head(4)
    samstrand_box2[['Transcript_ID_T1_fd','query_start_T1_fd', 'strand_symbol_1_fd', 'Transcript_ID_T2_fd' ,'query_start_T2_fd', 'strand_symbol_2_fd']]
    samstrand_box2['HEAD_plus_fd'] = samstrand_box2.apply(lambda x: x['Transcript_ID_T1_fd'] if (x['strand_symbol_1_fd'] == '+' and x['strand_symbol_2_fd'] == '+' \
                                                                                and x['query_start_T1_fd'] < x['query_start_T2_fd']) else x['Transcript_ID_T2_fd'] if (x['strand_symbol_1_fd'] == '+' and x['strand_symbol_2_fd'] == '+' \
                                                                                    and x['query_start_T1_fd'] > x['query_start_T2_fd'])  else "NONE", axis=1)

    samstrand_box2['HEAD_plus_rv'] = samstrand_box2.apply(lambda x: x['Transcript_ID_T1_rv'] if (x['strand_symbol_1_rv'] == '+' and x['strand_symbol_2_rv'] == '+' \
                                                                                    and x['query_start_T1_rv'] < x['query_start_T2_rv']) else x['Transcript_ID_T2_rv'] if (x['strand_symbol_1_rv'] == '+' and x['strand_symbol_2_rv'] == '+' \
                                                                                    and x['query_start_T1_rv'] > x['query_start_T2_rv']) else "NONE", axis=1)
    samstrand_box2['HEAD_minus_fd'] = samstrand_box2.apply(lambda x: x['Transcript_ID_T1_fd'] if (x['strand_symbol_1_fd'] == '-' and x['strand_symbol_2_fd'] == '-' \
                                                                                    and x['query_start_T1_fd'] > x['query_start_T2_fd']) else x['Transcript_ID_T2_fd'] if (x['strand_symbol_1_fd'] == '-' and x['strand_symbol_2_fd'] == '-' \
                                                                                    and x['query_start_T1_fd'] < x['query_start_T2_fd']) else "NONE", axis=1)
    samstrand_box2['HEAD_minus_rv'] = samstrand_box2.apply(lambda x: x['Transcript_ID_T1_rv'] if (x['strand_symbol_1_rv'] == '-' and x['strand_symbol_2_rv'] == '-' \
                                                                                    and x['query_start_T1_rv'] > x['query_start_T2_rv']) else x['Transcript_ID_T2_rv'] if (x['strand_symbol_1_rv'] == '-' and x['strand_symbol_2_rv'] == '-' \
                                                                                    and x['query_start_T1_rv'] < x['query_start_T2_rv']) else "NONE", axis=1)
    samstrand_box2["Head"]  = samstrand_box2[['HEAD_plus_fd', 'HEAD_plus_rv', 'HEAD_minus_fd', 'HEAD_minus_rv']].apply(lambda x: "".join(list(filter(lambda y: y != 'NONE', x.unique()))), axis=1)
    samstrand_box2.iloc[:,-5:]
    samstrand_box2.shape

    samstrand_box2[['Transcript_ID_T1_fd', 'Transcript_ID_T2_fd', 'Transcript_ID_T1_rv', 'Transcript_ID_T2_rv']] = samstrand_box2[['Transcript_ID_T1_fd', 'Transcript_ID_T2_fd', 'Transcript_ID_T1_rv', 'Transcript_ID_T2_rv']].astype(str)
    samstrand_box2["Tail"]  = samstrand_box2[['Transcript_ID_T1_fd', 'Transcript_ID_T2_fd', 'Transcript_ID_T1_rv', 'Transcript_ID_T2_rv']].apply(lambda x: list(filter(lambda y: y != 'nan', x.unique())), axis=1)

    samstrand_box2['Tail'] = samstrand_box2.apply(lambda x: "".join(list(filter(lambda y: y not in x['Head'], x['Tail']))), axis=1)
    # samstrand_box2.head(3)

    # samstrand_box2.iloc[:,-11:]

    samstrand_box2.drop(columns=[col for col in samstrand_box2.columns if 'HEAD_' in col], inplace=True)

    samstrand_box2_NaN = samstrand_box2[samstrand_box2.isna().any(axis=1)]
    samstrand_box2_NoNaN = samstrand_box2[~samstrand_box2.isna().any(axis=1)]
    samstrand_box2_NoNaN_subset = samstrand_box2_NoNaN.copy()

    samstrand_box2_NoNaN_subset['keep'] = samstrand_box2_NoNaN_subset.apply(lambda x: 'reverse' if x['overlap_query_fd'] == 0 and x['overlap_query_rv'] > 0 \
                                                                    else 'forward' if x['overlap_query_fd'] > 0 and x['overlap_query_rv'] == 0 \
                                                                    else 'reverse' if x['overlap_query_fd'] > 0 and x['overlap_query_rv'] > 0 and x['overlap_query_fd'] > x['overlap_query_rv'] \
                                                                    else 'forward' if x['overlap_query_fd'] > 0 and x['overlap_query_rv'] > 0 and x['overlap_query_fd'] < x['overlap_query_rv'] \
                                                                    else 'reverse' if x['gap_query_fd'] > 0 and x['gap_query_rv'] > 0 and x['gap_query_fd'] > x['gap_query_rv'] \
                                                                    else 'forward' if x['gap_query_fd'] > 0 and x['gap_query_rv'] > 0 and x['gap_query_fd'] < x['gap_query_rv'] \
                                                                    else 'Either' if x['overlap_query_fd'] == 0 and x['overlap_query_rv'] == 0 and x['gap_query_fd'] == 0 and x['gap_query_rv'] == 0 \
                                                                    else 'reverse' if x['overlap_query_fd'] == x['overlap_query_rv'] and x['gap_query_fd'] == 0 and x['gap_query_rv'] > 0 \
                                                                    else 'forward' if x['overlap_query_fd'] == x['overlap_query_rv'] and x['gap_query_fd'] > 0 and x['gap_query_rv'] == 0 \
                                                                    else "Either", axis=1)
    # samstrand_box2.head(3)
    samstrand_box2_NaN['keep'] = samstrand_box2_NaN.apply(lambda x: 'forward' if (x['Transcript_ID_T1_fd'] != 'nan' and x['Transcript_ID_T2_fd'] != 'nan') \
                                                            else 'reverse' if (x['Transcript_ID_T1_rv'] != 'nan' and x['Transcript_ID_T2_rv'] != 'nan') else 'Either', axis=1)

    samstrand_box2 = pd.concat([samstrand_box2_NoNaN_subset, samstrand_box2_NaN], axis=0)
    samstrand_box2.drop(columns=[col for col in samstrand_box2.columns if 'condition' in col], inplace=True)
    samstrand_box2.drop(columns=['keep'], inplace=True)
    # samstrand_box2

    # samstrand_box2.columns

    samstrand_box2['subject_start_HEAD_fd'] = samstrand_box2.apply(lambda x: x['subject_start_T1_fd'] if x['Transcript_ID_T1_fd'] == x['Head'] \
                                                        else x['subject_start_T2_fd'] if x['Transcript_ID_T2_fd'] == x['Head'] else "NONE", axis=1)
    samstrand_box2['subject_end_HEAD_fd'] = samstrand_box2.apply(lambda x: x['subject_end_T1_fd'] if x['Transcript_ID_T1_fd'] == x['Head'] \
                                                            else x['subject_end_T2_fd'] if x['Transcript_ID_T2_fd'] == x['Head'] else "NONE", axis=1)
    samstrand_box2['subject_start_HEAD_rv'] = samstrand_box2.apply(lambda x: x['subject_start_T1_rv'] if x['Transcript_ID_T1_rv'] == x['Head'] \
                                                            else x['subject_start_T2_rv'] if x['Transcript_ID_T2_rv'] == x['Head'] else "NONE", axis=1)
    samstrand_box2['subject_end_HEAD_rv'] = samstrand_box2.apply(lambda x: x['subject_end_T1_rv'] if x['Transcript_ID_T1_rv'] == x['Head'] \
                                                                else x['subject_end_T2_rv'] if x['Transcript_ID_T2_rv'] == x['Head'] else "NONE", axis=1)

    samstrand_box2['subject_start_TAIL_fd'] = samstrand_box2.apply(lambda x: x['subject_start_T1_fd'] if x['Transcript_ID_T1_fd'] == x['Tail'] \
                                                            else x['subject_start_T2_fd'] if x['Transcript_ID_T2_fd'] == x['Tail'] else "NONE", axis=1)
    samstrand_box2['subject_end_TAIL_fd'] = samstrand_box2.apply(lambda x: x['subject_end_T1_fd'] if x['Transcript_ID_T1_fd'] == x['Tail'] \
                                                            else x['subject_end_T2_fd'] if x['Transcript_ID_T2_fd'] == x['Tail'] else "NONE", axis=1)
    samstrand_box2['subject_start_TAIL_rv'] = samstrand_box2.apply(lambda x: x['subject_start_T1_rv'] if x['Transcript_ID_T1_rv'] == x['Tail'] \
                                                                else x['subject_start_T2_rv'] if x['Transcript_ID_T2_rv'] == x['Tail'] else "NONE", axis=1)
    samstrand_box2['subject_end_TAIL_rv'] = samstrand_box2.apply(lambda x: x['subject_end_T1_rv'] if x['Transcript_ID_T1_rv'] == x['Tail'] \
                                else x['subject_end_T2_rv'] if x['Transcript_ID_T2_rv'] == x['Tail'] else "NONE", axis=1)

    # samstrand_box2.head(3)
    samstrand_box2['strand_HEAD_fd'] = samstrand_box2.apply(lambda x: x['strand_symbol_1_fd'] if x['Transcript_ID_T1_fd'] == x['Head'] \
                                                        else x['strand_symbol_2_fd'] if x['Transcript_ID_T2_fd'] == x['Head'] else "NONE", axis=1)
    samstrand_box2['strand_HEAD_rv'] = samstrand_box2.apply(lambda x: x['strand_symbol_1_rv'] if x['Transcript_ID_T1_rv'] == x['Head'] \
                                                        else x['strand_symbol_2_rv'] if x['Transcript_ID_T2_rv'] == x['Head'] else "NONE", axis=1)
    samstrand_box2['strand_TAIL_fd'] = samstrand_box2.apply(lambda x: x['strand_symbol_1_fd'] if x['Transcript_ID_T1_fd'] == x['Tail'] \
                                                        else x['strand_symbol_2_fd'] if x['Transcript_ID_T2_fd'] == x['Tail'] else "NONE", axis=1)
    samstrand_box2['strand_TAIL_rv'] = samstrand_box2.apply(lambda x: x['strand_symbol_1_rv'] if x['Transcript_ID_T1_rv'] == x['Tail'] \
                                                        else x['strand_symbol_2_rv'] if x['Transcript_ID_T2_rv'] == x['Tail'] else "NONE", axis=1)
    # samstrand_box2.head(3)
    samstrand_box2['retain_HEAD'] = samstrand_box2.apply(lambda x: 'YES' if (x['strand_HEAD_fd'] != x['strand_HEAD_rv'] and (x['subject_start_HEAD_fd'] == x['subject_start_HEAD_rv'] or x['subject_end_HEAD_fd'] == x['subject_end_HEAD_rv'])) else 'NO', axis=1)
    samstrand_box2['retain_TAIL'] = samstrand_box2.apply(lambda x: 'YES' if (x['strand_TAIL_fd'] != x['strand_TAIL_rv'] and (x['subject_start_TAIL_fd'] == x['subject_start_TAIL_rv'] or x['subject_end_TAIL_fd'] == x['subject_end_TAIL_rv'])) else 'NO', axis=1)

    samstrand_box2['retain'] = samstrand_box2.apply(lambda x: 'HEAD' if x['retain_HEAD'] == 'YES' and x['retain_TAIL'] == 'NO' \
                                            else 'TAIL' if x['retain_TAIL'] == 'YES' and x['retain_HEAD'] == 'NO' \
                                            else 'BOTH' if x['retain_HEAD'] == 'YES' and x['retain_TAIL'] == 'YES' else 'NONE', axis=1)

    samstrand_box2['keep'] = samstrand_box2.apply(lambda x: 'reverse' if x['overlap_query_fd'] == 0 and x['overlap_query_rv'] > 0 \
                                                else 'forward' if x['overlap_query_fd'] > 0 and x['overlap_query_rv'] == 0 \
                                                else 'reverse' if x['overlap_query_fd'] > 0 and x['overlap_query_rv'] > 0 and x['overlap_query_fd'] > x['overlap_query_rv'] \
                                                else 'forward' if x['overlap_query_fd'] > 0 and x['overlap_query_rv'] > 0 and x['overlap_query_fd'] < x['overlap_query_rv'] \
                                                else 'reverse' if x['gap_query_fd'] > 0 and x['gap_query_rv'] > 0 and x['gap_query_fd'] > x['gap_query_rv'] \
                                                else 'forward' if x['gap_query_fd'] > 0 and x['gap_query_rv'] > 0 and x['gap_query_fd'] < x['gap_query_rv'] \
                                                else 'Either' if x['overlap_query_fd'] == 0 and x['overlap_query_rv'] == 0 and x['gap_query_fd'] == 0 and x['gap_query_rv'] == 0 \
                                                else 'reverse' if x['overlap_query_fd'] == x['overlap_query_rv'] and x['gap_query_fd'] == 0 and x['gap_query_rv'] > 0 \
                                                else 'forward' if x['overlap_query_fd'] == x['overlap_query_rv'] and x['gap_query_fd'] > 0 and x['gap_query_rv'] == 0 \
                                                else "Either", axis=1)

    # samstrand_box2.head(3)


    samstrand_box2.loc[samstrand_box2['keep'] == 'reverse', samstrand_box2.columns[samstrand_box2.columns.str.endswith('_fd')]] = np.nan
    samstrand_box2.loc[samstrand_box2['keep'] == 'forward', samstrand_box2.columns[samstrand_box2.columns.str.endswith('_rv')]] = np.nan
    samstrand_box2.loc[samstrand_box2['keep'] == 'Either', samstrand_box2.columns[samstrand_box2.columns.str.endswith('_rv')]] = np.nan

    samstrand_box2.drop(columns=['keep'], inplace=True)

    samstrand_box2.drop(columns=[col for col in samstrand_box2.columns if 'HEAD' in col or 'TAIL' in col], inplace=True)
    # samstrand_box2.head(3)


    samstrand_box2_fd = samstrand_box2[['query_ID'] + [col for col in samstrand_box2.columns if col.endswith('_fd')] + ['Head', 'Tail', 'retain']]
    samstrand_box2_rv = samstrand_box2[['query_ID'] + [col for col in samstrand_box2.columns if col.endswith('_rv')] + ['Head', 'Tail', 'retain']]

    samstrand_box2_fd.dropna(inplace=True)
    samstrand_box2_rv.dropna(inplace=True)
    # samstrand_box2.head(3)

    samstrand_box2_fd.columns = samstrand_box2_fd.columns.str.replace('_fd', '')
    samstrand_box2_rv.columns = samstrand_box2_rv.columns.str.replace('_rv', '')
    samstrand_box2 = pd.concat([samstrand_box2_fd, samstrand_box2_rv], axis=0)
    samstrand_box2.rename(columns={'query_ID': 'ID_query'}, inplace=True)
    samstrand_box2.drop(columns=[col for col in samstrand_box2.columns if col.startswith('query_') or col.startswith('expect_') or col.startswith('identity_') or col.startswith('strand_symbol_') or col.startswith('frame_id') or col.startswith('subject_start_T1') or col.startswith('subject_end_T1') or col.startswith('subject_start_T2') or col.startswith('subject_end_T2')], inplace=True)
    # samstrand_box2.head(3)
    samstrand_box2_flip = samstrand_box2[~samstrand_box2['Transcript_ID_T1'].isin(samstrand_box2['Head'])]
    samstrand_box2 = samstrand_box2[samstrand_box2['Transcript_ID_T1'].isin(samstrand_box2['Head'])]
    # samstrand_box2.head(3)
    #  This line copies the values of columns indexed 1 to 3 (second to fourth columns) from samstrand_box2_flip into a new DataFrame called temp_vals.
    temp_vals = samstrand_box2_flip.iloc[:, 1:4].copy()
    #
    samstrand_box2_flip.iloc[:, 1:4] = samstrand_box2_flip.iloc[:, 4:7].values
    samstrand_box2_flip.iloc[:, 4:7] = temp_vals.values
    samstrand_box2 = pd.concat([samstrand_box2, samstrand_box2_flip], axis=0)

    samstrand_box2.drop(columns=['Head', 'Tail'], inplace=True)
    samstrand_box2[['q_ID_count', 'ID_query']] = samstrand_box2['ID_query'].str.split('#', expand=True)
    # samstrand_box2.head(3)

    samstrand_box2.drop(columns=['q_ID_count'], inplace=True)
    samstrand_box2.rename(columns={'ID_query': 'query_ID'}, inplace=True)
    samstrand_box2.sort_values(['query_ID'], inplace=True)

    # # samstrand_box2.to_csv(dir_path + "splitReads.csv", index=False)
    # samstrand_box2.head(3)

    # samstrand_box2.head()

    df_blast_1 =samstrand.copy()
    # df_blast_1.shape


    df_blast_1[['query_start_T1_fd', 'query_end_T1_fd']] = df_blast_1['query_start_end_T1_fd'].str.split('-', expand=True)
    df_blast_1[['query_start_T2_fd', 'query_end_T2_fd']] = df_blast_1['query_start_end_T2_fd'].str.split('-', expand=True)
    df_blast_1[['query_start_T1_rv', 'query_end_T1_rv']] = df_blast_1['query_start_end_T1_rv'].str.split('-', expand=True)
    df_blast_1[['query_start_T2_rv', 'query_end_T2_rv']] = df_blast_1['query_start_end_T2_rv'].str.split('-', expand=True)
    # df_blast_1.head(4)

    df_blast_1[['subject_start_T1_fd', 'subject_end_T1_fd']] = df_blast_1['subject_start_end_T1_fd'].str.split('-', expand=True)
    df_blast_1[['subject_start_T2_fd', 'subject_end_T2_fd']] = df_blast_1['subject_start_end_T2_fd'].str.split('-', expand=True)
    df_blast_1[['subject_start_T1_rv', 'subject_end_T1_rv']] = df_blast_1['subject_start_end_T1_rv'].str.split('-', expand=True)
    df_blast_1[['subject_start_T2_rv', 'subject_end_T2_rv']] = df_blast_1['subject_start_end_T2_rv'].str.split('-', expand=True)

    for i in ['subject_start_T1_fd', 'subject_end_T1_fd', 'subject_start_T2_fd', 'subject_end_T2_fd', 'subject_start_T1_rv', 'subject_end_T1_rv', 'subject_start_T2_rv', 'subject_end_T2_rv']:
        df_blast_1[i] = pd.to_numeric(df_blast_1[i])

    for i in ['query_start_T1_fd', 'query_end_T1_fd', 'query_start_T2_fd', 'query_end_T2_fd', 'query_start_T1_rv', 'query_end_T1_rv', 'query_start_T2_rv', 'query_end_T2_rv']:
        df_blast_1[i] = pd.to_numeric(df_blast_1[i])

    df_blast_1['subject_start_T1_fd'], df_blast_1['subject_end_T1_fd'] = np.where(df_blast_1['subject_start_T1_fd'] > df_blast_1['subject_end_T1_fd'], (df_blast_1['subject_end_T1_fd'], df_blast_1['subject_start_T1_fd']), (df_blast_1['subject_start_T1_fd'], df_blast_1['subject_end_T1_fd']))
    df_blast_1['subject_start_T2_fd'], df_blast_1['subject_end_T2_fd'] = np.where(df_blast_1['subject_start_T2_fd'] > df_blast_1['subject_end_T2_fd'], (df_blast_1['subject_end_T2_fd'], df_blast_1['subject_start_T2_fd']), (df_blast_1['subject_start_T2_fd'], df_blast_1['subject_end_T2_fd']))
    df_blast_1['subject_start_T1_rv'], df_blast_1['subject_end_T1_rv'] = np.where(df_blast_1['subject_start_T1_rv'] > df_blast_1['subject_end_T1_rv'], (df_blast_1['subject_end_T1_rv'], df_blast_1['subject_start_T1_rv']), (df_blast_1['subject_start_T1_rv'], df_blast_1['subject_end_T1_rv']))
    df_blast_1['subject_start_T2_rv'], df_blast_1['subject_end_T2_rv'] = np.where(df_blast_1['subject_start_T2_rv'] > df_blast_1['subject_end_T2_rv'], (df_blast_1['subject_end_T2_rv'], df_blast_1['subject_start_T2_rv']), (df_blast_1['subject_start_T2_rv'], df_blast_1['subject_end_T2_rv']))
    df_blast_1.head(4)

    """df_blast_1.loc[df_blast_1['query_ID'].str.contains('SRR11404272.10008287')]"""

    def find_overlap_T1(row):
        range1_start, range1_end, range2_start, range2_end = row['subject_start_T1_fd'], row['subject_end_T1_fd'], row['subject_start_T1_rv'], row['subject_end_T1_rv']
        overlap_start = max(range1_start, range2_start)
        overlap_end = min(range1_end, range2_end)
        if np.isnan(row['subject_start_T1_fd']) or np.isnan(row['subject_start_T1_rv']) or np.isnan(row['subject_end_T1_fd']) or np.isnan(row['subject_end_T1_rv']):
            return np.nan
        return overlap_end - overlap_start +  1 if overlap_start <= overlap_end else  0


    df_blast_1['overlap_subject_T1'] = df_blast_1.apply(find_overlap_T1, axis=1)
    # df_blast_1.head(3)

    # df_blast_1.tail()

    # df_blast_1.to_csv('df_blast_1.txt',sep='\t',index=False)

    # df_blast_1 = pd.read_csv('df_blast_1.txt',delimiter='\t')

    # df_blast_1.columns

    # st = df_blast_1[['query_ID', 'Transcript_ID_T1_fd','query_start_end_T1_fd', 'strand_T1_fd', 'Transcript_ID_T2_fd', 'strand_T2_fd', 'Transcript_ID_T1_rv','strand_T1_rv', 'Transcript_ID_T2_rv','strand_T2_rv']]
    # samstrand_box2 = st.loc[(st['final_transcript_strand_1'] == '+') & (st['final_transcript_strand_2'] == '-') |(st['final_transcript_strand_1'] == '-') & (st['final_transcript_strand_2'] == '+') ]
    # samstrand_box2
    # st


    def find_overlap_T2(row):
        range1_start, range1_end, range2_start, range2_end = row['subject_start_T2_fd'], row['subject_end_T2_fd'], row['subject_start_T2_rv'], row['subject_end_T2_rv']
        overlap_start = max(range1_start, range2_start)
        overlap_end = min(range1_end, range2_end)
        if np.isnan(row['subject_start_T2_fd']) or np.isnan(row['subject_start_T2_rv']) or np.isnan(row['subject_end_T2_fd']) or np.isnan(row['subject_end_T2_rv']):
            return np.nan
        return overlap_end - overlap_start +  1 if overlap_start <= overlap_end else  0

    def find_overlap_query_fd(row):
        range1_start, range1_end, range2_start, range2_end = row['query_start_T1_fd'], row['query_end_T1_fd'], row['query_start_T2_fd'], row['query_end_T2_fd']
        overlap_start = max(range1_start, range2_start)
        overlap_end = min(range1_end, range2_end)
        if np.isnan(row['query_start_T1_fd']) or np.isnan(row['query_start_T2_fd']) or np.isnan(row['query_end_T1_fd']) or np.isnan(row['query_end_T2_fd']):
            return np.nan
        return overlap_end - overlap_start +  1 if overlap_start <= overlap_end else  0

    def find_overlap_query_rv(row):
        range1_start, range1_end, range2_start, range2_end = row['query_start_T1_rv'], row['query_end_T1_rv'], row['query_start_T2_rv'], row['query_end_T2_rv']
        overlap_start = max(range1_start, range2_start)
        overlap_end = min(range1_end, range2_end)
        if np.isnan(row['query_start_T1_rv']) or np.isnan(row['query_start_T2_rv']) or np.isnan(row['query_end_T1_rv']) or np.isnan(row['query_end_T2_rv']):
            return np.nan
        return overlap_end - overlap_start +  1 if overlap_start <= overlap_end else  0

    df_blast_1['overlap_subject_T2'] = df_blast_1.apply(find_overlap_T2, axis=1)
    df_blast_1['overlap_query_fd'] = df_blast_1.apply(find_overlap_query_fd, axis=1)
    df_blast_1['overlap_query_rv'] = df_blast_1.apply(find_overlap_query_rv, axis=1)
    # keep in mind that the overlaps in nromal and complex should be same
    #

    df_blast_1.columns

    df_blast_1.shape

    df_blast_1['gap_query_fd'] = df_blast_1.apply(lambda x: abs(min(x['query_end_T1_fd'], x['query_end_T2_fd']) - max(x['query_start_T1_fd'], x['query_start_T2_fd'])) if x['overlap_query_fd'] == 0 else 0, axis=1)
    df_blast_1['gap_query_rv'] = df_blast_1.apply(lambda x: abs(min(x['query_end_T1_rv'], x['query_end_T2_rv']) - max(x['query_start_T1_rv'], x['query_start_T2_rv'])) if x['overlap_query_rv'] == 0 else 0, axis=1)
    df_blast_1

    # replace NAN with 0
    df_blast_1['gap_query_fd'] = df_blast_1['gap_query_fd'].fillna(0)
    df_blast_1['gap_query_rv'] = df_blast_1['gap_query_rv'].fillna(0)

    # subtracting 1 from query if the vlaue if more than 0
    # df_blast_1['overlap_query_fd'] =
    df_blast_1['gap_query_fd'] = df_blast_1['gap_query_fd'].apply(lambda x: x - 1 if x > 0 else 0)
    df_blast_1['gap_query_rv'] = df_blast_1['gap_query_rv'].apply(lambda x: x - 1 if x > 0 else 0)
    df_blast_1.iloc[:20,-11:]



    df_blast_1['gap_subject_T1'] = df_blast_1.apply(lambda x: abs(min(x['subject_end_T1_fd'], x['subject_end_T1_rv']) - max(x['subject_start_T1_fd'], x['subject_start_T1_rv'])) if x['overlap_subject_T1'] == 0 else 0, axis=1)
    df_blast_1['gap_subject_T2'] = df_blast_1.apply(lambda x: abs(min(x['subject_end_T2_fd'], x['subject_end_T2_rv']) - max(x['subject_start_T2_fd'], x['subject_start_T2_rv'])) if x['overlap_subject_T2'] == 0 else 0, axis=1)
    df_blast_1.iloc[:, -9:] = df_blast_1.iloc[:, -9:].fillna(0)
    ##_____________________________________________________________________PARAMETER__________________OVERLAPGAP__________________________________________________________ (4% of 150)
    df_blast_1['gap_condition'] = df_blast_1.apply(lambda x: 'NO' if (x['gap_query_fd'] <= int(args.gap))  and x['gap_query_rv'] <= int(args.gap)  else 'YES', axis=1)
    # samstrand_box2.iloc[:,-9:]
    ##_____________________________________________________________________PARAMETER___________________OVERLAPQUERY_________________________________________________________ (30% of 150)
    df_blast_1['overlap_query_condition'] = df_blast_1.apply(lambda x: 'NO' if (x['overlap_query_fd'] <=  int(args.overlap_query) and x['overlap_query_rv'] <= int(args.overlap_query)) else 'YES', axis=1)
    # samstrand_box2.iloc[:,-11:]
    # samstrand_box2.apply(lambda x: 'NO' if (x['overlap_subject_T1'] == 0 and x['overlap_subject_T2'] == 0) else \
                                                            #    'NO' if ((x['overlap_subject_T1'] > 0 or x['overlap_subject_T2'] > 0) and (x['gap_subject_T1'] <= 5 and x['gap_subject_T2'] <= 5)) else 'YES', axis=1)
    # If both overlap_subject_T1 and overlap_subject_T2 are 0, and both gap_subject_T1 and gap_subject_T2 are less than or equal to 5, return 'NO'.
    # If either overlap_subject_T1 or overlap_subject_T2 is greater than 0, return 'NO'.
    ##_____________________________________________________________________PARAMETER________________OVERLAPSUBJECT/ ANCHOR LENGTH____________________________________________________________ (20% of 150)
    df_blast_1['overlap_subject_condition'] = df_blast_1.apply(lambda row: 'NO' if row['overlap_subject_T1'] == int(args.subject_query) and row['overlap_subject_T2'] == int(args.subject_query) and row['gap_subject_T1'] <=  int(args.gap_subject) and row['gap_subject_T2'] <= int(args.gap_subject) else \
        'NO' if row['overlap_subject_T1'] > int(args.subject_query) or row['overlap_subject_T2'] > int(args.subject_query) else 'YES', axis=1)

    # df_blast_1['gap_condition'] = df_blast_1.apply(lambda x: 'NO' if (x['gap_query_fd'] <= int(args.gap) and x['gap_query_rv'] <= int(args.gap)) else 'YES', axis=1)
    # df_blast_1['overlap_query_condition'] = df_blast_1.apply(lambda x: 'NO' if (x['overlap_query_fd'] <= int(args.overlap_query) and x['overlap_query_rv'] <= int(args.overlap_query)) else 'YES', axis=1)
    # df_blast_1['overlap_subject_condition'] = df_blast_1.apply(lambda x: 'NO' if (x['overlap_subject_T1'] == 0 and x['overlap_subject_T2'] == 0) else \
    #                                                         'NO' if ((x['overlap_subject_T1'] > 0 or x['overlap_subject_T2'] > 0) and (x['gap_subject_T1'] <= 5 and x['gap_subject_T2'] <= 5)) else 'YES', axis=1)
    # df_blast_1.columns

    df_blast_1['reads_skipped'] = df_blast_1.apply(lambda x: 'YES' if (x['overlap_query_condition'] == 'YES' or x['overlap_subject_condition'] == 'YES' or x['gap_condition'] == 'YES') else 'NO', axis=1)

    df_blast_skipped_3 = df_blast_1[df_blast_1['reads_skipped'] == 'YES']

    df_blast_1.drop(columns=[col for col in df_blast_1.columns if 'skipped' in col], inplace=True)
    df_blast_skipped_3.drop(columns=[col for col in df_blast_skipped_3.columns if 'skipped' in col], inplace=True)

    df_blast_1 = df_blast_1[~df_blast_1.isin(df_blast_skipped_3.to_dict(orient='list')).all(1)]

    df_blast_1['query_start_T1_fd'] = pd.to_numeric(df_blast_1['query_start_T1_fd'])
    df_blast_1['query_start_T2_fd'] = pd.to_numeric(df_blast_1['query_start_T2_fd'])

    df_blast_1['query_start_T1_rv'] = pd.to_numeric(df_blast_1['query_start_T1_rv'])
    df_blast_1['query_start_T2_rv'] = pd.to_numeric(df_blast_1['query_start_T2_rv'])
    df_blast_1.columns

    df_blast_1

    df_blast_1['HEAD_plus_fd'] = df_blast_1.apply(lambda x: x['Transcript_ID_T1_fd'] if (x['strand_symbol_1_fd'] == '+' and x['strand_symbol_2_fd'] == '+' \
                                                                                and x['query_start_T1_fd'] < x['query_start_T2_fd']) else x['Transcript_ID_T2_fd'] if (x['strand_symbol_1_fd'] == '+' and x['strand_symbol_2_fd'] == '+' \
                                                                                    and x['query_start_T1_fd'] > x['query_start_T2_fd'])  else "NONE", axis=1)


    df_blast_1['HEAD_plus_rv'] = df_blast_1.apply(lambda x: x['Transcript_ID_T1_rv'] if (x['strand_symbol_1_rv'] == '+' and x['strand_symbol_2_rv'] == '+' \
                                                                                    and x['query_start_T1_rv'] < x['query_start_T2_rv']) else x['Transcript_ID_T2_rv'] if (x['strand_symbol_1_rv'] == '+' and x['strand_symbol_2_rv'] == '+' \
                                                                                    and x['query_start_T1_rv'] > x['query_start_T2_rv']) else "NONE", axis=1)


    df_blast_1['HEAD_minus_fd'] = df_blast_1.apply(lambda x: x['Transcript_ID_T1_fd'] if (x['strand_symbol_1_fd'] == '-' and x['strand_symbol_2_fd'] == '-' \
                                                                                    and x['query_start_T1_fd'] > x['query_start_T2_fd']) else x['Transcript_ID_T2_fd'] if (x['strand_symbol_1_fd'] == '-' and x['strand_symbol_2_fd'] == '-' \
                                                                                    and x['query_start_T1_fd'] < x['query_start_T2_fd']) else "NONE", axis=1)
    df_blast_1['HEAD_minus_rv'] = df_blast_1.apply(lambda x: x['Transcript_ID_T1_rv'] if (x['strand_symbol_1_rv'] == '-' and x['strand_symbol_2_rv'] == '-' \
                                                                                    and x['query_start_T1_rv'] > x['query_start_T2_rv']) else x['Transcript_ID_T2_rv'] if (x['strand_symbol_1_rv'] == '-' and x['strand_symbol_2_rv'] == '-' \
                                                                                    and x['query_start_T1_rv'] < x['query_start_T2_rv']) else "NONE", axis=1)
    df_blast_1.shape

    df_blast_1

    df_blast_1["Head"]  = df_blast_1[['HEAD_plus_fd', 'HEAD_plus_rv', 'HEAD_minus_fd', 'HEAD_minus_rv']].apply(lambda x: "".join(list(filter(lambda y: y != 'NONE', x.unique()))), axis=1)

    df_blast_1[['Transcript_ID_T1_fd', 'Transcript_ID_T2_fd', 'Transcript_ID_T1_rv', 'Transcript_ID_T2_rv']] = df_blast_1[['Transcript_ID_T1_fd', 'Transcript_ID_T2_fd', 'Transcript_ID_T1_rv', 'Transcript_ID_T2_rv']].astype(str)
    df_blast_1["Tail"]  = df_blast_1[['Transcript_ID_T1_fd', 'Transcript_ID_T2_fd', 'Transcript_ID_T1_rv', 'Transcript_ID_T2_rv']].apply(lambda x: list(filter(lambda y: y != 'nan', x.unique())), axis=1)

    df_blast_1['Tail'] = df_blast_1.apply(lambda x: "".join(list(filter(lambda y: y not in x['Head'], x['Tail']))), axis=1)

    df_blast_1.drop(columns=[col for col in df_blast_1.columns if 'HEAD_' in col], inplace=True)

    df_blast_1_NaN = df_blast_1[df_blast_1.isna().any(axis=1)]
    df_blast_1_NoNaN = df_blast_1[~df_blast_1.isna().any(axis=1)]
    df_blast_1_NoNaN_subset = df_blast_1_NoNaN.copy()

    df_blast_1_NoNaN_subset['keep'] = df_blast_1_NoNaN_subset.apply(lambda x: 'reverse' if x['overlap_query_fd'] == 0 and x['overlap_query_rv'] > 0 \
                                                                    else 'forward' if x['overlap_query_fd'] > 0 and x['overlap_query_rv'] == 0 \
                                                                    else 'reverse' if x['overlap_query_fd'] > 0 and x['overlap_query_rv'] > 0 and x['overlap_query_fd'] > x['overlap_query_rv'] \
                                                                    else 'forward' if x['overlap_query_fd'] > 0 and x['overlap_query_rv'] > 0 and x['overlap_query_fd'] < x['overlap_query_rv'] \
                                                                    else 'reverse' if x['gap_query_fd'] > 0 and x['gap_query_rv'] > 0 and x['gap_query_fd'] > x['gap_query_rv'] \
                                                                    else 'forward' if x['gap_query_fd'] > 0 and x['gap_query_rv'] > 0 and x['gap_query_fd'] < x['gap_query_rv'] \
                                                                    else 'Either' if x['overlap_query_fd'] == 0 and x['overlap_query_rv'] == 0 and x['gap_query_fd'] == 0 and x['gap_query_rv'] == 0 \
                                                                    else 'reverse' if x['overlap_query_fd'] == x['overlap_query_rv'] and x['gap_query_fd'] == 0 and x['gap_query_rv'] > 0 \
                                                                    else 'forward' if x['overlap_query_fd'] == x['overlap_query_rv'] and x['gap_query_fd'] > 0 and x['gap_query_rv'] == 0 \
                                                                    else "Either", axis=1)

    df_blast_1_NaN['keep'] = df_blast_1_NaN.apply(lambda x: 'forward' if (x['Transcript_ID_T1_fd'] != 'nan' and x['Transcript_ID_T2_fd'] != 'nan') \
                                                            else 'reverse' if (x['Transcript_ID_T1_rv'] != 'nan' and x['Transcript_ID_T2_rv'] != 'nan') else 'Either', axis=1)

    df_blast_1 = pd.concat([df_blast_1_NoNaN_subset, df_blast_1_NaN], axis=0)
    df_blast_1.drop(columns=[col for col in df_blast_1.columns if 'condition' in col], inplace=True)
    df_blast_1.drop(columns=['keep'], inplace=True)

    df_blast_1['subject_start_HEAD_fd'] = df_blast_1.apply(lambda x: x['subject_start_T1_fd'] if x['Transcript_ID_T1_fd'] == x['Head'] \
                                                        else x['subject_start_T2_fd'] if x['Transcript_ID_T2_fd'] == x['Head'] else "NONE", axis=1)
    df_blast_1['subject_end_HEAD_fd'] = df_blast_1.apply(lambda x: x['subject_end_T1_fd'] if x['Transcript_ID_T1_fd'] == x['Head'] \
                                                            else x['subject_end_T2_fd'] if x['Transcript_ID_T2_fd'] == x['Head'] else "NONE", axis=1)
    df_blast_1['subject_start_HEAD_rv'] = df_blast_1.apply(lambda x: x['subject_start_T1_rv'] if x['Transcript_ID_T1_rv'] == x['Head'] \
                                                            else x['subject_start_T2_rv'] if x['Transcript_ID_T2_rv'] == x['Head'] else "NONE", axis=1)
    df_blast_1['subject_end_HEAD_rv'] = df_blast_1.apply(lambda x: x['subject_end_T1_rv'] if x['Transcript_ID_T1_rv'] == x['Head'] \
                                                                else x['subject_end_T2_rv'] if x['Transcript_ID_T2_rv'] == x['Head'] else "NONE", axis=1)

    df_blast_1['subject_start_TAIL_fd'] = df_blast_1.apply(lambda x: x['subject_start_T1_fd'] if x['Transcript_ID_T1_fd'] == x['Tail'] \
                                                            else x['subject_start_T2_fd'] if x['Transcript_ID_T2_fd'] == x['Tail'] else "NONE", axis=1)
    df_blast_1['subject_end_TAIL_fd'] = df_blast_1.apply(lambda x: x['subject_end_T1_fd'] if x['Transcript_ID_T1_fd'] == x['Tail'] \
                                                            else x['subject_end_T2_fd'] if x['Transcript_ID_T2_fd'] == x['Tail'] else "NONE", axis=1)
    df_blast_1['subject_start_TAIL_rv'] = df_blast_1.apply(lambda x: x['subject_start_T1_rv'] if x['Transcript_ID_T1_rv'] == x['Tail'] \
                                                                else x['subject_start_T2_rv'] if x['Transcript_ID_T2_rv'] == x['Tail'] else "NONE", axis=1)
    df_blast_1['subject_end_TAIL_rv'] = df_blast_1.apply(lambda x: x['subject_end_T1_rv'] if x['Transcript_ID_T1_rv'] == x['Tail'] \
                                else x['subject_end_T2_rv'] if x['Transcript_ID_T2_rv'] == x['Tail'] else "NONE", axis=1)


    df_blast_1['strand_HEAD_fd'] = df_blast_1.apply(lambda x: x['strand_symbol_1_fd'] if x['Transcript_ID_T1_fd'] == x['Head'] \
                                                        else x['strand_symbol_2_fd'] if x['Transcript_ID_T2_fd'] == x['Head'] else "NONE", axis=1)
    df_blast_1['strand_HEAD_rv'] = df_blast_1.apply(lambda x: x['strand_symbol_1_rv'] if x['Transcript_ID_T1_rv'] == x['Head'] \
                                                        else x['strand_symbol_2_rv'] if x['Transcript_ID_T2_rv'] == x['Head'] else "NONE", axis=1)
    df_blast_1['strand_TAIL_fd'] = df_blast_1.apply(lambda x: x['strand_symbol_1_fd'] if x['Transcript_ID_T1_fd'] == x['Tail'] \
                                                        else x['strand_symbol_2_fd'] if x['Transcript_ID_T2_fd'] == x['Tail'] else "NONE", axis=1)
    df_blast_1['strand_TAIL_rv'] = df_blast_1.apply(lambda x: x['strand_symbol_1_rv'] if x['Transcript_ID_T1_rv'] == x['Tail'] \
                                                        else x['strand_symbol_2_rv'] if x['Transcript_ID_T2_rv'] == x['Tail'] else "NONE", axis=1)

    df_blast_1['retain_HEAD'] = df_blast_1.apply(lambda x: 'YES' if (x['strand_HEAD_fd'] != x['strand_HEAD_rv'] and (x['subject_start_HEAD_fd'] == x['subject_start_HEAD_rv'] or x['subject_end_HEAD_fd'] == x['subject_end_HEAD_rv'])) else 'NO', axis=1)
    df_blast_1['retain_TAIL'] = df_blast_1.apply(lambda x: 'YES' if (x['strand_TAIL_fd'] != x['strand_TAIL_rv'] and (x['subject_start_TAIL_fd'] == x['subject_start_TAIL_rv'] or x['subject_end_TAIL_fd'] == x['subject_end_TAIL_rv'])) else 'NO', axis=1)

    df_blast_1['retain'] = df_blast_1.apply(lambda x: 'HEAD' if x['retain_HEAD'] == 'YES' and x['retain_TAIL'] == 'NO' \
                                            else 'TAIL' if x['retain_TAIL'] == 'YES' and x['retain_HEAD'] == 'NO' \
                                            else 'BOTH' if x['retain_HEAD'] == 'YES' and x['retain_TAIL'] == 'YES' else 'NONE', axis=1)

    df_blast_1['keep'] = df_blast_1.apply(lambda x: 'reverse' if x['overlap_query_fd'] == 0 and x['overlap_query_rv'] > 0 \
                                                else 'forward' if x['overlap_query_fd'] > 0 and x['overlap_query_rv'] == 0 \
                                                else 'reverse' if x['overlap_query_fd'] > 0 and x['overlap_query_rv'] > 0 and x['overlap_query_fd'] > x['overlap_query_rv'] \
                                                else 'forward' if x['overlap_query_fd'] > 0 and x['overlap_query_rv'] > 0 and x['overlap_query_fd'] < x['overlap_query_rv'] \
                                                else 'reverse' if x['gap_query_fd'] > 0 and x['gap_query_rv'] > 0 and x['gap_query_fd'] > x['gap_query_rv'] \
                                                else 'forward' if x['gap_query_fd'] > 0 and x['gap_query_rv'] > 0 and x['gap_query_fd'] < x['gap_query_rv'] \
                                                else 'Either' if x['overlap_query_fd'] == 0 and x['overlap_query_rv'] == 0 and x['gap_query_fd'] == 0 and x['gap_query_rv'] == 0 \
                                                else 'reverse' if x['overlap_query_fd'] == x['overlap_query_rv'] and x['gap_query_fd'] == 0 and x['gap_query_rv'] > 0 \
                                                else 'forward' if x['overlap_query_fd'] == x['overlap_query_rv'] and x['gap_query_fd'] > 0 and x['gap_query_rv'] == 0 \
                                                else "Either", axis=1)


    df_blast_1.to_csv('df_blast_1_2.txt',sep='\t',index=False)

    df_blast_1.loc[df_blast_1['keep'] == 'reverse', df_blast_1.columns[df_blast_1.columns.str.endswith('_fd')]] = np.nan
    df_blast_1.loc[df_blast_1['keep'] == 'forward', df_blast_1.columns[df_blast_1.columns.str.endswith('_rv')]] = np.nan
    df_blast_1.loc[df_blast_1['keep'] == 'Either', df_blast_1.columns[df_blast_1.columns.str.endswith('_rv')]] = np.nan

    df_blast_1.drop(columns=['keep'], inplace=True)

    df_blast_1.drop(columns=[col for col in df_blast_1.columns if 'HEAD' in col or 'TAIL' in col], inplace=True)

    df_blast_1_fd = df_blast_1[['query_ID'] + [col for col in df_blast_1.columns if col.endswith('_fd')] + ['Head', 'Tail', 'retain']]
    df_blast_1_rv = df_blast_1[['query_ID'] + [col for col in df_blast_1.columns if col.endswith('_rv')] + ['Head', 'Tail', 'retain']]

    df_blast_1_fd.dropna(inplace=True)
    df_blast_1_rv.dropna(inplace=True)

    df_blast_1_fd.columns = df_blast_1_fd.columns.str.replace('_fd', '')
    df_blast_1_rv.columns = df_blast_1_rv.columns.str.replace('_rv', '')
    df_blast_1 = pd.concat([df_blast_1_fd, df_blast_1_rv], axis=0)
    df_blast_1.rename(columns={'query_ID': 'ID_query'}, inplace=True)
    df_blast_1.drop(columns=[col for col in df_blast_1.columns if col.startswith('query_') or col.startswith('expect_') or col.startswith('identity_') or col.startswith('strand_symbol_') or col.startswith('frame_id') or col.startswith('subject_start_T1') or col.startswith('subject_end_T1') or col.startswith('subject_start_T2') or col.startswith('subject_end_T2')], inplace=True)

    df_blast_1_flip = df_blast_1[~df_blast_1['Transcript_ID_T1'].isin(df_blast_1['Head'])]
    df_blast_1 = df_blast_1[df_blast_1['Transcript_ID_T1'].isin(df_blast_1['Head'])]

    #  This line copies the values of columns indexed 1 to 3 (second to fourth columns) from df_blast_1_flip into a new DataFrame called temp_vals.
    temp_vals = df_blast_1_flip.iloc[:, 1:4].copy()
    #
    df_blast_1_flip.iloc[:, 1:4] = df_blast_1_flip.iloc[:, 4:7].values
    df_blast_1_flip.iloc[:, 4:7] = temp_vals.values
    df_blast_1 = pd.concat([df_blast_1, df_blast_1_flip], axis=0)

    df_blast_1.drop(columns=['Head', 'Tail'], inplace=True)
    df_blast_1[['q_ID_count', 'ID_query']] = df_blast_1['ID_query'].str.split('#', expand=True)



    df_blast_1.drop(columns=['q_ID_count'], inplace=True)
    df_blast_1.rename(columns={'ID_query': 'query_ID'}, inplace=True)
    df_blast_1.sort_values(['query_ID'], inplace=True)

    # df_blast_1.to_csv(dir_path + "splitReads.csv", index=False)
    df_blast_1

    len(df_blast_1)

    samstrand_box2

    # # for the breakpoint identification, it is important for the identification of breakpoint.

    samstrand_box2.head(2)

    samstrand_box2.head(3)

    # cmpl_all.to_csv('cmpl.txt',sep='\t', index=False)
    # samstrand_box2.to_csv('samstrand_box2',sep='\t', index=False)

    cmpl_all = cmpl_all.rename(columns={'headT': 'Transcript_ID_T1', 'tailT': 'Transcript_ID_T2'})
    cmpl_all[['id', 'query_ID']] = cmpl_all['query_ID'].str.split('#', expand=True)
    # cmpl_all.head(2)
    # # cmpl.head(2)

    samstrand_box2_strand =  pd.merge(samstrand_box2, cmpl_all, on=['query_ID' ,	'Transcript_ID_T1','Transcript_ID_T2'], how='left')
    samstrand_box2_strand = samstrand_box2_strand[['query_ID', 'Transcript_ID_T1', 'subject_start_end_T1','final_transcript_strand_1', 'Transcript_ID_T2', 'subject_start_end_T2', 'final_transcript_strand_2', 'overlap_query', 'gap_query', 'retain' ]]
    samstrand_box2_strand

    #  pd.merge(samstrand_box2, cmpl_all, on=['query_ID' ,'Transcript_ID_T1','Transcript_ID_T2'], how='left')

    # cmpl_all.head(2)
    #print(samstrand_box2.columns)
    #print(cmpl_all.columns)

    # pd.merge(samstrand_box2, cmpl_all,on=['query_ID','Transcript_ID_T1','Transcript_ID_T2'])

    # cmpl_all.loc[(cmpl_all['Transcript_ID_T1']=='AT3G28730.1') & (cmpl_all['Transcript_ID_T2']=='AT3G28740.1')]

    # cmpl_all.loc[(cmpl_all['Transcript_ID_T1']=='AT3G28740.1') & (cmpl_all['Transcript_ID_T2']=='AT3G28730.1')]





    # samstrand_box2.head(3)

    Nrm_all = Nrm_all.rename(columns={'headT': 'Transcript_ID_T1', 'tailT': 'Transcript_ID_T2'})
    Nrm_all[['id', 'query_ID']] = Nrm_all['query_ID'].str.split('#', expand=True)

    # Nrm_all = Nrm_all.rename(columns={'headT': 'Transcript_ID_T1', 'tailT': 'Transcript_ID_T2'})
    # Nrm_all[['id', 'query_ID']] = Nrm_all['query_ID'].str.split('#', expand=True)
    # Nrm_all.head(2)
    # # Nrm.head(2)

    Nrm.head(4)



    df_blast_1_strand =  pd.merge( df_blast_1, Nrm_all, on=['query_ID' ,	'Transcript_ID_T1','Transcript_ID_T2'], how='left')
    df_blast_1_strand = df_blast_1_strand[['query_ID', 'Transcript_ID_T1', 'subject_start_end_T1','final_transcript_strand_1', 'Transcript_ID_T2', 'subject_start_end_T2', 'final_transcript_strand_2', 'overlap_query', 'gap_query', 'retain' ]]
    df_blast_1_strand = df_blast_1_strand.drop_duplicates()
    df_blast_1_strand.shape

 

    df_blast_1.shape



    allsplit = pd.concat([samstrand_box2_strand,df_blast_1_strand])

else:

    df_blast_1 =samstrand.copy()
    df_blast_1.shape


    df_blast_1[['query_start_T1_fd', 'query_end_T1_fd']] = df_blast_1['query_start_end_T1_fd'].str.split('-', expand=True)
    df_blast_1[['query_start_T2_fd', 'query_end_T2_fd']] = df_blast_1['query_start_end_T2_fd'].str.split('-', expand=True)
    df_blast_1[['query_start_T1_rv', 'query_end_T1_rv']] = df_blast_1['query_start_end_T1_rv'].str.split('-', expand=True)
    df_blast_1[['query_start_T2_rv', 'query_end_T2_rv']] = df_blast_1['query_start_end_T2_rv'].str.split('-', expand=True)
    df_blast_1.head(4)

    df_blast_1[['subject_start_T1_fd', 'subject_end_T1_fd']] = df_blast_1['subject_start_end_T1_fd'].str.split('-', expand=True)
    df_blast_1[['subject_start_T2_fd', 'subject_end_T2_fd']] = df_blast_1['subject_start_end_T2_fd'].str.split('-', expand=True)
    df_blast_1[['subject_start_T1_rv', 'subject_end_T1_rv']] = df_blast_1['subject_start_end_T1_rv'].str.split('-', expand=True)
    df_blast_1[['subject_start_T2_rv', 'subject_end_T2_rv']] = df_blast_1['subject_start_end_T2_rv'].str.split('-', expand=True)

    for i in ['subject_start_T1_fd', 'subject_end_T1_fd', 'subject_start_T2_fd', 'subject_end_T2_fd', 'subject_start_T1_rv', 'subject_end_T1_rv', 'subject_start_T2_rv', 'subject_end_T2_rv']:
        df_blast_1[i] = pd.to_numeric(df_blast_1[i])

    for i in ['query_start_T1_fd', 'query_end_T1_fd', 'query_start_T2_fd', 'query_end_T2_fd', 'query_start_T1_rv', 'query_end_T1_rv', 'query_start_T2_rv', 'query_end_T2_rv']:
        df_blast_1[i] = pd.to_numeric(df_blast_1[i])

    df_blast_1['subject_start_T1_fd'], df_blast_1['subject_end_T1_fd'] = np.where(df_blast_1['subject_start_T1_fd'] > df_blast_1['subject_end_T1_fd'], (df_blast_1['subject_end_T1_fd'], df_blast_1['subject_start_T1_fd']), (df_blast_1['subject_start_T1_fd'], df_blast_1['subject_end_T1_fd']))
    df_blast_1['subject_start_T2_fd'], df_blast_1['subject_end_T2_fd'] = np.where(df_blast_1['subject_start_T2_fd'] > df_blast_1['subject_end_T2_fd'], (df_blast_1['subject_end_T2_fd'], df_blast_1['subject_start_T2_fd']), (df_blast_1['subject_start_T2_fd'], df_blast_1['subject_end_T2_fd']))
    df_blast_1['subject_start_T1_rv'], df_blast_1['subject_end_T1_rv'] = np.where(df_blast_1['subject_start_T1_rv'] > df_blast_1['subject_end_T1_rv'], (df_blast_1['subject_end_T1_rv'], df_blast_1['subject_start_T1_rv']), (df_blast_1['subject_start_T1_rv'], df_blast_1['subject_end_T1_rv']))
    df_blast_1['subject_start_T2_rv'], df_blast_1['subject_end_T2_rv'] = np.where(df_blast_1['subject_start_T2_rv'] > df_blast_1['subject_end_T2_rv'], (df_blast_1['subject_end_T2_rv'], df_blast_1['subject_start_T2_rv']), (df_blast_1['subject_start_T2_rv'], df_blast_1['subject_end_T2_rv']))
    df_blast_1.head(4)

    """df_blast_1.loc[df_blast_1['query_ID'].str.contains('SRR11404272.10008287')]"""

    def find_overlap_T1(row):
        range1_start, range1_end, range2_start, range2_end = row['subject_start_T1_fd'], row['subject_end_T1_fd'], row['subject_start_T1_rv'], row['subject_end_T1_rv']
        overlap_start = max(range1_start, range2_start)
        overlap_end = min(range1_end, range2_end)
        if np.isnan(row['subject_start_T1_fd']) or np.isnan(row['subject_start_T1_rv']) or np.isnan(row['subject_end_T1_fd']) or np.isnan(row['subject_end_T1_rv']):
            return np.nan
        return overlap_end - overlap_start +  1 if overlap_start <= overlap_end else  0


    df_blast_1['overlap_subject_T1'] = df_blast_1.apply(find_overlap_T1, axis=1)

    df_blast_1.tail()

    df_blast_1.to_csv('df_blast_1.txt',sep='\t',index=False)

    df_blast_1 = pd.read_csv('df_blast_1.txt',delimiter='\t')

    df_blast_1.columns

    # st = df_blast_1[['query_ID', 'Transcript_ID_T1_fd','query_start_end_T1_fd', 'strand_T1_fd', 'Transcript_ID_T2_fd', 'strand_T2_fd', 'Transcript_ID_T1_rv','strand_T1_rv', 'Transcript_ID_T2_rv','strand_T2_rv']]
    # samstrand_box2 = st.loc[(st['final_transcript_strand_1'] == '+') & (st['final_transcript_strand_2'] == '-') |(st['final_transcript_strand_1'] == '-') & (st['final_transcript_strand_2'] == '+') ]
    # samstrand_box2
    # st







    def find_overlap_T2(row):
        range1_start, range1_end, range2_start, range2_end = row['subject_start_T2_fd'], row['subject_end_T2_fd'], row['subject_start_T2_rv'], row['subject_end_T2_rv']
        overlap_start = max(range1_start, range2_start)
        overlap_end = min(range1_end, range2_end)
        if np.isnan(row['subject_start_T2_fd']) or np.isnan(row['subject_start_T2_rv']) or np.isnan(row['subject_end_T2_fd']) or np.isnan(row['subject_end_T2_rv']):
            return np.nan
        return overlap_end - overlap_start +  1 if overlap_start <= overlap_end else  0

    def find_overlap_query_fd(row):
        range1_start, range1_end, range2_start, range2_end = row['query_start_T1_fd'], row['query_end_T1_fd'], row['query_start_T2_fd'], row['query_end_T2_fd']
        overlap_start = max(range1_start, range2_start)
        overlap_end = min(range1_end, range2_end)
        if np.isnan(row['query_start_T1_fd']) or np.isnan(row['query_start_T2_fd']) or np.isnan(row['query_end_T1_fd']) or np.isnan(row['query_end_T2_fd']):
            return np.nan
        return overlap_end - overlap_start +  1 if overlap_start <= overlap_end else  0

    def find_overlap_query_rv(row):
        range1_start, range1_end, range2_start, range2_end = row['query_start_T1_rv'], row['query_end_T1_rv'], row['query_start_T2_rv'], row['query_end_T2_rv']
        overlap_start = max(range1_start, range2_start)
        overlap_end = min(range1_end, range2_end)
        if np.isnan(row['query_start_T1_rv']) or np.isnan(row['query_start_T2_rv']) or np.isnan(row['query_end_T1_rv']) or np.isnan(row['query_end_T2_rv']):
            return np.nan
        return overlap_end - overlap_start +  1 if overlap_start <= overlap_end else  0

    df_blast_1['overlap_subject_T2'] = df_blast_1.apply(find_overlap_T2, axis=1)
    df_blast_1['overlap_query_fd'] = df_blast_1.apply(find_overlap_query_fd, axis=1)
    df_blast_1['overlap_query_rv'] = df_blast_1.apply(find_overlap_query_rv, axis=1)
    # keep in mind that the overlaps in nromal and complex should be same
    #

    df_blast_1.columns

    df_blast_1.shape

    df_blast_1['gap_query_fd'] = df_blast_1.apply(lambda x: abs(min(x['query_end_T1_fd'], x['query_end_T2_fd']) - max(x['query_start_T1_fd'], x['query_start_T2_fd'])) if x['overlap_query_fd'] == 0 else 0, axis=1)
    df_blast_1['gap_query_rv'] = df_blast_1.apply(lambda x: abs(min(x['query_end_T1_rv'], x['query_end_T2_rv']) - max(x['query_start_T1_rv'], x['query_start_T2_rv'])) if x['overlap_query_rv'] == 0 else 0, axis=1)
    df_blast_1

    # replace NAN with 0
    df_blast_1['gap_query_fd'] = df_blast_1['gap_query_fd'].fillna(0)
    df_blast_1['gap_query_rv'] = df_blast_1['gap_query_rv'].fillna(0)

    # subtracting 1 from query if the vlaue if more than 0
    # df_blast_1['overlap_query_fd'] =
    df_blast_1['gap_query_fd'] = df_blast_1['gap_query_fd'].apply(lambda x: x - 1 if x > 0 else 0)
    df_blast_1['gap_query_rv'] = df_blast_1['gap_query_rv'].apply(lambda x: x - 1 if x > 0 else 0)
    df_blast_1.iloc[:20,-11:]



    df_blast_1['gap_subject_T1'] = df_blast_1.apply(lambda x: abs(min(x['subject_end_T1_fd'], x['subject_end_T1_rv']) - max(x['subject_start_T1_fd'], x['subject_start_T1_rv'])) if x['overlap_subject_T1'] == 0 else 0, axis=1)
    df_blast_1['gap_subject_T2'] = df_blast_1.apply(lambda x: abs(min(x['subject_end_T2_fd'], x['subject_end_T2_rv']) - max(x['subject_start_T2_fd'], x['subject_start_T2_rv'])) if x['overlap_subject_T2'] == 0 else 0, axis=1)
    df_blast_1.iloc[:, -9:] = df_blast_1.iloc[:, -9:].fillna(0)
    ##_____________________________________________________________________PARAMETER__________________OVERLAPGAP__________________________________________________________ (4% of 150)
    df_blast_1['gap_condition'] = df_blast_1.apply(lambda x: 'NO' if (x['gap_query_fd'] <= int(args.gap))  and x['gap_query_rv'] <= int(args.gap)  else 'YES', axis=1)
    # samstrand_box2.iloc[:,-9:]
    ##_____________________________________________________________________PARAMETER___________________OVERLAPQUERY_________________________________________________________ (30% of 150)
    df_blast_1['overlap_query_condition'] = df_blast_1.apply(lambda x: 'NO' if (x['overlap_query_fd'] <=  int(args.overlap_query) and x['overlap_query_rv'] <= int(args.overlap_query)) else 'YES', axis=1)
    # samstrand_box2.iloc[:,-11:]
    # samstrand_box2.apply(lambda x: 'NO' if (x['overlap_subject_T1'] == 0 and x['overlap_subject_T2'] == 0) else \
                                                            #    'NO' if ((x['overlap_subject_T1'] > 0 or x['overlap_subject_T2'] > 0) and (x['gap_subject_T1'] <= 5 and x['gap_subject_T2'] <= 5)) else 'YES', axis=1)
    # If both overlap_subject_T1 and overlap_subject_T2 are 0, and both gap_subject_T1 and gap_subject_T2 are less than or equal to 5, return 'NO'.
    # If either overlap_subject_T1 or overlap_subject_T2 is greater than 0, return 'NO'.
    ##_____________________________________________________________________PARAMETER________________OVERLAPSUBJECT/ ANCHOR LENGTH____________________________________________________________ (20% of 150)
    df_blast_1['overlap_subject_condition'] = df_blast_1.apply(lambda row: 'NO' if row['overlap_subject_T1'] == int(args.subject_query) and row['overlap_subject_T2'] == int(args.subject_query) and row['gap_subject_T1'] <=  int(args.gap_subject) and row['gap_subject_T2'] <= int(args.gap_subject) else \
        'NO' if row['overlap_subject_T1'] > int(args.subject_query) or row['overlap_subject_T2'] > int(args.subject_query) else 'YES', axis=1)

    df_blast_1['reads_skipped'] = df_blast_1.apply(lambda x: 'YES' if (x['overlap_query_condition'] == 'YES' or x['overlap_subject_condition'] == 'YES' or x['gap_condition'] == 'YES') else 'NO', axis=1)

    df_blast_skipped_3 = df_blast_1[df_blast_1['reads_skipped'] == 'YES']

    df_blast_1.drop(columns=[col for col in df_blast_1.columns if 'skipped' in col], inplace=True)
    df_blast_skipped_3.drop(columns=[col for col in df_blast_skipped_3.columns if 'skipped' in col], inplace=True)

    df_blast_1 = df_blast_1[~df_blast_1.isin(df_blast_skipped_3.to_dict(orient='list')).all(1)]

    df_blast_1['query_start_T1_fd'] = pd.to_numeric(df_blast_1['query_start_T1_fd'])
    df_blast_1['query_start_T2_fd'] = pd.to_numeric(df_blast_1['query_start_T2_fd'])

    df_blast_1['query_start_T1_rv'] = pd.to_numeric(df_blast_1['query_start_T1_rv'])
    df_blast_1['query_start_T2_rv'] = pd.to_numeric(df_blast_1['query_start_T2_rv'])
    df_blast_1.columns

    df_blast_1

    df_blast_1['HEAD_plus_fd'] = df_blast_1.apply(lambda x: x['Transcript_ID_T1_fd'] if (x['strand_symbol_1_fd'] == '+' and x['strand_symbol_2_fd'] == '+' \
                                                                                and x['query_start_T1_fd'] < x['query_start_T2_fd']) else x['Transcript_ID_T2_fd'] if (x['strand_symbol_1_fd'] == '+' and x['strand_symbol_2_fd'] == '+' \
                                                                                    and x['query_start_T1_fd'] > x['query_start_T2_fd'])  else "NONE", axis=1)


    df_blast_1['HEAD_plus_rv'] = df_blast_1.apply(lambda x: x['Transcript_ID_T1_rv'] if (x['strand_symbol_1_rv'] == '+' and x['strand_symbol_2_rv'] == '+' \
                                                                                    and x['query_start_T1_rv'] < x['query_start_T2_rv']) else x['Transcript_ID_T2_rv'] if (x['strand_symbol_1_rv'] == '+' and x['strand_symbol_2_rv'] == '+' \
                                                                                    and x['query_start_T1_rv'] > x['query_start_T2_rv']) else "NONE", axis=1)


    df_blast_1['HEAD_minus_fd'] = df_blast_1.apply(lambda x: x['Transcript_ID_T1_fd'] if (x['strand_symbol_1_fd'] == '-' and x['strand_symbol_2_fd'] == '-' \
                                                                                    and x['query_start_T1_fd'] > x['query_start_T2_fd']) else x['Transcript_ID_T2_fd'] if (x['strand_symbol_1_fd'] == '-' and x['strand_symbol_2_fd'] == '-' \
                                                                                    and x['query_start_T1_fd'] < x['query_start_T2_fd']) else "NONE", axis=1)
    df_blast_1['HEAD_minus_rv'] = df_blast_1.apply(lambda x: x['Transcript_ID_T1_rv'] if (x['strand_symbol_1_rv'] == '-' and x['strand_symbol_2_rv'] == '-' \
                                                                                    and x['query_start_T1_rv'] > x['query_start_T2_rv']) else x['Transcript_ID_T2_rv'] if (x['strand_symbol_1_rv'] == '-' and x['strand_symbol_2_rv'] == '-' \
                                                                                    and x['query_start_T1_rv'] < x['query_start_T2_rv']) else "NONE", axis=1)
    df_blast_1.shape

    df_blast_1

    df_blast_1["Head"]  = df_blast_1[['HEAD_plus_fd', 'HEAD_plus_rv', 'HEAD_minus_fd', 'HEAD_minus_rv']].apply(lambda x: "".join(list(filter(lambda y: y != 'NONE', x.unique()))), axis=1)

    df_blast_1[['Transcript_ID_T1_fd', 'Transcript_ID_T2_fd', 'Transcript_ID_T1_rv', 'Transcript_ID_T2_rv']] = df_blast_1[['Transcript_ID_T1_fd', 'Transcript_ID_T2_fd', 'Transcript_ID_T1_rv', 'Transcript_ID_T2_rv']].astype(str)
    df_blast_1["Tail"]  = df_blast_1[['Transcript_ID_T1_fd', 'Transcript_ID_T2_fd', 'Transcript_ID_T1_rv', 'Transcript_ID_T2_rv']].apply(lambda x: list(filter(lambda y: y != 'nan', x.unique())), axis=1)

    df_blast_1['Tail'] = df_blast_1.apply(lambda x: "".join(list(filter(lambda y: y not in x['Head'], x['Tail']))), axis=1)

    df_blast_1.drop(columns=[col for col in df_blast_1.columns if 'HEAD_' in col], inplace=True)

    df_blast_1_NaN = df_blast_1[df_blast_1.isna().any(axis=1)]
    df_blast_1_NoNaN = df_blast_1[~df_blast_1.isna().any(axis=1)]
    df_blast_1_NoNaN_subset = df_blast_1_NoNaN.copy()

    df_blast_1_NoNaN_subset['keep'] = df_blast_1_NoNaN_subset.apply(lambda x: 'reverse' if x['overlap_query_fd'] == 0 and x['overlap_query_rv'] > 0 \
                                                                    else 'forward' if x['overlap_query_fd'] > 0 and x['overlap_query_rv'] == 0 \
                                                                    else 'reverse' if x['overlap_query_fd'] > 0 and x['overlap_query_rv'] > 0 and x['overlap_query_fd'] > x['overlap_query_rv'] \
                                                                    else 'forward' if x['overlap_query_fd'] > 0 and x['overlap_query_rv'] > 0 and x['overlap_query_fd'] < x['overlap_query_rv'] \
                                                                    else 'reverse' if x['gap_query_fd'] > 0 and x['gap_query_rv'] > 0 and x['gap_query_fd'] > x['gap_query_rv'] \
                                                                    else 'forward' if x['gap_query_fd'] > 0 and x['gap_query_rv'] > 0 and x['gap_query_fd'] < x['gap_query_rv'] \
                                                                    else 'Either' if x['overlap_query_fd'] == 0 and x['overlap_query_rv'] == 0 and x['gap_query_fd'] == 0 and x['gap_query_rv'] == 0 \
                                                                    else 'reverse' if x['overlap_query_fd'] == x['overlap_query_rv'] and x['gap_query_fd'] == 0 and x['gap_query_rv'] > 0 \
                                                                    else 'forward' if x['overlap_query_fd'] == x['overlap_query_rv'] and x['gap_query_fd'] > 0 and x['gap_query_rv'] == 0 \
                                                                    else "Either", axis=1)

    df_blast_1_NaN['keep'] = df_blast_1_NaN.apply(lambda x: 'forward' if (x['Transcript_ID_T1_fd'] != 'nan' and x['Transcript_ID_T2_fd'] != 'nan') \
                                                            else 'reverse' if (x['Transcript_ID_T1_rv'] != 'nan' and x['Transcript_ID_T2_rv'] != 'nan') else 'Either', axis=1)

    df_blast_1 = pd.concat([df_blast_1_NoNaN_subset, df_blast_1_NaN], axis=0)
    df_blast_1.drop(columns=[col for col in df_blast_1.columns if 'condition' in col], inplace=True)
    df_blast_1.drop(columns=['keep'], inplace=True)

    df_blast_1['subject_start_HEAD_fd'] = df_blast_1.apply(lambda x: x['subject_start_T1_fd'] if x['Transcript_ID_T1_fd'] == x['Head'] \
                                                        else x['subject_start_T2_fd'] if x['Transcript_ID_T2_fd'] == x['Head'] else "NONE", axis=1)
    df_blast_1['subject_end_HEAD_fd'] = df_blast_1.apply(lambda x: x['subject_end_T1_fd'] if x['Transcript_ID_T1_fd'] == x['Head'] \
                                                            else x['subject_end_T2_fd'] if x['Transcript_ID_T2_fd'] == x['Head'] else "NONE", axis=1)
    df_blast_1['subject_start_HEAD_rv'] = df_blast_1.apply(lambda x: x['subject_start_T1_rv'] if x['Transcript_ID_T1_rv'] == x['Head'] \
                                                            else x['subject_start_T2_rv'] if x['Transcript_ID_T2_rv'] == x['Head'] else "NONE", axis=1)
    df_blast_1['subject_end_HEAD_rv'] = df_blast_1.apply(lambda x: x['subject_end_T1_rv'] if x['Transcript_ID_T1_rv'] == x['Head'] \
                                                                else x['subject_end_T2_rv'] if x['Transcript_ID_T2_rv'] == x['Head'] else "NONE", axis=1)

    df_blast_1['subject_start_TAIL_fd'] = df_blast_1.apply(lambda x: x['subject_start_T1_fd'] if x['Transcript_ID_T1_fd'] == x['Tail'] \
                                                            else x['subject_start_T2_fd'] if x['Transcript_ID_T2_fd'] == x['Tail'] else "NONE", axis=1)
    df_blast_1['subject_end_TAIL_fd'] = df_blast_1.apply(lambda x: x['subject_end_T1_fd'] if x['Transcript_ID_T1_fd'] == x['Tail'] \
                                                            else x['subject_end_T2_fd'] if x['Transcript_ID_T2_fd'] == x['Tail'] else "NONE", axis=1)
    df_blast_1['subject_start_TAIL_rv'] = df_blast_1.apply(lambda x: x['subject_start_T1_rv'] if x['Transcript_ID_T1_rv'] == x['Tail'] \
                                                                else x['subject_start_T2_rv'] if x['Transcript_ID_T2_rv'] == x['Tail'] else "NONE", axis=1)
    df_blast_1['subject_end_TAIL_rv'] = df_blast_1.apply(lambda x: x['subject_end_T1_rv'] if x['Transcript_ID_T1_rv'] == x['Tail'] \
                                else x['subject_end_T2_rv'] if x['Transcript_ID_T2_rv'] == x['Tail'] else "NONE", axis=1)


    df_blast_1['strand_HEAD_fd'] = df_blast_1.apply(lambda x: x['strand_symbol_1_fd'] if x['Transcript_ID_T1_fd'] == x['Head'] \
                                                        else x['strand_symbol_2_fd'] if x['Transcript_ID_T2_fd'] == x['Head'] else "NONE", axis=1)
    df_blast_1['strand_HEAD_rv'] = df_blast_1.apply(lambda x: x['strand_symbol_1_rv'] if x['Transcript_ID_T1_rv'] == x['Head'] \
                                                        else x['strand_symbol_2_rv'] if x['Transcript_ID_T2_rv'] == x['Head'] else "NONE", axis=1)
    df_blast_1['strand_TAIL_fd'] = df_blast_1.apply(lambda x: x['strand_symbol_1_fd'] if x['Transcript_ID_T1_fd'] == x['Tail'] \
                                                        else x['strand_symbol_2_fd'] if x['Transcript_ID_T2_fd'] == x['Tail'] else "NONE", axis=1)
    df_blast_1['strand_TAIL_rv'] = df_blast_1.apply(lambda x: x['strand_symbol_1_rv'] if x['Transcript_ID_T1_rv'] == x['Tail'] \
                                                        else x['strand_symbol_2_rv'] if x['Transcript_ID_T2_rv'] == x['Tail'] else "NONE", axis=1)

    df_blast_1['retain_HEAD'] = df_blast_1.apply(lambda x: 'YES' if (x['strand_HEAD_fd'] != x['strand_HEAD_rv'] and (x['subject_start_HEAD_fd'] == x['subject_start_HEAD_rv'] or x['subject_end_HEAD_fd'] == x['subject_end_HEAD_rv'])) else 'NO', axis=1)
    df_blast_1['retain_TAIL'] = df_blast_1.apply(lambda x: 'YES' if (x['strand_TAIL_fd'] != x['strand_TAIL_rv'] and (x['subject_start_TAIL_fd'] == x['subject_start_TAIL_rv'] or x['subject_end_TAIL_fd'] == x['subject_end_TAIL_rv'])) else 'NO', axis=1)

    df_blast_1['retain'] = df_blast_1.apply(lambda x: 'HEAD' if x['retain_HEAD'] == 'YES' and x['retain_TAIL'] == 'NO' \
                                            else 'TAIL' if x['retain_TAIL'] == 'YES' and x['retain_HEAD'] == 'NO' \
                                            else 'BOTH' if x['retain_HEAD'] == 'YES' and x['retain_TAIL'] == 'YES' else 'NONE', axis=1)

    df_blast_1['keep'] = df_blast_1.apply(lambda x: 'reverse' if x['overlap_query_fd'] == 0 and x['overlap_query_rv'] > 0 \
                                                else 'forward' if x['overlap_query_fd'] > 0 and x['overlap_query_rv'] == 0 \
                                                else 'reverse' if x['overlap_query_fd'] > 0 and x['overlap_query_rv'] > 0 and x['overlap_query_fd'] > x['overlap_query_rv'] \
                                                else 'forward' if x['overlap_query_fd'] > 0 and x['overlap_query_rv'] > 0 and x['overlap_query_fd'] < x['overlap_query_rv'] \
                                                else 'reverse' if x['gap_query_fd'] > 0 and x['gap_query_rv'] > 0 and x['gap_query_fd'] > x['gap_query_rv'] \
                                                else 'forward' if x['gap_query_fd'] > 0 and x['gap_query_rv'] > 0 and x['gap_query_fd'] < x['gap_query_rv'] \
                                                else 'Either' if x['overlap_query_fd'] == 0 and x['overlap_query_rv'] == 0 and x['gap_query_fd'] == 0 and x['gap_query_rv'] == 0 \
                                                else 'reverse' if x['overlap_query_fd'] == x['overlap_query_rv'] and x['gap_query_fd'] == 0 and x['gap_query_rv'] > 0 \
                                                else 'forward' if x['overlap_query_fd'] == x['overlap_query_rv'] and x['gap_query_fd'] > 0 and x['gap_query_rv'] == 0 \
                                                else "Either", axis=1)


    #df_blast_1.to_csv('df_blast_1_2.txt',sep='\t',index=False)
    df_blast_1.loc[df_blast_1['keep'] == 'reverse', df_blast_1.columns[df_blast_1.columns.str.endswith('_fd')]] = np.nan
    df_blast_1.loc[df_blast_1['keep'] == 'forward', df_blast_1.columns[df_blast_1.columns.str.endswith('_rv')]] = np.nan
    df_blast_1.loc[df_blast_1['keep'] == 'Either', df_blast_1.columns[df_blast_1.columns.str.endswith('_rv')]] = np.nan
    df_blast_1.drop(columns=['keep'], inplace=True)
    df_blast_1.drop(columns=[col for col in df_blast_1.columns if 'HEAD' in col or 'TAIL' in col], inplace=True)
    df_blast_1_fd = df_blast_1[['query_ID'] + [col for col in df_blast_1.columns if col.endswith('_fd')] + ['Head', 'Tail', 'retain']]
    df_blast_1_rv = df_blast_1[['query_ID'] + [col for col in df_blast_1.columns if col.endswith('_rv')] + ['Head', 'Tail', 'retain']]
    df_blast_1_fd.dropna(inplace=True)
    df_blast_1_rv.dropna(inplace=True)
    df_blast_1_fd.columns = df_blast_1_fd.columns.str.replace('_fd', '')
    df_blast_1_rv.columns = df_blast_1_rv.columns.str.replace('_rv', '')
    df_blast_1 = pd.concat([df_blast_1_fd, df_blast_1_rv], axis=0)
    df_blast_1.rename(columns={'query_ID': 'ID_query'}, inplace=True)
    df_blast_1.drop(columns=[col for col in df_blast_1.columns if col.startswith('query_') or col.startswith('expect_') or col.startswith('identity_') or col.startswith('strand_symbol_') or col.startswith('frame_id') or col.startswith('subject_start_T1') or col.startswith('subject_end_T1') or col.startswith('subject_start_T2') or col.startswith('subject_end_T2')], inplace=True)
    df_blast_1_flip = df_blast_1[~df_blast_1['Transcript_ID_T1'].isin(df_blast_1['Head'])]
    df_blast_1 = df_blast_1[df_blast_1['Transcript_ID_T1'].isin(df_blast_1['Head'])]
    #  This line copies the values of columns indexed 1 to 3 (second to fourth columns) from df_blast_1_flip into a new DataFrame called temp_vals.
    temp_vals = df_blast_1_flip.iloc[:, 1:4].copy()
    df_blast_1_flip.iloc[:, 1:4] = df_blast_1_flip.iloc[:, 4:7].values
    df_blast_1_flip.iloc[:, 4:7] = temp_vals.values
    df_blast_1 = pd.concat([df_blast_1, df_blast_1_flip], axis=0)
    df_blast_1.drop(columns=['Head', 'Tail'], inplace=True)
    df_blast_1[['q_ID_count', 'ID_query']] = df_blast_1['ID_query'].str.split('#', expand=True)
    df_blast_1.drop(columns=['q_ID_count'], inplace=True)
    df_blast_1.rename(columns={'ID_query': 'query_ID'}, inplace=True)
    df_blast_1.sort_values(['query_ID'], inplace=True)
    Nrm_all = Nrm_all.rename(columns={'headT': 'Transcript_ID_T1', 'tailT': 'Transcript_ID_T2'})
    Nrm_all[['id', 'query_ID']] = Nrm_all['query_ID'].str.split('#', expand=True)    
    df_blast_1_strand =  pd.merge( df_blast_1, Nrm_all, on=['query_ID' ,	'Transcript_ID_T1','Transcript_ID_T2'], how='left')
    df_blast_1_strand = df_blast_1_strand[['query_ID', 'Transcript_ID_T1', 'subject_start_end_T1','final_transcript_strand_1', 'Transcript_ID_T2', 'subject_start_end_T2', 'final_transcript_strand_2', 'overlap_query', 'gap_query', 'retain' ]]
    df_blast_1_strand = df_blast_1_strand.drop_duplicates()
    df_blast_1_strand.shape
    allsplit = df_blast_1_strand 


allsplit = allsplit.drop_duplicates()
# Tobe Added
## I think it's equally import to know that find out the total number of the bases from the total that are mapping to the transcripts

# TotalReadLength = args.Rlength

# Function to calculate differences and subtract overlap_query
def calculate_difference(value, overlap_query):
    try:
        start, end = map(int, value.split('-'))
        return abs(start - end) - overlap_query
    except:
        return None

# Apply the function to the columns
allsplit['T1_dif'] = allsplit.apply(lambda row: calculate_difference(row['subject_start_end_T1'], row['overlap_query']), axis=1)
allsplit['T2_dif'] = allsplit.apply(lambda row: calculate_difference(row['subject_start_end_T2'], row['overlap_query']), axis=1)
# # Calculate the percentage of the total read length
# allsplit['T1_perc'] = (allsplit['T1_dif'] * 100) / TotalReadLength
# allsplit['T2_perc'] = (allsplit['T2_dif'] * 100) / TotalReadLength
# Filter the rows where the percentage exceeds 40%
# Uniq reads
allsplit2 = allsplit[(allsplit['T1_dif'] > int(args.AnchorLength)) | (allsplit['T2_dif'] > int(args.AnchorLength))]
allsplit2.to_csv(args.output,sep='\t',index=False)


del df_blast_1
del allsplit
del samstrand_box2_strand
del Nrm_all
del samstrand_box2
del df_blast_skipped_3
del normalsplitRead
del Nrm
del cmpl
del cmpl_all
del cmpl_inv
del st
del samstrand
del samstrand_box2_fd
del samstrand_box2_rv
del samstrand_box2_flip
del samstrand_box2_NoNaN
del samstrand_box2_NoNaN_subset

end_time = time.time()
usage = resource.getrusage(resource.RUSAGE_SELF)
max_memory = usage.ru_maxrss  # Memory in kilobytes

print(f"Execution Time: {end_time - start_time} seconds")
print(f"Max Memory Usage: {max_memory / 1024} MB")
