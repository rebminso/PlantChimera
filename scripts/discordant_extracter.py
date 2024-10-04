# ############################################################################################################################
#                                               Bedpe ReConstruction
# ############################################################################################################################


# discordant.bedpe
import pandas as pd
import re
import argparse
import resource
import time

start_time = time.time()

parser = argparse.ArgumentParser(description='Bedpe Reconstruction')
parser.add_argument('input',help='bedpe')
parser.add_argument('output',help='modified_bedpe')
args = parser.parse_args()


df_bedpe = pd.read_csv(args.input, sep="\t", header=None)
df_bedpe.columns = ["chr1", "start1", "end1", "chr2", "start2", "end2", "name", "score", "strand1", "strand2"]

df_bedpe = df_bedpe[df_bedpe["chr1"] != "."]
# df_bedpe[["chr1", "chr1_2"]] = df_bedpe[0].str.split(".", expand=True)
#df_bedpe[["chr1", "chr1_2"]] = df_bedpe["chr1"].str.split(".", expand=True)
#df_bedpe[["chr2", "chr2_2"]] = df_bedpe["chr2"].str.split(".", expand=True)
#df_bedpe_pt = df_bedpe[df_bedpe["chr1"] == "Pt"]
#df_bedpe_mt = df_bedpe[df_bedpe["chr1"] == "Mt"]
df_bedpe = df_bedpe[df_bedpe["chr1"] != df_bedpe["chr2"]]
#df_bedpe["chr1"] = df_bedpe["chr1"] + "." + df_bedpe["chr1_2"]
#df_bedpe["chr2"] = df_bedpe["chr2"] + "." + df_bedpe["chr2_2"]
#df_bedpe = df_bedpe.drop(["chr1_2", "chr2_2"], axis=1)
# df_bedpe
# df_bedpe

safile = df_bedpe[['name','chr1','chr2']]
safile = safile.rename(columns={'name': 'read_id', 'chr1': 'transcript1_id', 'chr2': 'transcript2_id'})
safile.to_csv(args.output,sep='\t',index=False)
# safile.shape

#clean up memory
del df_bedpe
#del df_bedpe_mt
#del df_bedpe_pt
del safile

end_time = time.time()
usage = resource.getrusage(resource.RUSAGE_SELF)
max_memory = usage.ru_maxrss  # Memory in kilobytes

print(f"Execution Time: {end_time - start_time} seconds")
print(f"Max Memory Usage: {max_memory / 1024} MB")
