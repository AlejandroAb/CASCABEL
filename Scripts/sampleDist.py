import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

#df = pd.read_csv(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/seqs_fw_rev_filtered.dist.txt", sep="\t", header=None, names=['samples', 'counts'])
df = pd.read_csv(snakemake.input[0], sep="\t", header=None, names=['samples', 'counts'])
y_pos = np.arange(len(df['samples']))
plt.figure(figsize=(15,10))
plt.bar(y_pos, df['counts'])
plt.xticks(y_pos, df['samples'], color='blue', rotation='vertical')
plt.tick_params(axis='x', labelsize=7)
plt.savefig(snakemake.output[0])
#plt.savefig(snakemake.wildcards.PROJECT+"/runs/"+snakemake.wildcards.run+"/"+snakemake.wildcards.sample+"_data/seqs_fw_rev_filtered.dist.png")
