import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

#distribution="NIOZ66_bc_correction/runs/bc_correction/NIOZ66_data/seqs_fw_rev_filtered.dist.txt"
#output="NIOZ66_bc_correction/runs/bc_correction/report_files/seqs_fw_rev_filtered.NIOZ66.dist.png"
distribution=sys.argv[1]
output=sys.argv[2]

df = pd.read_csv(distribution, sep="\t", header=None, names=['samples', 'counts'])
#df = pd.read_csv(snakemake.input[0], sep="\t", header=None, names=['samples', 'counts'])
y_pos = np.arange(len(df['samples']))
plt.figure(figsize=(15,10))
plt.bar(y_pos, df['counts'])
plt.xticks(y_pos, df['samples'], color='blue', rotation='vertical')
plt.tick_params(axis='x', labelsize=7)
plt.savefig(output)

