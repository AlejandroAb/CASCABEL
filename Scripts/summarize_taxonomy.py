import argparse
from argparse import RawTextHelpFormatter
import pandas as pd 
#import numphy as np 

parser = argparse.ArgumentParser(description="""
This Script summarize the taxonomy in a OTU table.\n 
""",formatter_class=RawTextHelpFormatter)  #without RawTextHelpFormatter description ommit the new lines

#parser.add_argument('-m','--mapping-file', metavar='file', required = True, type=str, help='Mapping file to be validated.')
parser.add_argument('otu_table', metavar='otu_table',  type=str, help='Mapping file to be validated.')
parser.add_argument('-l', '--levels', required = False,dest="levels", type=int,  default=6,help='Number of taxonomy levels to output (bottom-up)')
parser.add_argument('-o','--out', required = False, type=str, default='otus', help='Optional directory + prefix for the output file. i.e. /my_dir/otus')
parser.add_argument('-t', '--tax-delim', required = False,dest="tax_delim", type=str,  default=';',help='Taxonomy delimiter in OTU table')
parser.add_argument('-th', '--tax-header', required = False,dest="tax_header", type=str,  default='taxonomy',help='Name of the column containing')
parser.add_argument('--print-row-names', required = False,dest="print_names", action=argparse.BooleanOptionalAction,default=False, help='The row-names for the summarized tables consist of the concatenation of the summarized taxonomy.')

args = parser.parse_args()

otus=pd.read_csv(args.otu_table,sep="\t", header=0)
l=1
for l in range(args.levels+1):
	tlevels=[]
	for i in range(l):
		tlevels.append("L_"+str(i+1))
	if l>=1:
			#Fixed header as an attribute
			#otus[tlevels]=otus.taxonomy.str.split(args.tax_delim, n=l, expand=True).iloc[0:,0:l]
			try:
				otus[tlevels]=otus[args.tax_header].str.split(args.tax_delim, n=l, expand=True).iloc[0:,0:l]
			except ValueError as err:
				print(err)
				print("The number of taxonomy levels in the OTU table is less than the required number of levels to summarize level (-l): "+str(l))
			else:
				t=otus.groupby(tlevels,as_index=False).sum()
				subset=t.loc[:, ~t.columns.isin([args.tax_header, '#OTUID'])]
				if args.print_names:
					rownames=subset.loc[:, subset.columns.isin(tlevels)].apply(lambda x: ';'.join(x.dropna().astype(str).values), axis=1).values
					subset.index=rownames
					subset.to_csv(args.out+"_L"+str(i+1)+".txt", sep='\t', index=True)
				else:
					subset.to_csv(args.out+"_L"+str(i+1)+".txt", sep='\t', index=False)
#otu_table="test_otuTable.txt"
#otus=pd.read_csv(otu_table,sep="\t", header=0)
#tlevels=["0","1"]
#otus[tlevels]=otus.taxonomys.str.split(";", n=2, expand=True).iloc[0:,0:2]
#t=otus.groupby(tlevels,as_index=False).sum()
