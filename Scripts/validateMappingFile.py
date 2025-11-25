##!/usr/bin/python
# Alejandro Abdala 
# Last recorded update: 05-2025
#

#"""
#This Script validates the mapping file.\n 
#Validation fields / rules
#  * Sample Name: Must be at the first column, start with a hash "#" and it could be called:
#   - SampleID (QIIME legacy)
#   - SampleName
#   - Sample
#   - ID
#   - For the moment there is no validation on the column content.
#  * Barcodes:  Does not matter the column number. All barcodes are thought to be from 5'to 3'.
#  * Captions on header names are not validated (not case sensitive). 
#    Valid names are:
#   - BarcodeSequence; Barcode (QIIME legacy)  only for Single end
#   - Forward_barcode; fwbarcode; FW_barcode; ForwardBarcode. Are valid names for the forward barcode on paired-end mode.
#   - Reverse_barcode; rvbarcode; RV_barcode; ReverseBarcode. Are valid names for the reverse barcode on paired-end mode.
#"""

import os
import subprocess
import shutil
import argparse
from argparse import RawTextHelpFormatter

#Global flags
DMX = False #demultiplex
PRM = False #checkprimers
LYT = 'PE' #Layout PE or SE


# Arguments you can use when running this script.
parser = argparse.ArgumentParser(description="""
This Script validates the mapping file.\n 
Validation fields / rules
  * Sample Name: Must be at the first column, start with a hash "#" and it could be called:
   - SampleID (QIIME legacy)
   - SampleName
   - Sample
   - ID
   - For the moment there is no validation on the column content.
  * Barcodes:  Does not matter the column number. All barcodes are thought to be from 5'to 3'.
  * Captions on header names are not validated (not case sensitive). 
    Valid names are:
   - BarcodeSequence; Barcode (QIIME legacy)  only for Single end
   - Forward_barcode; fwbarcode; FW_barcode; ForwardBarcode. Are valid names for the forward barcode on paired-end mode.
   - Reverse_barcode; rvbarcode; RV_barcode; ReverseBarcode. Are valid names for the reverse barcode on paired-end mode.
""",formatter_class=RawTextHelpFormatter)  #without RawTextHelpFormatter description ommit the new lines

#parser.add_argument('-m','--mapping-file', metavar='file', required = True, type=str, help='Mapping file to be validated.')
parser.add_argument('mapping_file', metavar='mapping_file',  type=str, help='Mapping file to be validated.')
parser.add_argument('--demultiplex-check', required = False,dest="DMX", action=argparse.BooleanOptionalAction,default=True, help='Whether this mapping file is going to be used for demultiplexing the library.')
#boolena parser: https://stackoverflow.com/questions/15008758/parsing-boolean-values-with-argparse
parser.add_argument('--primer-check', required = False,dest="PRM", type=str, action=argparse.BooleanOptionalAction,help='Whether this mapping file is going to be used for primer removal.')
parser.add_argument('-l','--layout', choices=['PE', 'SE'],dest="LYT", required = False, default='PE', type=str, help='Library layout for demultiplexing. Valid values are:\n "PE" for paired-end.\n "SE" for single-end.')
parser.add_argument('-o','--out', required = False, type=str, default='out.log', help='Output file with result of the validation')
parser.add_argument('-f','--forward-barcodes',  required = False, type=str, default='fw_barcodes.fasta', help='Output fasta file for forward barcodes or single-end barcodes')
parser.add_argument('-r','--reverse-barcodes',  required = False, type=str, default='rv_barcodes.fasta',help='Output fasta file for reverse barcodes')
#parser.add_argument('-s','--single-barcodes',  required = False, type=str, default='se_barcodes.fasta',help='Output fasta file for single-end barcodes')

args = parser.parse_args()
#For debuggin var names:
#print(vars(args))
# Check if the arguments required are supplied by the user (will do it automatically in the pipeline).
if args.mapping_file is None or not os.path.exists(args.mapping_file):
    print("Mapping file is missing!")
    exit(-1)

#Validation headers and init indexes. 
sampleHeader = ['SampleID','SampleName','Sample','ID']
idx_sample = -1 
barcodeHeader_se = ['BarcodeSequence', 'Barcode']
idx_barcode_se = -1
barcodeHeader_fw = ['fwbarcode', 'FW_barcode', 'ForwardBarcode', 'Forward_barcode']
idx_barcode_fw = -1
barcodeHeader_rv = ['rvbarcode', 'RV_barcode', 'ReverseBarcode', 'Reverse_barcode']
idx_barcode_rv = -1
primerHeader_fw = ['LinkerPrimerSequence','ForwardPrimerSequence','FW_primer','ForwardPrimer','FWPrimer']
idx_primer_fw = -1
primerHeader_rv = ['ReverseLinkerPrimerSequence','RvLinkerPrimerSequence', 'ReversePrimer','ReversePrimerSequence','RVPrimer']
idx_primer_rv = -1
fw_barcodes={}
rv_barcodes={}
se_barcodes={}
c_count=-1
def isDNAseq(str):
   validation = set(str)
   if validation.issubset(DNAchars):
      return True
   else:
      return False
      
DNAchars = set('ATGCXURYSWKMVDHVNIatgcxuryswkmvdhni')
with open(args.mapping_file, 'r') as mappingFile:
    l=0
    for line in mappingFile:
        l=l+1;
        columns = line.replace('#','').split('\t')
        #the header is always at row 1 and must contain these first 3 fields (qiime specs):
        #SampleID BarcodeSequence LinkerPrimerSequence Description
        if l==1 :
            c=0
            #Find target headers
            for col in columns:
                if col.casefold() in (name.casefold() for name in sampleHeader):
                    idx_sample=c
                elif col.casefold() in (name.casefold() for name in barcodeHeader_se):
                    idx_barcode_se=c
                elif col.casefold() in (name.casefold() for name in barcodeHeader_fw):
                    idx_barcode_fw=c
                elif col.casefold() in (name.casefold() for name in barcodeHeader_rv):
                    idx_barcode_rv=c
                elif col.casefold() in (name.casefold() for name in primerHeader_fw):
                    idx_primer_fw=c
                elif col.casefold() in (name.casefold() for name in primerHeader_rv):
                    idx_primer_rv=c
                
                c=c+1
            c_count=c #column count for the log
            if idx_sample == -1:
                print("\033[91m ERROR: Sample ID column is missing on the mapping file: "+ args.mapping_file +" \033[0m")
                exit(1)  
            if idx_sample != 0:
                print("\033[91m ERROR: Sample ID must be the at the first column \033[0m")
                exit(1)  
            if (args.DMX and args.LYT == 'PE') and  (idx_barcode_fw == -1 or idx_barcode_rv == -1):
                print ("\033[91m ERROR: For demultiplexing a paired-end library, you need to supply both forward and reverse barcode sequences.\033[0m")
                exit(1)
            if (args.PRM and args.LYT == 'PE') and  (idx_primer_fw == -1 or idx_primer_rv == -1):
                print ("\033[91m ERROR: Primer information is missing.\033[0m")
                exit(1)
                if idx_primer_fw == -1:
                    print ("\033[91m  - Forward primer not found.\033[0m")	
                if idx_primer_rv == -1:
                    print ("\033[91m  - Reverse primer not found.\033[0m")
            if args.PRM and args.LYT == 'SE' and  idx_primer_fw == -1 :
                print ("\033[91m ERROR: Primer information is missing.\033[0m")
                print ("\033[91m  - Forward primer not found.\033[0m")
        elif not line.startswith("#") and len(line) > 10:
            if args.DMX and args.LYT == 'PE':
                if not isDNAseq(columns[idx_barcode_fw]):
            	    print("\033[91m ERROR: Not valid DNA sequence for FORWARD BARCODE: "+ columns[idx_barcode_fw] + "\033[0m")
            	    exit(1)
                if not isDNAseq(columns[idx_barcode_rv]):
            	    print("\033[91m ERROR: Not valid DNA sequence for REVERSE BARCODE: "+ columns[idx_barcode_fw] + "\033[0m")
            	    exit(1)
                if columns[idx_sample] not in fw_barcodes:
                #check_key = fw_barcodes.setdefault(columns[idx_sample], columns[idx_barcode_fw])
                    fw_barcodes[columns[idx_sample]]=columns[idx_barcode_fw]
                    rv_barcodes[columns[idx_sample]]=columns[idx_barcode_rv]
                else:
	                print("\033[91m ERROR: Duplicated Sample ID: "+ columns[idx_sample] +". Line: "+ str(l) +"  \033[0m")
	                exit(1)
               
            if args.DMX and args.LYT == 'SE': 
                if not isDNAseq(columns[idx_barcode_se]):
            	    print("\033[91m ERROR: Not valid DNA sequence for BARCODE:"+ columns[idx_barcode_se] + "\033[0m")
            	    exit(1)               
                if columns[idx_sample] not in se_barcodes:
                    se_barcodes[columns[idx_sample]]=columns[idx_barcode_se]
                else:	
                    print("\033[91m ERROR: Duplicated Sample ID: "+ columns[idx_sample] +" \033[0m")
                    exit(1)    
            if args.PRM and args.LYT == 'PE': 
                if not isDNAseq(columns[idx_primer_fw]):
                    print("\033[91m ERROR: Not valid DNA sequence for FORWARD PRIMER:"+ columns[idx_primer_fw] + "\033[0m")
                    exit(1)
                if not isDNAseq(columns[idx_primer_rv]):
                    print("\033[91m ERROR: Not valid DNA sequence for REVERSE PRIMER:"+ columns[idx_primer_rv] + "\033[0m")
                    exit(1)
    
if args.DMX and args.LYT == 'PE':
    with open(args.forward_barcodes, 'w') as fw, open(args.reverse_barcodes, 'w') as rv:
        for sample in fw_barcodes:
            fw.write(">"+sample+"\n"+fw_barcodes[sample]+"\n")
            if sample in rv_barcodes:
                rv.write(">"+sample+"\n"+rv_barcodes[sample]+"\n")
            else:
                print("\033[91m ERROR: No reverse barcode for sample: "+ sample +" \033[0m")
                exit(1)

if args.DMX and args.LYT == 'SE':
    with open(args.forward_barcodes, 'w') as se:
        for sample, seq in se_barcodes:
            se.write(">"+sample+"\n"+seq+"\n")
            
with open(args.out, 'w') as out:
    out.write("Mapping file: " + args.mapping_file+"\n")
    out.write("Column count: " + str(c_count)+"\n")
    out.write("Sample ID column: " + str(idx_sample+1)+"\n")
    if args.DMX and args.LYT == 'PE':
        out.write("Number of samples: " + str(len(fw_barcodes))+"\n")	
        out.write("Barcodes check:"+"\n")
        out.write("  -Forward barcode column: " +str(idx_barcode_fw+1)+"\n")
        out.write("  -Reverse barcode column: " +str(idx_barcode_rv+1)+"\n")
    if args.DMX and args.LYT == 'SE':
        out.write("Number of samples: " + str(len(se_barcodes))+"\n")
        out.write("Barcodes check:"+"\n")
        out.write("  -barcode column: " +str(idx_barcode_se+1)+"\n")
    if args.PRM and args.LYT == 'PE':
        out.write("Primers check:"+"\n")
        out.write("  -Forward primer column: " +str(idx_primer_fw+1)+"\n")
        out.write("  -Reverse primer column: "+str(idx_primer_rv+1)+"\n")    	
