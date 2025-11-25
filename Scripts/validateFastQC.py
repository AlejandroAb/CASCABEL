#!/usr/bin/python
# Alejandro Abdala 
# Last recorded update: 05-2025
#
"""
Script to validate fastqc txt output
"""

import argparse
from argparse import RawTextHelpFormatter
import sys

parser = argparse.ArgumentParser(description="""
This Script validates FastQC txt output.\n
If there are more errors than expected i.e. snakemake.config['fastQC']['qcLimit'])
Cascabel's interactive method stops the workflow
and waits for user validation.
This script only makes sense within Cascabel and therefore some variables/arguments still come
from snakemake built-in variables.
""",formatter_class=RawTextHelpFormatter)  #without RawTextHelpFormatter description ommit the new lines

#parser.add_argument('-m','--mapping-file', metavar='file', required = True, type=str, help='Mapping file to be validated.')
parser.add_argument('-f','--forward-log',  required = False, type=str, default='', help='FastQC "summary.txt" forward or single-end output file')
parser.add_argument('-r','--reverse-log',  required = False, type=str, default='', help='FastQC "summary.txt" reverse output file')
parser.add_argument('-fr','--forward-raw',  required = False, type=str, default='', help='Forward or single-end raw file')
parser.add_argument('-rr','--reverse-raw',  required = False, type=str, default='', help='Reverse raw file')
parser.add_argument('-t','--qc-threshold',  required = False, type=int, default=3, help='Maximum number of warning/errors before entering to interactive mode.')
parser.add_argument('--interactive', required = False,dest="interactive", action=argparse.BooleanOptionalAction,default=False, help='Whether if the script behaviour is interactive or not')

parser.add_argument('-o','--out',  required = False, type=str, default='', help='Output file')

args = parser.parse_args()

#from sys import stdin

def checkFastQlog(fqc_log,raw_file,out,append):
    if raw_file is None or not os.path.exists(args.forward_log):
        print("Reference for raw file is missing")
    fails=0
    strFails=""
    with open(fqc_log, 'r') as qc:
        for l in qc:
                #tmpLine = l.split('\t') #If we want to know the specific 'fails'
            if "FAIL" in l:
                fails+=1
                strFails+=l
            elif "WARN" in l:
                strFails+=l
    if fails > args.qc_threshold: # snakemake.config["fastQC"]["qcLimit"]:
            #print("\x1b[6;30;42m" + "FastQC reports to many fails on raw file: " + input[i+4] + "\x1b[0m")
        print("\033[91m" + "FastQC reports too many fails on raw data file: " + raw_file + "\033[0m")
        print(strFails);
        print("We suggest to review the full FastQC report before continuing")
        print("\033[93m" +"Do you want to continue anyway y/n?"+ "\033[0m")
        user_input = stdin.readline() #READS A LINE
        user_input = user_input[:-1]
        #user_input = stdin.read(1)
        #if user_input == "Y" or user_input == "y":
        if user_input.upper() == "Y" or user_input.upper() == "YES":
            print("\033[92m" +"The flow goes on!"+ "\033[0m")
            with open(out, append) as tmplog:
                tmplog.write("Sequences are not best quality, user continue with the analysis")
                tmplog.close()
        else:
            print("Aborting workflow...")
            exit(1)
    else:
        with open(out, append) as tmplog:
            tmplog.write("Sequences pass fastQC")



if  not args.interactive:
    print("\033[93m" +"Interactive mode off \033[0m")
    print("\033[93m" +"We suggest to review FastQC reports at: "+args.forward_log+" \033[0m")
    #for i in range(0,2):
    with open(args.out, "w") as tmplog:
            tmplog.write("Interactive mode. QC validation skipped")
            tmplog.close()
else:
    if args.forward_log is not None and os.path.exists(args.forward_log):
        checkFastQlog(args.forward_log, args.forward_raw, args.out, 'w')
    if args.reverse_log is not None and os.path.exists(args.reverse_log):
        checkFastQlog(args.reverse_log, args.reverse_raw, arg.out, 'a')