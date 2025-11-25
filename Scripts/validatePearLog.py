import os
import subprocess
import shutil
import gzip
import argparse
from argparse import RawTextHelpFormatter
from sys import stdin

parser = argparse.ArgumentParser(description="""
This Script validates validates the read assembly process.
It count all the *log files on the (argument) directory
""",formatter_class=RawTextHelpFormatter)  #without RawTextHelpFormatter description ommit the new lines

#parser.add_argument('-m','--mapping-file', metavar='file', required = True, type=str, help='Mapping file to be validated.')
parser.add_argument('log_directory', metavar='log_directory',  type=str, help='Directory containing log files from PEAR')
#parser.add_argument('-r','--raw-read', dest="raw_data", required = True, type=str, help='Raw file used for the demultiplexing process (it could be R1 or R2).')
#parser.add_argument('--is-raw-', required = False,dest="interactive", action=argparse.BooleanOptionalAction,default=True, help='Whether if the script behaviour is interactive or not')
#parser.add_argument('-o','--out', dest="out", required = False, default='summary_demux.txt', type=str, help='Output file with demultiplexing results per sample')
parser.add_argument('--interactive', required = False,dest="interactive", action=argparse.BooleanOptionalAction,default=False, help='Whether if the script behaviour is interactive or not')
parser.add_argument('-o','--out', dest="out", required = False, default='pear_validation.out', type=str, help='Output file with demultiplexing results per sample')
parser.add_argument('-l','--log', dest="log", required = False, default='pear_validation.tsv', type=str, help='Log file of this script')
parser.add_argument('-m','--min', dest="min", required = False, default='80', type=float, help='Minimum percentage needed to bypass the warning. Default 80')

#parser.add_argument('-m','--min', dest="min", required = False, default='80', type=str, help='')


args = parser.parse_args()


tot_assembled = 0
tot_discarded = 0
tot_not_assembled = 0
total = 0
sample_log = "Sample\tReads\tAssembled\tDiscarded\tNot_Assembled\n"
with open(args.out, 'w') as out:
    out.write("Sample\tReads\tAssembled\tDiscarded\tNot_Assembled\n")
    for file in os.listdir(args.log_directory):
        sample=file[:-4]
        if file.endswith(".log") and sample != "pear":
            with open(args.log_directory+"/"+file, 'r') as logfile:            
                assembled = 0
                tot_sample = 0
                discarded = 0
                not_assembled =0
                for line in logfile:
                    if line.startswith("Assembled reads ...") :#and snakemake.config["UNPAIRED_DATA_PIPELINE"] != "T":
                        l =line.split(" ")
                        try:
                            assembled = int(l[3].replace(",",""))
                            tot_assembled = tot_assembled + assembled
                            tot_sample = int(l[5].replace(",",""))
                            total = total + tot_sample  
                        except ValueError:
                            print("Error trying to cast: "+ line)
                    if line.startswith("Discarded reads ...") :
                        l =line.split(" ")
                        try:
                            discarded = int(l[3].replace(",",""))
                            tot_discarded = tot_discarded + discarded
                        except ValueError:
                            print("Error trying to cast: "+ line)
                    if line.startswith("Not assembled reads ...") :#and snakemake.config["UNPAIRED_DATA_PIPELINE"] != "T":
                        l =line.split(" ")
                        try:
                            not_assembled = int(l[4].replace(",",""))
                            tot_not_assembled = tot_not_assembled + not_assembled
                        except ValueError:
                            print("Error trying to cast: "+ line)
                           
                if sample != "pear":
                    out.write(sample+"\t"+str(tot_sample)+"\t"+str(assembled)+"\t"+str(discarded)+"\t"+str(not_assembled)+"\n")
                    sample_log = sample_log + sample+"\t"+str(tot_sample)+"\t"+str(assembled)+"\t"+str(discarded)+"\t"+str(not_assembled)+"\n" 
    prc=(tot_assembled/total)*100
    prc_not=(tot_not_assembled/total)*100
    prc_dis=(tot_discarded/total)*100
    #"{0:.2f}".format(avg)
    print("***Read assembly validation***")
    print("* Total reads:"+str(total))
    print("* Assembled reads:"+ str(tot_assembled) + " ("+"{0:.2f}".format(prc)+"%)") 
    print("* Not assembled reads: "+ str(tot_not_assembled)+ " ("+"{0:.2f}".format(prc_not)+"%)")
    print("* Discarded reads: "+ str(tot_discarded)+ " ("+"{0:.2f}".format(prc_dis)+"%)")
    with open(args.log, "w") as tmplog:
        tmplog.write("***Read assembly validation***\n")
        tmplog.write("* Total reads: "+str(total)+"\n")
        tmplog.write("* Assembled reads: "+ str(tot_assembled) + " ("+"{0:.2f}".format(prc)+"%)\n") 
        tmplog.write("* Not assembled reads: "+ str(tot_not_assembled)+ " ("+"{0:.2f}".format(prc_not)+"%)\n")
        tmplog.write("* Discarded reads: "+ str(tot_discarded)+ " ("+"{0:.2f}".format(prc_dis)+"%)\n")
        tmplog.close()    
    if not args.interactive:
        print("\033[93m" +"Interactive mode off \033[0m")
        print("\033[93m" +"We suggest to review Pear log files at: "+ args.log_directory+ "\033[0m")

    elif (prc < args.min):
    	go_on=True
    	while go_on:
            print("\033[91m" + "Peared percentage is not good enough ("+"{0:.2f}".format(prc)+"%)!\nMinimum expected: "+str(args.min)+"%\nMore info: " + args.log + "\033[0m")
            print("\033[91mPress 1 to see the assembly reads summary\033[0m")
            #print("We suggest to try with different parameters or with the FLASH program (info on config.yaml)")
            print("\033[93m" +"Do you want to continue anyway y/n?"+ "\033[0m")
            user_input = stdin.readline() #READS A LINE
            user_input = user_input[:-1]
            if user_input.upper() == "Y" or user_input.upper() == "YES":
                print("\033[92m" +"The flow goes on!"+ "\033[0m")
                with open(args.log, "a") as tmplog:
                    tmplog.write("User continue." )
                go_on=False
            elif user_input == "1":
                print(sample_log)
            elif user_input.upper() == "N" or user_input.upper() == "NO":
                print("Aborting workflow...")
                exit(1)

    else:
        with open(args.log, "a") as tmplog:
            tmplog.write("Peared OK")
            
   
