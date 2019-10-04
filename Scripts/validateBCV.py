import os
from sys import stdin
htmlFile="Not Found"
logFile="Not Found"
for file in os.listdir(snakemake.params[0]):
    if file.endswith(".log"):
        logFile=snakemake.params[0]+file
    elif file.endswith(".html"):
        htmlFile=snakemake.params[0]+file
user_bc_length = snakemake.config["ext_bc"]["bc_length"]
tot_length = 0;
if "--bc2_len" in user_bc_length:
    tmp_l = user_bc_length.index("--bc2_len")
    ll = int(user_bc_length[(tmp_l+10):])
    tmp_l1 = user_bc_length.index("bc1_len")
    ll1 = int(user_bc_length[(tmp_l1+8):tmp_l])
    tot_length = ll + ll1
elif "--bc1_len" in user_bc_length:
    tmp_l = user_bc_length.index("--bc1_len")
    ll = int(user_bc_length[(tmp_l+10):])
    tot_length = ll
else:
    print("\033[91m" + "Expected --bc1_len # at split:barcode_type into configuration file\033[0m")
  #bc_length: "--bc1_len 6 --bc2_len 6"
  #bc_length: "--bc1_len 12"
isWrong = False
message = "Barcode validation OK"
with open(snakemake.input[1]) as mappingFile:
    for line in mappingFile:
        #line.encode('utf-8').strip()
        if not line.startswith("#"):
            columns = line.split('\t')
            try:
                if(len(columns[1]) != tot_length):
                    print("\033[91m" + "The total length between ext_bc:bc_length and barcodes in mapping file differs!\033[0m")
                    print("\033[93m" + "ext_bc:bc_length:"+str(tot_length)+"\033[0m")
                    print("\033[93m" + "mapping barcode:"+str(len(columns[1]))+"\033[0m")
                    print("\033[92m" + "Please correct configuration file!\033[0m")
                    print("\033[91m" + "Aborting workflow!\033[0m")
                    exit(1)
                    break
                else:
                    print("\033[92m" + "Total length between extract_ba:bc_length and barcodes in mapping file: OK\033[0m")
                    break
            except IndexError:
                print("\033[91m" + "Error parsing file. We coul not validate barcode length\033[0m")
with open(logFile) as bcvlog:
    for line in bcvlog:
        if not "No errors or warnings found in mapping file" in line:
            print("\033[91m" + "Validation mapping file contains some warnings or errors: " + logFile + "\033[0m")
            print("Please take a look on complete report at: "+ htmlFile)
            print("\033[93m" +"If continue, maybe an error will be thrown during extract_bc rule. Do you want to continue anyway y/n?"+ "\033[0m")
            user_input = stdin.readline() #READS A LINE
            user_input = user_input[:-1]
            if user_input.upper() == "Y" or user_input.upper() == "YES":
                print("\033[92m" +"The flow goes on!"+ "\033[0m")
                with open(snakemake.output[0], "w") as tmplog:
                    tmplog.write("Error on barcode validation mapping, user continue...")
                    tmplog.close()
                break
            else:
                print("Aborting workflow...")
                logfile.close()
                exit(1)
        else:
            with open(snakemake.output[0], "w") as tmplog:
                tmplog.write("Barcode validation log OK")
                tmplog.close()
