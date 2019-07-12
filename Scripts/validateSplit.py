from sys import stdin
import shutil
with open(snakemake.input[0]) as splitlog:
    for line in splitlog:
        if "Barcode not in mapping file:" in line:
            try:
                seqsNotBC = int(line[line.find(":")+1:])
                splitlog.close()
                break
            except ValueError:
                print("Error trying to cast: "+ line[line.find(":")+1:])
            #print("The flow will continue")
            #break
        #if (seqsNotBC < config["split"]["barcode_type"]):
if (seqsNotBC > 0):
    print("\033[91m" + "Barcode not in mapping file " + str(seqsNotBC) + "\033[0m")
    print("\033[91m" + "Sequences with those barcodes are dismiss \033[0m")
    print("Please take a look on complete log at: "+ snakemake.input[0])
    print("You can try to run split_libraries_fastq.py with --retain_unassigned_reads")
    print("\033[93m" +"Do you want to continue y/n?"+ "\033[0m")
    user_input = stdin.readline() #READS A LINE
    user_input = user_input[:-1]
    if user_input.upper() == "Y" or user_input.upper() == "YES":
        print("\033[92m" +"The flow goes on!"+ "\033[0m")
        with open(snakemake.output[0], "w") as tmplog:
            tmplog.write("Split warning dismissed, user continue...")
            tmplog.close()
    else:
        print("Aborting workflow...")
        print("Cleaning files...")
        shutil.rmtree(snakemake.params[0])
        shutil.rmtree(snakemake.params[1])
        exit(1)
else:
    with open(snakemake.output[0], "w") as tmplog:
        tmplog.write("Split validation log OK")
        tmplog.close()
