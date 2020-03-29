from sys import stdin
if snakemake.config["interactive"] == "F":
    print("\033[93m" +"Interactive mode off \033[0m")
    print("\033[93m" +"We suggest to review the full Pear log at: "+ snakemake.input[0]+ "\033[0m")
    with open(snakemake.output[0], "w") as tmplog:
        tmplog.write("Interactive mode. Pear validation skipped")
        tmplog.close()
else:
    with open(snakemake.input[0]) as logfile:
        ok = False
        for line in logfile:
            if line.startswith("Assembled reads") and snakemake.config["UNPAIRED_DATA_PIPELINE"] != "T":
                try:
                    peared = float(line[line.find("(")+1:line.find("%")])
                except ValueError:
                    print("Error trying to cast: "+ line[line.find("(")+1:line.find("%")])

                if (peared < float(snakemake.config["pear"]["prcpear"])):
                    print("\033[91m" + "Peared percentage is not good enough ("+str(peared)+"%) Minimum expected: "+str(snakemake.config["pear"]["prcpear"])+"%\nMore info: " + snakemake.input[0] + "\033[0m")
                    #print("We suggest to try with different parameters or with the FLASH program (info on config.yaml)")
                    print("\033[93m" +"Do you want to continue anyway y/n?"+ "\033[0m")
                    user_input = stdin.readline() #READS A LINE
                    user_input = user_input[:-1]
                    if user_input.upper() == "Y" or user_input.upper() == "YES":
                        print("\033[92m" +"The flow goes on!"+ "\033[0m")
                        with open(snakemake.output[0], "w") as tmplog:
                            tmplog.write(str(peared)+"\n"+"User continue." )
                            tmplog.close()
                            logfile.close()
                            break
                    else:
                        print("Aborting workflow...")
                        logfile.close()
                        exit(1)
                else:
                    with open(snakemake.output[0], "w") as tmplog:
                        tmplog.write("Peared OK")
                        tmplog.close()
                    print("Pear OK: "+ str(peared))
                    break
            elif line.startswith("Not assembled reads") and snakemake.config["UNPAIRED_DATA_PIPELINE"] == "T":
                try:
                    peared = float(line[line.find("(")+1:line.find("%")])
                except ValueError:
                    print("Error trying to cast: "+ line[line.find("(")+1:line.find("%")])

                print("\033[92m****Un-assembled flow: Working with not assembled reads***\033[$0m")
                print("\033[91mUn-assembled percentage: "+str(peared)+"%\nMore info: " + str(snakemake.input[0]) + "\033[$0m")
                print("\033[93m" +"Do you want to continue?  y/n:"+ "\033[0m")
                user_input = stdin.readline() #READS A LINE
                user_input = user_input[:-1]
                if user_input.upper() == "Y" or user_input.upper() == "YES":
                    print("\033[92m" +"The flow goes on!"+ "\033[0m")
                    with open(snakemake.output[0], "w") as tmplog:
                        tmplog.write(str(peared)+"\n"+"User continue." )
                        tmplog.close()
                        logfile.close()
                        break
                else:
                    print("Aborting workflow...")
                    logfile.close()
                    exit(1)

