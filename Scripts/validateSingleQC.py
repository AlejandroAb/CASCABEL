from sys import stdin
fails=0
strFails=""
with open(snakemake.input[0]) as qc:
    for l in qc:
    #tmpLine = l.split('\t') #If we want to know the specific 'fails'
        if "FAIL" in l:
            fails+=1
            strFails+=l
        elif "WARN" in l:
            strFails+=l
    if fails > snakemake.config["qcLimit"]:
        #print("\x1b[6;30;42m" + "FastQC reports to many fails on raw file: " + input[i+4] + "\x1b[0m")
        print("\033[91m" + "FastQC reports too many fails on peared file: " + snakemake.input[2] + "\033[0m")
        print(strFails);
        print("We suggest to review the full FastQC report before continuing: "+ snakemake.input[1])
        print("\033[93m" +"Do you want to continue anyway y/n?"+ "\033[0m")
        user_input = stdin.readline() #READS A LINE
        user_input = user_input[:-1]
        #user_input = stdin.read(1)
        #if user_input == "Y" or user_input == "y":
        if user_input.upper() == "Y" or user_input.upper() == "YES":
            print("\033[92m" +"The flow goes on!"+ "\033[0m")
            with open(snakemake.output[0], "w") as tmplog:
                tmplog.write("Sequences are not best quality, user continue with the analysis")
                tmplog.close()
        else:
            print("Aborting workflow...")
            exit(1)
    else:
        with open(snakemake.output[0], "w") as tmplog:
            tmplog.write("Sequences pass fastQC")
            tmplog.close()
