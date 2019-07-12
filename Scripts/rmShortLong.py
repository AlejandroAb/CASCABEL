import os
from sys import stdin
minm=0
firstq=0
median=0
mean=0
thirdq=0
maxm=0
mode=0
histogram_txt="#reads  length\n"
with open(snakemake.input[2]) as hist_txt:
    for line in hist_txt:
        histogram_txt += line
    hist_txt.close()
with open(snakemake.input[0]) as hist:
    for line in hist:
        tmpLine = line.split(' ') #
        try:
            minm = float(tmpLine[0])
            firstq = float(tmpLine[1])
            median = float(tmpLine[2])
            mean = float(tmpLine[3])
            thirdq = float(tmpLine[4])
            maxm = float(tmpLine[5])
            mode = float(tmpLine[6])
            break
        except ValueError:
            print("Error trying to cast: "+ line)
    hist.close()
    if median > 0 and snakemake.config["interactive"] != "F":
        print("\033[91m This step can remove too short and too long reads \033[0m")
        print("\033[92m LIBRARY: "+snakemake.wildcards.sample+" \033[0m")
        print("\033[93m Sequence distribution: Min.       1st Qu.       Median       Mean       Mode       3rd Qu.       Max.  \033[0m")
        print("\033[92m                       "+str(minm) + "       "+str(firstq) + "       "+str(median) + "       "+str(mean) + "       "+str(mode) + "       "+str(thirdq) + "       "+str(maxm) + " \033[0m" )
        print("\033[93m You can see the histogram chart at: " + snakemake.input[1] + " \033[0m")
        print("\033[93m Please enter the option which fits better for your data: \033[0m")
        print("\033[93m 1. Use values from the configuration file: length < "+str(snakemake.config["rm_reads"]["shorts"])+" and length > "+str(snakemake.config["rm_reads"]["longs"])+ "\033[0m")
        print("\033[93m 2. Use values from median + /-"+str(snakemake.config["rm_reads"]["offset"])+": length < " + str(int(median)-snakemake.config["rm_reads"]["offset"]) + " and length > "+ str(int(median)+snakemake.config["rm_reads"]["offset"])  +" \033[0m")
        print("\033[93m 3. Specify new values! \033[0m")
        print("\033[93m 4. Print sequence length histogram \033[0m")
        print("\033[93m 5. Do not remove any sequence \033[0m")
        print("\033[93m 6. Interrupt workflow \033[0m")
        user_input="0"
        while (user_input != "1" and user_input !=  "2" and user_input != "3" and  user_input != "5" and user_input != "6"):
            if user_input == "4":
                print(histogram_txt)
            print("\033[92m Enter your option: \033[0m")
            user_input = stdin.readline() #READS A LINE
            user_input = user_input[:-1]
        if user_input == "1":
            shorts = snakemake.config["rm_reads"]["shorts"]
            longs = snakemake.config["rm_reads"]["longs"]
        elif user_input == "2":
            shorts = int(median)-snakemake.config["rm_reads"]["offset"]
            longs = int(median)+snakemake.config["rm_reads"]["offset"]
        elif user_input == "3":
            ss=-1
            while ss == -1:
                print("\033[92m Please enter the shortest length allowed: \033[0m")
                ui = stdin.readline() #READS A LINE
                ui = ui[:-1]
                try:
                    ss = int(ui)
                    shorts = ss
                except ValueError:
                    print ("Please enter a valid number")
                    ss = -1
            ll=-1
            while ll == -1:
                print("\033[92m Please enter the longest length allowed: \033[0m")
                ui = stdin.readline() #READS A LINE
                ui = ui[:-1]
                try:
                    ll = int(ui)
                    longs = ll
                except ValueError:
                    print ("Please enter a valid number")
                    ll = -1
        elif user_input == "5":
            shorts = 0
            longs = int(maxm) + 1
        elif user_input == "6":
            print("Aborting workflow...")
            exit(1)

        os.system("awk '!/^>/ { next } { getline seq } length(seq) > " + str(shorts) + " && length(seq) < " + str(longs) + " { print $0 \"\\n\" seq }' " + snakemake.input[3] + " > " + snakemake.output[0])
        with open(snakemake.output[1], "a") as tmplog:
            tmplog.write(snakemake.input[0] + "\t" + str(shorts) + "\t" + str(longs) + "\n")
            tmplog.close()
        #print("awk '!/^>/ { next } { getline seq } length(seq) > " + str(shorts) + " && length(seq) < " + str(longs) + " { print $0 \"\\n\" seq }' " + snakemake.input[0] + " > "+ snakemake.output[0])
        #os.system("awk '!/^>/ {{ next }} {{ getline seq }} length(seq) >= {config[rm_reads][shorts]} && length(seq) <= {config[rm_reads][longs]} {{ print $0 \"\\n\" seq }}' " + input[0] + " > {output}")

    elif median > 0 and snakemake.config["interactive"] == "F":
        if snakemake.config["rm_reads"]["non_interactive_behaviour"] == "AVG":
            shorts = int(median)-snakemake.config["rm_reads"]["offset"]
            longs = int(median)+snakemake.config["rm_reads"]["offset"]
        elif snakemake.config["rm_reads"]["non_interactive_behaviour"] == "CFG":
            shorts = snakemake.config["rm_reads"]["shorts"]
            longs = snakemake.config["rm_reads"]["longs"]
        elif snakemake.config["rm_reads"]["non_interactive_behaviour"] == "NONE":
            shorts = 0
            longs = int(maxm) + 1
        else:
            print("\033[91m" +"Invalid option for [rm_reads][non_interactive_behaviour] values at --configfile  \033[0m")
            print("\033[92m" +"Valid options are: AVG or GFG \033[0m")
            print("Aborting workflow...")
            exit(1)
        os.system("awk '!/^>/ { next } { getline seq } length(seq) > " + str(shorts) + " && length(seq) < " + str(longs) + " { print $0 \"\\n\" seq }' " + snakemake.input[3] + " > " + snakemake.output[0])
        print("\033[93m" +"Interactive mode off \033[0m")
        print("\033[92m LIBRARY: "+snakemake.wildcards.sample+" \033[0m")
        print("\033[93m Sequence distribution: Min.       1st Qu.       Median       Mean       Mode       3rd Qu.       Max.  \033[0m")
        print("\033[92m                       "+str(minm) + "       "+str(firstq) + "       "+str(median) + "       "+str(mean) + "       "+str(mode) + "       "+str(thirdq) + "       "+str(maxm) + " \033[0m" )
        print("\033[93m You can see the histogram chart at: " + snakemake.input[1] + " \033[0m")
        if snakemake.config["rm_reads"]["non_interactive_behaviour"] == "AVG":
            print("\033[93m" +"Removing sequences based on median ("+str(median)+") + / - "+str(snakemake.config["rm_reads"]["offset"])+": length > " + str(int(median)-snakemake.config["rm_reads"]["offset"]) + " and length < "+ str(int(median)+snakemake.config["rm_reads"]["offset"]) + "\033[0m")
            with open(snakemake.output[1], "a") as tmplog:
                tmplog.write("Interactive mode. remove short & long\n")
                tmplog.write(snakemake.input[0] + "\t" + str(shorts) + "\t" + str(longs) + "\n")
                tmplog.close()
        elif snakemake.config["rm_reads"]["non_interactive_behaviour"] == "NONE":
            print("\033[93mconfig value = NONE. Skipping length filtering...\033[0m")
            with open(snakemake.output[1], "a") as tmplog:
                tmplog.write("Interactive mode. remove short & long\n")
                tmplog.write(snakemake.input[0] + "\t0\tAll\n")
                tmplog.close()
        else:
            print("\033[93m" +"Removing sequences based on configuration file values: length > " + str(snakemake.config["rm_reads"]["shorts"]) + " and length < "+ str(snakemake.config["rm_reads"]["longs"]) + "\033[0m")
            with open(snakemake.output[1], "a") as tmplog:
                tmplog.write("Interactive mode. remove short & long\n")
                tmplog.write("Removing sequences based on configuration file values: length > " + str(snakemake.config["rm_reads"]["shorts"]) + " and length < "+ str(snakemake.config["rm_reads"]["longs"])+ "\n")
                tmplog.close()
