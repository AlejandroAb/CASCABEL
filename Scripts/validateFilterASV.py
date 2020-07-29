import os
import subprocess
from sys import stdin
import shutil


sample_counts="#samples\treads.in\treads.out\n"

with open(snakemake.input[0]) as filter_summary:
    l=0
    summ=0
    samples=0
    for line in filter_summary:
        l+=1
        tmpLine = line.split('\t')
        if len(tmpLine) > 2:
          try:
              summ+=(float(tmpLine[2])/float(tmpLine[1]))*100
              sample_counts += line
              samples+=1
          except ValueError:
              summ+=0
    avg=float(summ/samples)
filter_summary.close()



if snakemake.config["interactive"] == "F":
    print("\033[93m" +"Interactive mode off \033[0m")
    print("\033[93m" + "Total number of samples: " + str(l) + "\033[0m")
    print("\033[93m" + "Average percentage of reads passing filters: " + "{0:.2f}".format(avg) + "% \033[0m")
    print("\033[93m" +"We suggest to review the filter log at: "+ snakemake.input[0]+ "\033[0m")
    with open(snakemake.output[0], "w") as tmplog:
        tmplog.write("Interactive mode off.\n")
        tmplog.write("Total number of samples: " + str(l)+ "\n")
        tmplog.write("Average percentage of reads passing filters: " + "{0:.2f}".format(avg)+"%")
        tmplog.close()
    exit(0)
else:

    if avg > 90:
        print("\033[92m" + "Total number of samples: " + str(samples) + "\033[0m")
        print("\033[92m" + "Average percentage of reads passing filters: " +"{0:.2f}".format(avg) +"% \033[0m")
        with open(snakemake.output[0], "w") as tmplog:
            tmplog.write( "Total number of samples: " + str(samples)+  "\n")
            tmplog.write( "Average percentage of reads passing filters: " + "{0:.2f}".format(avg)+"%")
            tmplog.close()
        print("\033[93mContinuing workflow... \033[0m")
        exit(0)
    else:
        print("\033[92m" + "Total number of samples: " + str(samples) + "\033[0m")
        print("\033[92m" + "Average percentage of reads passing filters: \033[91m" + "{0:.2f}".format(avg) + "% \033[0m")
        print("\033[92mPlease take a look into complete log file at: \033[93m "+ snakemake.input[0] + " \033[0m")
        #print("\033[92mFind the dada2 quality plots at: \033[93m "+ snakemake.input[0] + " \033[0m")
        print("\033[93m If too few reads are passing the filter, consider relaxing maxEE, \033[0m")
        print("\033[93m perhaps especially on the reverse reads, and reducing the truncLen  \033[0m")
        print("\033[93m to remove low quality tails. Remember though, when choosing truncLen  \033[0m")
        print("\033[93m for paired-end reads you must maintain overlap after truncation in  \033[0m")
        print("\033[93m order to merge them later.  \033[0m")
        print("\033[92m What would you like to do? \033[0m")
        print("\033[92m 1. Continue with the workflow \033[0m")
        print("\033[92m 2. Interrupt the workflow \033[0m")
        print("\033[92m 3. Print the number of reads \033[0m")
        print("\033[92m Enter your option: \033[0m")

        user_input = stdin.readline() #READS A LINE
        user_input = " ".join(user_input.split())
        while (user_input != "1" and user_input !=  "2"):
            if user_input == "3":
                print(sample_counts)
            print("\033[92m What would you like to do? \033[0m")
            print("\033[92m 1. Continue with the workflow \033[0m")
            print("\033[92m 2. Interrupt the workflow \033[0m")
            print("\033[92m 3. Print the number of reads \033[0m")
            print("\033[92m Enter your option: \033[0m")
            user_input = stdin.readline() #READS A LINE
            user_input = " ".join(user_input.split())
        if user_input == "1":
            print("\033[93mContinuing workflow... \033[0m")
            with open(snakemake.output[0], "w") as tmplog:
                tmplog.write( "Total number of samples: " + str(samples) + "\n")
                tmplog.write( "Average percentage of reads passing filters: " + "{0:.2f}".format(avg)+"%")
                tmplog.close()
            exit(0)
        elif user_input == "2":
            print("\033[91mAborting workflow... \033[0m")
            exit(1)



