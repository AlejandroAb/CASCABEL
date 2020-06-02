import os

with open(snakemake.input[0]) as files:
    for library in files:
        tmpLine = library.split('\t') #
        try:
            lib = tmpLine[0]
            fw = tmpLine[1]
            rv = tmpLine[2]
            mapp = ""
            if len(tmpLine) > 3:
                mapp = tmpLine[3].rstrip()
            if lib.lower() == snakemake.wildcards.sample.lower():
                if len(mapp)>1 :
                    print("Scripts/init_sample.sh "+snakemake.wildcards.PROJECT+" " + lib +" "+mapp+" "+fw +" " +rv)
                    os.system("Scripts/init_sample.sh "+snakemake.wildcards.PROJECT+" " + lib +" "+mapp+" "+fw +" " +rv)
                else:
                    os.system("Scripts/init_sample_dmx.sh "+snakemake.wildcards.PROJECT+" " + lib +" "+fw +" " +rv)
                files.close()
                exit(0)
        except ValueError:
          print("Error trying to cast: "+ line)
    print("\033[92m There is no entry for LIBRARY: "+ snakemake.wildcards.sample + " in file: " + snakemake.input[0] + " \033[0m")
    print("\033[91m Exiting Cascabel \033[0m")




