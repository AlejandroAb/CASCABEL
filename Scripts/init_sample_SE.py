import os

with open(snakemake.input[0]) as files:
    for library in files:
        if not library.startswith("#"):
            tmpLine = library.split('\t') #
            try:
                lib = tmpLine[0]
                fw = tmpLine[1]
                mapp = ""
                if len(tmpLine) > 2:
                    mapp = tmpLine[2].rstrip()
                if lib.lower() == snakemake.wildcards.sample.lower():
                    if len(mapp)>1 :
                        os.system("Scripts/init_sample_SE.sh "+snakemake.wildcards.PROJECT+" " + lib +" "+mapp+" "+fw)
                    else:
                        os.system("Scripts/init_sample_dmx_SE.sh "+snakemake.wildcards.PROJECT+" " + lib +" "+fw)
                    files.close()
                    exit(0)
            except ValueError:
              print("Error trying to cast: "+ line)
    print("\033[92m There is no entry for LIBRARY: "+ snakemake.wildcards.sample + " in file: " + snakemake.input[0] + " \033[0m")
    print("\033[91m Exiting Cascabel \033[0m")




