#This scripts take the outfile from count_seqs.py and creates a rsText string
#with the expected format and information to be included in the final report

def parseCounts(ouput_counts):
    with open(ouput_counts) as counts:
        #n_calls = sum(1 for l in vcf if not l.startswith("#"))
        countTxt="Following you can see the final read counts: \n\n"
        fqVersion = ""
        for l in counts:
            tmpLine = l.split(' ')
            #print(tmpLine)
            if len(tmpLine)>1:
                filePath = tmpLine[3].split("/")
                fType = "undefined"
                if "raw" in tmpLine[3]:
                    fType = "raw data"
                elif "fw_rev" in tmpLine[3]:
                    fType = "accepted clean reads"
                elif "assembled" in tmpLine[3]:
                    fType = "Peared file"
                elif "split" in tmpLine[3]:
                    fType = "Split library file"
                elif "Total" in tmpLine[3]:
                    fType = "TOTAL"

                fileName = filePath[-1]
                if fType != "TOTAL":
                    countTxt += "* **File**: " + fileName + " **reads**: " + tmpLine[0] + " **type**: :green:`" + fType + "` \n\n"
    return countTxt
