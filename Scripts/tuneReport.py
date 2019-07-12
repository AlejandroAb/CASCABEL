#import pdfkit
with open(snakemake.input[0]) as reportfile:
    with open(snakemake.output[0], "w") as newreportfile:
        newReport = ""
        findSD=False
        tables=0;
        for line in reportfile:
            #if "Sample distribution" in line:
            #    findSD = True
            if "id=\"sample-distribution\"" in line:
                #newReport = line.replace("\>"," style=\"float: left; width:100%; margin-right: 5px; \"\>",1)
                newReport = "<div class=\"section\" id=\"sample-distribution\" style=\"float: left; width:100%; margin-right: 5px; \">"
                findSD = True
            elif line.startswith("body {"):
                newReport = "p.cmmd{ text-align: left; padding: 10px; border-style:solid; border-color:#99AAC7; max-width:95%; margin-left:2%; height: auto; background-color: #010101; color: white; word-wrap:normal; }\n commd{text-align: left;} \n span.red{color:red;}\nspan.green{color:#008800;}\n"
                newReport += ".zui-table {border: solid 1px #DDDDDD; border-collapse: collapse; border-spacing: 0; font: normal 12px Arial, sans-serif;} .zui-table thead th { background-color: #EFEFEF; border: solid 1px #DDEEEE; color: #336B6B; padding: 10px; text-align: left; text-shadow: 1px 1px 1px #fff;} .zui-table tbody td { border: solid 1px #DDEEEE; color: #333; padding: 10px; text-shadow: 1px 1px 1px #fff; }\n"
                newReport += line
            elif line.startswith("div#metadata {"):
                newReport = "div.document p.cmmd{ text-align: left;}\n"
                newReport += line
            elif "class=\"commd\"" in line:
                newReport = line.replace("<p>","<p class=\"cmmd\">",1)
            elif "class=\"docutils\"" in line and not findSD:
                newReport = line.replace("docutils","zui-table",1)
            elif "class=\"docutils\"" in line and findSD:#this is for the sample distribution table
                #tables+=1
                #if tables < 4:
                newReport = line.replace("\"docutils\"","\"zui-table\"  style=\"float: left; margin-right: 5px; \"",1)
                #else:
                #    newReport = line.replace("\"docutils\"","\"zui-table\"  style=\"float: right; margin-right: 5px; \"",1)
            elif "<colgroup" in line or "<col width" in line or "</colgroup" in line:
                newReport = "" #skip print those lines
            else :
                newReport = line
            newreportfile.write(newReport)
        reportfile.close()
        newreportfile.close()
#    pdfkit.from_file(snakemake.input[0], 'out.pdf')
