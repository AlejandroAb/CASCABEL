import os
import subprocess
import math
from functools import reduce


#https://stackoverflow.com/questions/11347505/what-are-some-approaches-to-outputting-a-python-data-structure-to-restructuredte
#where grid: [['Name', 'Favorite Food', 'Favorite Subject'],['Joe', 'Hamburgers', 'Cars'], ['Jill', 'Salads', 'American Idol'], ['Sally', 'Tofu', 'Math']]
def make_table(grid):
    cell_width = 2+max(reduce(lambda x,y: x+y, [[len(item) for item in row] for row in grid], []))
    num_cols = len(grid[0])
    rst = table_div(num_cols, cell_width, 0)
    header_flag = 1
    for row in grid:
        rst = rst + '| ' + '| '.join([normalize_cell(x, cell_width-1) for x in row]) + '|' +'\n'
        rst = rst + table_div(num_cols, cell_width, header_flag)
        header_flag = 0
    return rst

def table_div(num_cols, col_width, header_flag):
    if header_flag == 1:
        return num_cols*('+' + (col_width)*'=') + '+\n'
    else:
        return num_cols*('+' + (col_width)*'-') + '+\n'

def normalize_cell(string, length):
    return string + ((length - len(string)) * ' ')

def readBenchmark( benchFile ):
    try:
       with open(benchFile) as bfile:
           txt = "**Benchmark info:**\n\n"
           #prcAssignedOtus = "**" + "{:.2f}".format(prcAssigned) + "%**"
           fileData = []
           headers = []
           data =[]
           i=0
           for l in bfile:
               l.rstrip('\n')
               tmpLine = l.split('\t')
               #len(tmpLine -1)
               #txt += l + "\n\n"
               if i == 0:
                   y=0
                   for h in tmpLine:
                       if h != "h:m:s" and y != 1:
                           #headers[y] = h
                           headers.append(h.rstrip())
                       y+=1
               else:
                   y=0
                   for d in tmpLine:
                       if y == 0:#h:m:s
                           #data[y] = "{:.2f}".format(float(d))
                           data.append("{:.2f}".format(float(d)))
                       elif y != 1 :
                           #data[y] = d
                           data.append(d.rstrip())
                       y += 1
               i+=1

           #fileData[0] = headers
           #fileData[1] = data
           fileData.append(headers)
           fileData.append(data)
           txt += make_table(fileData)
           txt +="\n"
    except FileNotFoundError:
        print("Benchmark data not found: " + benchFile)
        txt = "Benchmark data not found: " + benchFile
    return txt

def readSampleDist( distFile, num_reads, num_samples ):
    cols = 1;
    if num_samples > 10:
        cols = 4;

    elements_per_tbl = math.ceil(num_samples/cols)
    try:
       with open(distFile) as bfile:
           #txt = "**Sample distribution:**\n\n"
           txt=""
           #prcAssignedOtus = "**" + "{:.2f}".format(prcAssigned) + "%**"
           fileData = []
           headers = []
           data =[]
           i=0
           headers.append("Sample")
           headers.append("Seqs")
           headers.append("prc.")
           fileData.append(headers)
           for l in bfile:
               l.rstrip('\n')
               tmpLine = l.split('\t')
               if len(tmpLine) > 1:
                   data.append(tmpLine[0].rstrip())
                   data.append(tmpLine[1].rstrip())
                   try:
                       prc = "{:.2f}".format((float(tmpLine[1].rstrip())/num_reads)*100)
                   except Exception as e:
                       prc = "ND"
                   data.append(prc)
                   fileData.append(data)
                   i+=1
                   data=[]
               if i  >= elements_per_tbl :
                   txt += make_table(fileData)
                   fileData = []
                   fileData.append(headers)
                   txt+="\n\n"
                   i=0
           if len(fileData) > 1:
               txt += make_table(fileData)
    except FileNotFoundError:
        print("Distribution file not found: " + distFile)
        txt = "Distribution file not found: " + distFile
    return txt

def countFasta(fastaFile, isFastq):
    if isFastq:
        reads = subprocess.run( ["cat " + fastaFile + " | wc -l"],stdout=subprocess.PIPE, shell=True)
    else:
        reads = subprocess.run( ["grep '^>' " + fastaFile + " | wc -l"],stdout=subprocess.PIPE, shell=True)
    totReads =  reads.stdout.decode('utf-8').strip()
    try:
        total = int(totReads)
        if isFastq:
            total = total/4
    except ValueError:
        print("Error counting entries: " + fastaFile)
        total = -1;

    return total

def countTxt(inFile):
    reads = subprocess.run( ["cat " + inFile + " | wc -l"],stdout=subprocess.PIPE, shell=True)
    totReads =  reads.stdout.decode('utf-8').strip()
    try:
        total = int(totReads)
    except ValueError:
        print("Error reading parsing file: " + inFile)
        total = -1;

    return total


def countFastaGZ(fastaFile, isFastq):
    if isFastq:
        reads = subprocess.run( ["zcat " + fastaFile + " | wc -l"],stdout=subprocess.PIPE, shell=True)
    else:
        reads = subprocess.run( ["zcat " + fastaFile + " | grep '^>' | wc -l"],stdout=subprocess.PIPE, shell=True)
    totReads =  reads.stdout.decode('utf-8').strip()
    try:
        total = int(totReads)
        if isFastq:
            total = total/4
    except ValueError:
        print("Error counting entries: " + fastaFile)
        total = -1;

    return total

def readSampleDist2( distFile, num_reads, num_samples ):
    elements_per_tbl = math.ceil(num_samples/4)
    try:
       with open(distFile) as bfile:
           txt = "**Sample distribution:**\n\n"
           #prcAssignedOtus = "**" + "{:.2f}".format(prcAssigned) + "%**"
           fileData = []
           headers = []
           data =[]
           i=0
           if num_samples > 10:
               headers.append("Sample")
               headers.append("Seqs")
               headers.append("prc.")
               headers.append("")
               headers.append("Sample")
               headers.append("Seqs")
               headers.append("prc.")
               headers.append("")
               headers.append("Sample")
               headers.append("Seqs")
               headers.append("prc.")
               headers.append("")
               headers.append("Sample")
               headers.append("Seqs")
               headers.append("prc.")
               headers.append("")
               fileData.append(headers)
           for l in bfile:
               l.rstrip('\n')
               tmpLine = l.split('\t')
               if len(tmpLine) > 1:
                   data.append(tmpLine[1].rstrip())
                   data.append(tmpLine[0].rstrip())
                   try:
                       prc = "{:.2f}".format((float(tmpLine[0].rstrip())/num_reads)*100)
                   except Exception as e:
                       prc = "ND"
                   data.append(prc)
                   fileData.append(data)
                   i+=1
                   data=[]
               if i  >= elements_per_tbl :
                   txt += make_table(fileData)
                   fileData = []
                   fileData.append(headers)
                   txt+="\n\n"
                   i=0
           if len(fileData) > 1:
               txt += make_table(fileData)
    except FileNotFoundError:
        print("Distribution file not found: " + distFile)
        txt = "Distribution file not found: " + distFile
    print(txt)
    return txt
