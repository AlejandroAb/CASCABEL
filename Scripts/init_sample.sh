#!/bin/bash
#needs 5 arguments
#$1 name of the project
#$2 name of the sample
#$3 metadata
#$4 fw reads o just reads
#$5 rv reads
#set -e
# Any subsequent(*) commands which fail will cause the shell script to exit immediately

if  [ $# -lt 5 ]; then
    echo -e "This program needs 5 arguments:"
    echo -e "\tArg1: Project name\n\tArg2: library name\n\tArg3: Full path to barcode file\n\tArg4: Full path to forward reads\n\tArg5: Full path to reverse reads"
   exit 1
fi

if [[ $3 == ..* ]]; then
    echo "Please specify full path to barcode file"
    exit 1
fi
if [[ $4 == ../* ]]; then
    echo "Please specify full path to forward reads"
    exit 1
fi
if [[ $5 == ../* ]]; then
    echo "Please specify full path to reverse reads"
    exit 1
fi

#-f for files

if [ ! -d "$1" ]; then
    mkdir $1
    echo "Project folder created..."
fi

cd $1

if [ ! -d samples ]; then
    mkdir samples
    echo "Samples folder created..."
fi

if [ ! -d metadata ]; then
    mkdir metadata
    echo "Barcode folder created..."
fi

cd metadata

if ln -fs  $3 sampleList_mergedBarcodes_$2.txt  ; then
    echo "Barcode list successfuly linked..."
else
    echo "Problems linking barcode list, make sure that file exists: "$3
    echo "Aborting...!"
    exit 1
fi

cd ..
cd samples

if [ ! -d $2 ]; then
    mkdir $2
    echo "Sample folder created..."
fi

cd $2

if [ ! -d "rawdata" ]; then
   mkdir rawdata
   echo "rawdata folder created..."
fi

#cd data
#mkdir rawdata
cd rawdata

if [[ $4 == *.gz ]] ; then
  if   ln -s $4 fw.fastq.gz ; then
     echo "Forward reads successfuly linked..."
  else
      echo "Problems linking forward reads, make sure that file exists: "$4
      #echo "Aborting...!"
      #exit 1
  fi
else
  if   ln -s $4 fw.fastq ; then
     echo "Forward reads successfuly linked..."
  else
      echo "Problems linking forward reads, make sure that file exists: "$4
      #echo "Aborting...!"
      #exit 1
  fi
fi

if [[ $5 == *.gz ]] ; then
  if   ln -s $5 rv.fastq.gz ; then
      echo "Reverse reads successfuly linked..."
  else
      echo "Problems linking reverse reads, make sure that file exists: "$5
      #echo "Aborting...!"
      #exit 1
  fi
else
  if   ln -s $5 rv.fastq ; then
      echo "Reverse reads successfuly linked..."
  else
      echo "Problems linking reverse reads, make sure that file exists: "$5
      #echo "Aborting...!"
      #exit 1
  fi
fi

echo "Sample "$2" structure successfuly created!"
