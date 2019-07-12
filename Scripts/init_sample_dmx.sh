#!/bin/bash
#needs 3 arguments
#$1 name of the project
#$2 name of the sample
#$3 fw reads o just reads
#$4 rv reads
#set -e
# Any subsequent(*) commands which fail will cause the shell script to exit immediately

if  [ $# -lt 4 ]; then
    echo -e "This program needs 4 arguments:"
    echo -e "\tArg1: Project name\n\tArg2: Sample name\n\tArg3: Full path to forward reads\n\tArg4: Full path to reverse reads"
   exit 1
fi

if [[ $3 == ../* ]]; then
    echo "Please specify full path to forward reads"
    exit 1
fi
if [[ $4 == ../* ]]; then
    echo "Please specify full path to reverse reads"
    exit 1
fi

#-f para files

if [ ! -d "$1" ]; then
    mkdir $1
    echo "Project folder created..."
fi

cd $1

if [ ! -d samples ]; then
    mkdir samples
    echo "Samples folder created..."
fi

cd samples

if [ ! -d $2 ]; then
    mkdir $2
    echo "Sample folder created..."
fi

cd $2

if [ ! -d rawdata ]; then
    mkdir rawdata
    echo "Rawdata folder created..."
fi

cd rawdata

if [[ $3 == *.gz ]] ; then
  echo $3 " Ends with gz"
  if   ln -s $3 fw.fastq.gz ; then
     echo "Gun zipped forward reads successfuly linked..."
  else
      echo "Problems linking gun zipped forward reads, make sure that file exists: "$3
      #echo "Aborting...!"
      #exit 1
  fi
else
  echo $3 " Ends with fastq"
  if   ln -s $3 fw.fastq ; then
     echo "Forward reads successfuly linked..."
  else
      echo "Problems linking forward reads, make sure that file exists: "$3
      #echo "Aborting...!"
      #exit 1
  fi
fi

if [[ $4 == *.gz ]] ; then
  if   ln -s $4 rv.fastq.gz ; then
      echo "Reverse reads successfuly linked..."
  else
      echo "Problems linking reverse reads, make sure that file exists: "$4
      #echo "Aborting...!"
      #exit 1
  fi
else
  if   ln -s $4 rv.fastq ; then
      echo "Reverse reads successfuly linked..."
  else
      echo "Problems linking reverse reads, make sure that file exists: "$4
      #echo "Aborting...!"
      #exit 1
  fi
fi
echo "Sample "$2" structure successfuly created!"
