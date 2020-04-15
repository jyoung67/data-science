#! /bin/bash
# call this script with two parameters; 
# 1 - file path to directory
# 2 - path to adapter sequence filter

count=`ls -1 $1/*.trim.fastq 2>/dev/null | wc -l`

if [ $count == 0 ]
then
echo "Processing fastq files within directory $1."
cd $1

for f in *.fastq; do
    echo "File -> $f"
    b=`echo $f | sed 's/\.fastq$//'` 
    echo "Working on -> $b"
    TrimmomaticSE $b.fastq $b.trim.fastq ILLUMINACLIP:$2:4:10:3 SLIDINGWINDOW:30:20 MINLEN:25  
    rm $b.fastq
    
done



else
echo "Unable to process fastq files within direct $1, because they have already been processed."
fi 