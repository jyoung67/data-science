#! /bin/bash
# call this script with one parameter; 

# 1 - directory to download SRA files to

mkdir -p $1

fastq-dump SRR1763770 --outdir $1
fastq-dump SRR1763771 --outdir $1
fastq-dump SRR1763772 --outdir $1
fastq-dump SRR1763773 --outdir $1
fastq-dump SRR1763774 --outdir $1
fastq-dump SRR1763775 --outdir $1
fastq-dump SRR1763776 --outdir $1
fastq-dump SRR1763777 --outdir $1
fastq-dump SRR1763778 --outdir $1
fastq-dump SRR1763779 --outdir $1
fastq-dump SRR1763780 --outdir $1
fastq-dump SRR1763781 --outdir $1

cd $1

mv SRR1763770.fastq BC_29.fastq
mv SRR1763771.fastq BC_20.fastq
mv SRR1763772.fastq BC_25.fastq
mv SRR1763773.fastq BC_21.fastq
mv SRR1763774.fastq BC_24.fastq
mv SRR1763775.fastq BC_17.fastq
mv SRR1763776.fastq BC_28.fastq
mv SRR1763777.fastq BC_23.fastq
mv SRR1763778.fastq BC_22.fastq
mv SRR1763779.fastq BC_19.fastq
mv SRR1763780.fastq BC_30.fastq
mv SRR1763781.fastq BC_26.fastq










    

