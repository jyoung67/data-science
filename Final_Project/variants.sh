#! /bin/bash
# call this script with one parameter; reference genome optional

# 1 - reference genome basename
# 2 - base directory

cd $2

cp /mnt/c/Users/young/Documents/GitHub/Genomics/Lab05/NC_007898.fasta $2
cp /mnt/c/Users/young/Documents/GitHub/Genomics/Lab05/NC_007898.gff $2

bowtie2-build $1.fasta $1

for f in *.fastq; do
    echo "File -> $f"
    b=`echo $f | sed 's/\.fastq$//'` 
    echo "Working on -> $b"
    bowtie2 -x $1 -U $b.fastq -S $b.sam
    samtools view -uS $b.sam | samtools sort - -o $b.bam
    samtools index $b.bam
    samtools mpileup -uf $1.fasta $b.bam | bcftools call -c -v --ploidy 1 > $b.vcf
    /mnt/c/Users/young/Documents/UNCC_Spring_2020/Genomics/labs/Lab05/tools/IGV_2.8.0/IGV_2.8.0/igvtools index $b.vcf
    
done
