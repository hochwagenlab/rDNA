##### Alignment and variant call########
# June 2020                            #
# Daniel Sultanov                      #
#                                      #
# set the parameters to run on cluster #
# SCRIPT STARTS BELOW                  #
########################################

#!/bin/bash
#SBATCH --job-name=
#SBATCH --output=
#SBATCH --error=
#SBATCH --verbose
#SBATCH --array=
#SBATCH --time=
#SBATCH --nodes=
#SBATCH --cpus-per-task=
#SBATCH --mem=
#SBATCH --mail-type=

##Create folders to store results
#make sure folders with these names do not exist prior
#mkdir raw_vcf
#mkdir filtered_vcf_plain
#mkdir filtered_vcf
#mkdir rDNA_coverage

##load modules
#check for current version and change accordingly
module purge
module load sra-tools/2.10.5
module load bowtie2/2.3.4.3
module load samtools/intel/1.3.1
module load lofreq_star/2.1.3.1
module load bedtools/intel/2.27.1

##set variables
#SraAccList_indexed.txt is tab-delimeted indexed list of accession names, i.e.: 1    ERR123456
#caution: there should not be extra spaces between '='
sample=$(awk -v taskID=$SLURM_ARRAY_TASK_ID '$1==taskID {print $2}' SraAccList_indexed.txt)

##download fastq files
fastq-dump -I --split-files $sample

##map
#strict mapping (need UNmapped reads)
#the SK1 genome used in this script is from Yue et al. Change /path/to/genome/index/prefix accordingly
#it DOES NOT contain annotated rDNA
bowtie2 -5 1 -N 0 -p 8 \
--un-conc ${sample}\_strictunmap.fastq \
-x /path/to/genome/index/basename \
-1 ${sample}\_1.fastq \
-2 ${sample}\_2.fastq \
-S ${sample}\_strictmap.sam

#relaxed mapping
bowtie2 -5 1 -N 1 -p 8 \
--un-conc ${sample}\_rDNA_rich.fastq \
-x /path/to/genome/index/basename \
-1 ${sample}\_strictunmap.1.fastq \
-2 ${sample}\_strictunmap.2.fastq \
-S ${sample}\_loosemap.sam

#map these reads to the rDNA prototype; with relaxed settings. Change /path/to/rDNA/prototype/index/prefix accordingly
bowtie2 -5 1 -N 1 -p 8 \
-x /path/to/rDNA/prototype/index/basename \
-1 ${sample}\_rDNA_rich.1.fastq \
-2 ${sample}\_rDNA_rich.2.fastq \
-S ${sample}\_rDNA.sam

#convert to BAM, sort
samtools view -Sbh -F 12 ${sample}\_rDNA.sam > ${sample}\_rDNA.bam
samtools sort -@ 8 -o ${sample}\_rDNA.sort.bam -O 'bam' ${sample}\_rDNA.bam

#remove HERE the initial bam file
rm ${sample}\_rDNA.bam

#add dindel. Change /path/to/rDNA_prototype_fasta_file accordingly
lofreq indelqual --dindel \
-f /path/to/rDNA_prototype_fasta_file \
-o ${sample}\_rDNA.bam \
${sample}\_rDNA.sort.bam

##calculate coverage
#NB: -g is ignored when BAM input is provided so do not need to specify it here actually
#Change /path/to/rDNA_prototype_fasta_file accordingly
bedtools genomecov -d \
-ibam ${sample}\_rDNA.bam \
-g path/to/rDNA_prototype_fasta_file > ${sample}\_rDNA_coverage.txt

##call variants
#Change /path/to/rDNA_prototype_fasta_file accordingly
lofreq call --call-indels -f path/to/rDNA_prototype_fasta_file \
-o ${sample}\_rDNA.vcf \
${sample}\_rDNA.bam

##process and filter .vcf
#NOTE: this does not filter from: 1) SB and 2) variants embedded in homopolymer regions. You can do this later using R or other language/program
#create column names for the future processed vcf file (including the 'ERR' column)
awk 'BEGIN {printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "ACCESSION","CHROM","POS","ID","REF","ALT","QUAL","FILTER","DP","AF","SB","DP4","INDEL","HRUN")}' > ${sample}\_filtered_vcf.txt

#reformat .vcf
awk '/^[^#]/ {print $0}' ${sample}\_rDNA.vcf | awk -F ";" '$1=$1' OFS="\t" > ${sample}\_temp1.txt
awk -F "\t" 'BEGIN {s = 8; e = 13} {for (i=s; i<=e; i++) sub(".*=","",$i) ; print}' OFS="\t" ${sample}\_temp1.txt > ${sample}\_temp2.txt

#subset POS>10 and POS<9100; also for HRUN<4 and AF>0.005
awk '$2 > 10 && $2 < 9100 && $13 < 4 && $9 > 0.005 {print}' ${sample}\_temp2.txt > ${sample}\_temp3.txt

#filter out long (>5 nt) REF and ALT for GC content (>0.6)
awk '{b = length($4); if (b<5) print $0; else if ((gsub(/G/,"H",$4)+gsub(/C/,"D",$4))/b<0.55) print $0; else print("")}' ${sample}\_temp3.txt > ${sample}\_temp4.txt
awk '{b = length($5); if (b<5) print $0; else if ((gsub(/G/,"H",$5)+gsub(/C/,"D",$5))/b<0.55) print $0; else print("")}' ${sample}\_temp4.txt > ${sample}\_temp5.txt
awk  'gsub(/H/,"G",$4); $1=$1' OFS="\t" ${sample}\_temp5.txt | awk 'gsub(/D/,"C",$4); $1=$1' OFS="\t" | awk 'gsub(/H/,"G",$5); $1=$1' OFS="\t" | awk 'gsub(/D/,"C",$5); $1=$1' OFS="\t" | awk '!seen[$0]++' > ${sample}\_temp6.txt

#add accession name as one of the columns. This file does not contain column names
awk -F="\t" -v isolate=$sample '{print isolate"\t"$0}' ${sample}\_temp6.txt > ${sample}\_filtered_vcf_plain.txt

#same with column names
cat ${sample}\_filtered_vcf_plain.txt >> ${sample}\_filtered_vcf.txt

##move files
mv ${sample}\_rDNA.vcf raw_vcf
mv ${sample}\_filtered_vcf_plain.txt filtered_vcf_plain
mv ${sample}\_filtered_vcf.txt filtered_vcf
mv ${sample}\_rDNA_coverage.txt rDNA_coverage

#remove unnecessary files
rm ${sample}\_strictmap.sam
rm ${sample}\_*unmap*
rm ${sample}\_loosemap.sam
rm ${sample}\_strictunmap.*.fastq
rm ${sample}\_rDNA_rich.*.fastq
rm ${sample}\_rDNA.sort.bam
rm ${sample}\_temp*.txt
rm ${sample}\_*fastq
rm ${sample}\_rDNA.sam
