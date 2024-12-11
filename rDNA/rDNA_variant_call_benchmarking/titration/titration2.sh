#!/bin/bash
#SBATCH --verbose
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8GB
#SBATCH --time=2:00:00
#SBATCH --job-name=rDNA_SNPs
#SBATCH --mail-type=END
#SBATCH --mail-user=ds5639@nyu.edu
#SBATCH --output=/scratch/ds5639/job.out
#SBATCH --error=/scratch/ds5639/job.err
######
#load modules
module purge
module load bowtie2/2.3.4.3
module load samtools/intel/1.3.1
module load lofreq_star/2.1.3.1 
##map to rDNA
#SK2:S288C(%) 100:0 50:50 10:90 1:99 0.1:99.9 0.01:99.99
bowtie2 -5 1 -N 1 -p 8 -x /scratch/ds5639/rDNA_project/Yue_2017_SK1_hiseq/SRR4074258_1/rDNA_S288C_ref/rDNA_S288C.index -1 SK1_S288C_100_0.r1.fq -2 SK1_S288C_100_0.r2.fq -S SK1_S288C_100_0.sam
bowtie2 -5 1 -N 1 -p 8 -x /scratch/ds5639/rDNA_project/Yue_2017_SK1_hiseq/SRR4074258_1/rDNA_S288C_ref/rDNA_S288C.index -1 SK1_S288C_50_50.r1.fq -2 SK1_S288C_50_50.r2.fq -S SK1_S288C_50_50.sam
bowtie2 -5 1 -N 1 -p 8 -x /scratch/ds5639/rDNA_project/Yue_2017_SK1_hiseq/SRR4074258_1/rDNA_S288C_ref/rDNA_S288C.index -1 SK1_S288C_10_90.r1.fq -2 SK1_S288C_10_90.r2.fq -S SK1_S288C_10_90.sam
bowtie2 -5 1 -N 1 -p 8 -x /scratch/ds5639/rDNA_project/Yue_2017_SK1_hiseq/SRR4074258_1/rDNA_S288C_ref/rDNA_S288C.index -1 SK1_S288C_1_99.r1.fq -2 SK1_S288C_1_99.r2.fq -S SK1_S288C_1_99.sam
bowtie2 -5 1 -N 1 -p 8 -x /scratch/ds5639/rDNA_project/Yue_2017_SK1_hiseq/SRR4074258_1/rDNA_S288C_ref/rDNA_S288C.index -1 SK1_S288C_01_99.r1.fq -2 SK1_S288C_01_99.r2.fq -S SK1_S288C_01_99.sam
bowtie2 -5 1 -N 1 -p 8 -x /scratch/ds5639/rDNA_project/Yue_2017_SK1_hiseq/SRR4074258_1/rDNA_S288C_ref/rDNA_S288C.index -1 SK1_S288C_001_99.r1.fq -2 SK1_S288C_001_99.r2.fq -S SK1_S288C_001_99.sam
#convert to BAM, sort, index
samtools view -Sbh -F 12 SK1_S288C_100_0.sam > SK1_S288C_100_0.bam
samtools view -Sbh -F 12 SK1_S288C_50_50.sam > SK1_S288C_50_50.bam
samtools view -Sbh -F 12 SK1_S288C_10_90.sam > SK1_S288C_10_90.bam
samtools view -Sbh -F 12 SK1_S288C_1_99.sam > SK1_S288C_1_99.bam
samtools view -Sbh -F 12 SK1_S288C_01_99.sam > SK1_S288C_01_99.bam
samtools view -Sbh -F 12 SK1_S288C_001_99.sam > SK1_S288C_001_99.bam
samtools sort -o SK1_S288C_100_0.sort.bam -O 'bam' -T 'temp' SK1_S288C_100_0.bam
samtools sort -o SK1_S288C_50_50.sort.bam -O 'bam' -T 'temp' SK1_S288C_50_50.bam
samtools sort -o SK1_S288C_10_90.sort.bam -O 'bam' -T 'temp' SK1_S288C_10_90.bam
samtools sort -o SK1_S288C_1_99.sort.bam -O 'bam' -T 'temp' SK1_S288C_1_99.bam
samtools sort -o SK1_S288C_01_99.sort.bam -O 'bam' -T 'temp' SK1_S288C_01_99.bam
samtools sort -o SK1_S288C_001_99.sort.bam -O 'bam' -T 'temp' SK1_S288C_001_99.bam
samtools index SK1_S288C_100_0.sort.bam
samtools index SK1_S288C_50_50.sort.bam
samtools index SK1_S288C_10_90.sort.bam
samtools index SK1_S288C_1_99.sort.bam
samtools index SK1_S288C_01_99.sort.bam
samtools index SK1_S288C_001_99.sort.bam
#call lofreq variant caller
#first, add dindel quality (needs to call indels downstream; can skip if don't need it)
lofreq indelqual --dindel -f /scratch/ds5639/rDNA_project/Yue_2017_SK1_hiseq/SRR4074258_1/rDNA_S288C_ref/rDNA_S288C.fa -o SK1_S288C_100_0_dindel.sort.bam SK1_S288C_100_0.sort.bam
lofreq indelqual --dindel -f /scratch/ds5639/rDNA_project/Yue_2017_SK1_hiseq/SRR4074258_1/rDNA_S288C_ref/rDNA_S288C.fa -o SK1_S288C_50_50_dindel.sort.bam SK1_S288C_50_50.sort.bam
lofreq indelqual --dindel -f /scratch/ds5639/rDNA_project/Yue_2017_SK1_hiseq/SRR4074258_1/rDNA_S288C_ref/rDNA_S288C.fa -o SK1_S288C_10_90_dindel.sort.bam SK1_S288C_10_90.sort.bam
lofreq indelqual --dindel -f /scratch/ds5639/rDNA_project/Yue_2017_SK1_hiseq/SRR4074258_1/rDNA_S288C_ref/rDNA_S288C.fa -o SK1_S288C_1_99_dindel.sort.bam SK1_S288C_1_99.sort.bam
lofreq indelqual --dindel -f /scratch/ds5639/rDNA_project/Yue_2017_SK1_hiseq/SRR4074258_1/rDNA_S288C_ref/rDNA_S288C.fa -o SK1_S288C_01_99_dindel.sort.bam SK1_S288C_01_99.sort.bam
lofreq indelqual --dindel -f /scratch/ds5639/rDNA_project/Yue_2017_SK1_hiseq/SRR4074258_1/rDNA_S288C_ref/rDNA_S288C.fa -o SK1_S288C_001_99_dindel.sort.bam SK1_S288C_001_99.sort.bam
#call lofreq
lofreq call --call-indels -f /scratch/ds5639/rDNA_project/Yue_2017_SK1_hiseq/SRR4074258_1/rDNA_S288C_ref/rDNA_S288C.fa -o SK1_S288C_100_0.vcf SK1_S288C_100_0_dindel.sort.bam
lofreq call --call-indels -f /scratch/ds5639/rDNA_project/Yue_2017_SK1_hiseq/SRR4074258_1/rDNA_S288C_ref/rDNA_S288C.fa -o SK1_S288C_50_50.vcf SK1_S288C_50_50_dindel.sort.bam
lofreq call --call-indels -f /scratch/ds5639/rDNA_project/Yue_2017_SK1_hiseq/SRR4074258_1/rDNA_S288C_ref/rDNA_S288C.fa -o SK1_S288C_10_90.vcf SK1_S288C_10_90_dindel.sort.bam
lofreq call --call-indels -f /scratch/ds5639/rDNA_project/Yue_2017_SK1_hiseq/SRR4074258_1/rDNA_S288C_ref/rDNA_S288C.fa -o SK1_S288C_1_99.vcf SK1_S288C_1_99_dindel.sort.bam
lofreq call --call-indels -f /scratch/ds5639/rDNA_project/Yue_2017_SK1_hiseq/SRR4074258_1/rDNA_S288C_ref/rDNA_S288C.fa -o SK1_S288C_01_99.vcf SK1_S288C_01_99_dindel.sort.bam
lofreq call --call-indels -f /scratch/ds5639/rDNA_project/Yue_2017_SK1_hiseq/SRR4074258_1/rDNA_S288C_ref/rDNA_S288C.fa -o SK1_S288C_001_99.vcf SK1_S288C_001_99_dindel.sort.bam