#!/bin/bash
#SBATCH --verbose
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8GB
#SBATCH --time=3:00:00
#SBATCH --job-name=simulation
#SBATCH --mail-type=END
#SBATCH --mail-user=ds5639@nyu.edu
#SBATCH --output=/scratch/ds5639/job11.out
#SBATCH --error=/scratch/ds5639/job11.err
######
#6/8/20
#simulation 11
#60000X rDNA coverage (==400X genomic coverage)
#load modules
module purge
module load numpy/python2.7/intel/1.14.0
module load samtools/intel/1.9
module load bowtie2/2.3.4.3
module load samtools/intel/1.3.1
module load lofreq_star/2.1.3.1 
#run simulator; AF= 0.005; control parameters: c, M, and o
python neat-genreads-master/genReads.py \
-r /scratch/ds5639/rDNA_project/Yue_2017_SK1_hiseq/SRR4074258_1/rDNA_S288C_ref/rDNA_S288C.fa \
-R 150 \
-o out9 \
--bam \
--vcf \
--pe-model fraglen.p \
-e seq_error.p \
--gc-model gcmodel.p \
-p 1 \
-M 0.015 \
-c 300 \
-t targeted_rDNA_2.bed \
-to 0.4 \
--rng 123

#explicitly specify -M 0 here
python neat-genreads-master/genReads.py \
-r /scratch/ds5639/rDNA_project/Yue_2017_SK1_hiseq/SRR4074258_1/rDNA_S288C_ref/rDNA_S288C.fa \
-R 150 \
-o out10 \
--bam \
--pe-model fraglen.p \
-e seq_error.p \
--gc-model gcmodel.p \
-p 1 \
-M 0 \
-c 59700 \
-t targeted_rDNA_2.bed \
-to 0.4 \
--rng 456

#merge jobs
neat-genreads-master/mergeJobs.py -i out9 out10 -o simulation11 -s /share/apps/samtools/1.9/intel/bin/samtools --no-job

#mapping
bowtie2 -5 1 -N 1 -p 8 -x /scratch/ds5639/rDNA_project/Yue_2017_SK1_hiseq/SRR4074258_1/rDNA_S288C_ref/rDNA_S288C.index -1 simulation11_read1.fq -2 simulation11_read2.fq -S simulation11_pipeline.sam
#convert to bam, sort, index
samtools view -Sbh -F 12 simulation11_pipeline.sam > simulation11_pipeline.bam
samtools sort -o simulation11_pipeline.sort.bam -O 'bam' -T 'temp' simulation11_pipeline.bam
rm simulation11_pipeline.bam
samtools index simulation11_pipeline.sort.bam
#call lofreq variant caller
lofreq indelqual --dindel -f /scratch/ds5639/rDNA_project/Yue_2017_SK1_hiseq/SRR4074258_1/rDNA_S288C_ref/rDNA_S288C.fa -o simulation11_pipeline.dindel.bam simulation11_pipeline.sort.bam
lofreq call --call-indels -f /scratch/ds5639/rDNA_project/Yue_2017_SK1_hiseq/SRR4074258_1/rDNA_S288C_ref/rDNA_S288C.fa -o simulation11_pipeline.vcf simulation11_pipeline.dindel.bam

#compare two .vcf files
python /scratch/ds5639/rDNA_project/Yue_2017_SK1_hiseq/SRR4074258_1/Liti_data/benchmark/simulation/neat-genreads-master/utilities/vcf_compare_OLD.py -r /scratch/ds5639/rDNA_project/Yue_2017_SK1_hiseq/SRR4074258_1/rDNA_S288C_ref/rDNA_S288C.fa -g simulation11_golden.vcf -w simulation11_pipeline.vcf -o simulation11 -a 0.002 --vcf-out --incl-fail --no-plot

mkdir simulation11
mv out9* simulation11
mv out10* simulation11
mv simulation11_* simulation11