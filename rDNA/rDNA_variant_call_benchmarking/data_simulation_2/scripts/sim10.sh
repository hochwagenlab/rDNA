#!/bin/bash
#SBATCH --verbose
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8GB
#SBATCH --time=2:00:00
#SBATCH --job-name=simulation
#SBATCH --mail-type=END
#SBATCH --mail-user=ds5639@nyu.edu
#SBATCH --output=/scratch/ds5639/job10.out
#SBATCH --error=/scratch/ds5639/job10.err
######
#6/8/20
#simulation 10
#15000X rDNA coverage (==100X genomic coverage)
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
-o out7 \
--bam \
--vcf \
--pe-model fraglen.p \
-e seq_error.p \
--gc-model gcmodel.p \
-p 1 \
-M 0.015 \
-c 75 \
-t targeted_rDNA_2.bed \
-to 0.4 \
--rng 123

#explicitly specify -M 0 here
python neat-genreads-master/genReads.py \
-r /scratch/ds5639/rDNA_project/Yue_2017_SK1_hiseq/SRR4074258_1/rDNA_S288C_ref/rDNA_S288C.fa \
-R 150 \
-o out8 \
--bam \
--pe-model fraglen.p \
-e seq_error.p \
--gc-model gcmodel.p \
-p 1 \
-M 0 \
-c 14925 \
-t targeted_rDNA_2.bed \
-to 0.4 \
--rng 456

#merge jobs
neat-genreads-master/mergeJobs.py -i out7 out8 -o simulation10 -s /share/apps/samtools/1.9/intel/bin/samtools --no-job

#mapping
bowtie2 -5 1 -N 1 -p 8 -x /scratch/ds5639/rDNA_project/Yue_2017_SK1_hiseq/SRR4074258_1/rDNA_S288C_ref/rDNA_S288C.index -1 simulation10_read1.fq -2 simulation10_read2.fq -S simulation10_pipeline.sam
#convert to bam, sort, index
samtools view -Sbh -F 12 simulation10_pipeline.sam > simulation10_pipeline.bam
samtools sort -o simulation10_pipeline.sort.bam -O 'bam' -T 'temp' simulation10_pipeline.bam
rm simulation10_pipeline.bam
samtools index simulation10_pipeline.sort.bam
#call lofreq variant caller
lofreq indelqual --dindel -f /scratch/ds5639/rDNA_project/Yue_2017_SK1_hiseq/SRR4074258_1/rDNA_S288C_ref/rDNA_S288C.fa -o simulation10_pipeline.dindel.bam simulation10_pipeline.sort.bam
lofreq call --call-indels -f /scratch/ds5639/rDNA_project/Yue_2017_SK1_hiseq/SRR4074258_1/rDNA_S288C_ref/rDNA_S288C.fa -o simulation10_pipeline.vcf simulation10_pipeline.dindel.bam

#compare two .vcf files
python /scratch/ds5639/rDNA_project/Yue_2017_SK1_hiseq/SRR4074258_1/Liti_data/benchmark/simulation/neat-genreads-master/utilities/vcf_compare_OLD.py -r /scratch/ds5639/rDNA_project/Yue_2017_SK1_hiseq/SRR4074258_1/rDNA_S288C_ref/rDNA_S288C.fa -g simulation10_golden.vcf -w simulation10_pipeline.vcf -o simulation10 -a 0.002 --vcf-out --incl-fail --no-plot

mkdir simulation10
mv out7* simulation10
mv out8* simulation10
mv simulation10_* simulation10