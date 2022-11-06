####           Sequencing models for rDNA reads simultation             ####
# Original: Jun 2020; modified: Feb 2022                                   #
# Daniel Sultanov                                                          #
# - NEAT-genReads (see 10.1371/journal.pone.0167047)                       #
# - rDNA sequence: rDNA_repeat_S288c.fsa                                   #
# - S288c rDNA alignment SRR4074255.bam                                    #
# - S288c Illumina paired-end reads SRR4074255_1.fastq; SRR4074255_2.fastq #
############################################################################

#To create S288c (SRA: SRR4074255) rDNA alignment see "Alignment and variant call"

##Load modules
module load numpy/python2.7/intel/1.14.0
module load bedtools/intel/2.27.1
module load samtools/intel/1.9

##Calculate genome coverage
#.bam files must be sorted by position
bedtools genomecov -d -ibam SRR4074255.bam > SRR4074255.coverage

##Compute GC%
#Note: if it throws an error make sure the rDNA reference header has the same name  than strings in the bed coverage file (S288C);
#      i.e. in rDNA_repeat_S288c.fsa file, change ">S288C rDNA repeat coordinates 458433 to 467569" to ">S288C"
python neat-genreads-master/utilities/computeGC.py -r rDNA_repeat_S288c.fsa -i SRR4074255.coverage -o gcmodel.p

##Create fragment length model (creates fraglen.p in the working directory)
samtools view SRR4074255.bam | python neat-genreads-master/utilities/computeFraglen.py

#generate sequencing error model
python neat-genreads-master/utilities/genSeqErrorModel.py -i SRR4074255_1.fastq -i2 SRR4074255_2.fastq -o seq_error.p


