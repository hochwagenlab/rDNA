##### In silico titration of reads #####
# June 2020                            #
# Daniel Sultanov                      #
########################################

#Data from Yue et al 2017, 10.1038/ng.3847
#S288c reads SRA accession - SRR4074255
#SK1 reads SRA accession - SRR4074258
#Raw reads are mapped to the S288C rDNA reference prototype beforehand (See the alignment strategy in "Alignment and variant call"):
#SRR4074255_rDNA.mapped.sam - mapped S288c rDNA reads
#SRR4074258_rDNA.mapped.sam - mapped SK1 rDNA reads
#(SK1 genome SRA accession - ERP014555; rDNA prototype - rDNA_S288c.fsa)

########################################

##load modules
module load samtools/intel/1.3.1

##Extract mapped reads from S288C
#exclude R1 & R2 are unmapped (BOTH) (note: if mapped with the "Alignment and variant call" script this is already included there). 12 = 4(read unmapped) + 8(mate unmapped)
#samtools view -h -F 12 <input file> > SRR4074255_rDNA.mapped.sam
#can check stats here with samtools flagstat
#convert to .bam with no header
samtools view -Sb SRR4074255_rDNA.mapped.sam > SRR4074255_rDNA.mapped.bam
#sorting by read name
samtools sort -n SRR4074255_rDNA.mapped.bam > SRR4074255_rDNA.mapped.readsort.bam
#extracting reads
samtools fastq -1 SRR4074255_rDNA.mapped.r1.fq -2 SRR4074255_rDNA.mapped.r2.fq -n SRR4074255_rDNA.mapped.readsort.bam
##remap the reads to the prototype and call variants using "Alignment and variant call". It is also a good sanity check
#the resultant variants matches the ones from the initian sam file. AF are virtually the same (they only differ by .00001 or smth like that which is negligible)

##Extract mapped reads from SK1
#if piping use -u instead of -b; it asks samtools to output an uncompressed BAM. Option `-u' is similar to `-b', but is preferred in piping because it saves time on compression/decompression.
samtools view -Sb SRR4074258_rDNA.mapped.sam > SRR4074258_rDNA.mapped.bam
samtools sort -n SRR4074258_rDNA.mapped.bam > SRR4074258_rDNA.mapped.readsort.bam
samtools fastq -1 SRR4074258_rDNA.mapped.r1.fq -2 SRR4074258_rDNA.mapped.r2.fq -n SRR4074258_rDNA.mapped.readsort.bam

## Titrate SK1:S288c reads
#SK1:S288C (%:%) 100:0 50:50 10:90 1:99 0.1:99.9 0.01:99.99
#number of lines in SRR4074255_rDNA.mapped.r1.fq is 2800144 (info: read name, read seq, + , quality)
#number of lines in SRR4074258_rDNA.mapped.r1.fq is 3214596
#to keep the same coverage, all the .fq files will have 2800144 lines (=length(s288c SRR4074255_rDNA.mapped.r*.fq files))

#SK1:S288C 100%:0%
head -n 2800144 SRR4074258_rDNA.mapped.r1.fq >> SK1_S288C_100_0.r1.fq
head -n 2800144 SRR4074258_rDNA.mapped.r2.fq >> SK1_S288C_100_0.r2.fq

#SK1:S288C 50%:50% Additional randomizing is not necessary since read lengths are already 'random' bc files are ordered just by read name
#for r1
head -n 1400072 SRR4074255_rDNA.mapped.r1.fq >> s288c.50.r1.fq
head -n 1400072 SRR4074258_rDNA.mapped.r1.fq >> sk1.50.r1.fq
cat sk1.50.r1.fq >> s288c.50.r1.fq
mv s288c.50.r1.fq SK1_S288C_50_50.r1.fq
#for r2
head -n 1400072 SRR4074255_rDNA.mapped.r2.fq >> s288c.50.r2.fq
head -n 1400072 SRR4074258_rDNA.mapped.r2.fq >> sk1.50.r2.fq
cat sk1.50.r2.fq >> s288c.50.r2.fq
mv s288c.50.r2.fq SK1_S288C_50_50.r2.fq

#SK1:S288C 10%:90%
head -n 2520128 SRR4074255_rDNA.mapped.r1.fq >> SK1_S288C_10_90.r1.fq
head -n 280016 SRR4074258_rDNA.mapped.r1.fq >> sk1.10.r1.fq
cat sk1.10.r1.fq >> SK1_S288C_10_90.r1.fq

head -n 2520128 SRR4074255_rDNA.mapped.r2.fq >> SK1_S288C_10_90.r2.fq
head -n 280016 SRR4074258_rDNA.mapped.r2.fq >> sk1.10.r2.fq
cat sk1.10.r2.fq >> SK1_S288C_10_90.r2.fq

#SK1:S288C 1%:99%
head -n 2772136 SRR4074255_rDNA.mapped.r1.fq >> SK1_S288C_1_99.r1.fq
head -n 28008 SK1/SRR4074258_rDNA.mapped.r1.fq >> sk1.1.r1.fq
cat sk1.1.r1.fq >> SK1_S288C_1_99.r1.fq

head -n 2772136 SRR4074255_rDNA.mapped.r2.fq >> SK1_S288C_1_99.r2.fq
head -n 28008 SRR4074258_rDNA.mapped.r2.fq >> sk1.1.r2.fq
cat sk1.1.r2.fq >> SK1_S288C_1_99.r2.fq

#SK1:S288C 0.5%:99.5%
head -n 2786140 SRR4074255_rDNA.mapped.r1.fq >> SK1_S288C_05_99.r1.fq
head -n 14004 SRR4074258_rDNA.mapped.r1.fq >> sk1.05.r1.fq
cat sk1.05.r1.fq >> SK1_S288C_05_99.r1.fq

head -n 2786140 SRR4074255_rDNA.mapped.r2.fq >> SK1_S288C_05_99.r2.fq
head -n 14004 SRR4074258_rDNA.mapped.r2.fq >> sk1.05.r2.fq
cat sk1.05.r2.fq >> SK1_S288C_05_99.r2.fq

#SK1:S288C 0.1%:99.9%
head -n 2797344 SRR4074255_rDNA.mapped.r1.fq >> SK1_S288C_01_99.r1.fq
head -n 2800 SRR4074258_rDNA.mapped.r1.fq >> sk1.01.r1.fq
cat sk1.01.r1.fq >> SK1_S288C_01_99.r1.fq

head -n 2797344 SRR4074255_rDNA.mapped.r2.fq >> SK1_S288C_01_99.r2.fq
head -n 2800 SRR4074258_rDNA.mapped.r2.fq >> sk1.01.r2.fq
cat sk1.01.r2.fq >> SK1_S288C_01_99.r2.fq

#SK1:S288C 0.05%:99.95%
head -n 2798744 SRR4074255_rDNA.mapped.r1.fq >> SK1_S288C_005_99.r1.fq
head -n 1400 SRR4074258_rDNA.mapped.r1.fq >> sk1.005.r1.fq
cat sk1.005.r1.fq >> SK1_S288C_005_99.r1.fq

head -n 2798744 SRR4074255_rDNA.mapped.r2.fq >> SK1_S288C_005_99.r2.fq
head -n 1400 SRR4074258_rDNA.mapped.r2.fq >> sk1.005.r2.fq
cat sk1.005.r2.fq >> SK1_S288C_005_99.r2.fq

#SK1:S288C 0.01%:99.99%
head -n 2799864 SRR4074255_rDNA.mapped.r1.fq >> SK1_S288C_001_99.r1.fq
head -n 280 SRR4074258_rDNA.mapped.r1.fq >> sk1.001.r1.fq
cat sk1.001.r1.fq >> SK1_S288C_001_99.r1.fq

head -n 2799864 SRR4074255_rDNA.mapped.r2.fq >> SK1_S288C_001_99.r2.fq
head -n 280 SRR4074258_rDNA.mapped.r2.fq >> sk1.001.r2.fq
cat sk1.001.r2.fq >> SK1_S288C_001_99.r2.fq
