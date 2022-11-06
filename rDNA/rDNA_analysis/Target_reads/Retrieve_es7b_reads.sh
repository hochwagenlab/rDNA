####                     Retrieve reads for expansion segment 7b (25S)                    ####
# Original: Mar 2021; modified: Feb 2022                                                     #
# Daniel Sultanov                                                                            #
# Prerequisites:                                                                             #
# -mapped rDNA reads from ERR1309511 (ERR1309511_rDNA.bam; see "Alignment and variant call") #
# -annotation files for region of interest (.bed; see below how to create)                   #
##############################################################################################

######################################### ES7b in 25S #########################################
#25S rRNA coordinates: 524-568
#From pairwise analysis I found 5 cooccurring varinats in 7 industrial isolatesin ES7b;
#6 of the variants have VF > 99%, and ERR1309511 has VF ~28%
#Check if all the variants are present in the same set of rDNA copies (in the same short reads)
#Note: this script works fine but there are probably more efficient and scalable way to do this
###############################################################################################

#################### Make annotation .bed files ####################
#format: CHROM  START       END
#        S288C  0-based     non-inclusive

#for the whole ES7b rDNA coordinates (1-base inclusive): 5319-5363 :
#es7b.bed
#S288C   5318  5363

#for region within ES7b with each variant position separately (1 nt):
#es5322.bed
#S288C  5321    5322
#es5338.bed
#S288C  5337    5338
#es5343.bed
#S288C  5342    5343
#es5363.bed
#S288C  5318    5363
####################################################################

##Load modules
module purge
module load samtools/intel/1.3.1

##Retrieve reads
#for the whole ES27b
samtools view -h -L es7b.bed ERR1309511_rDNA.bam > ERR1309511_es7b.sam
#for variant positions - do sequentially
samtools view -h -L es5322.bed ERR1309511_es7b.sam > 5322.sam
samtools view -h -L es5338.bed 5322.sam > 5338.sam
samtools view -h -L es5343.bed 5338.sam > 5343.sam
samtools view -h -L es5349.bed 5343.sam > 5349.sam
samtools view -h -L es5360.bed 5349.sam > ERR1309511_es7b_var_pos.sam
#Result: ERR1309511_es7b_var_pos.sam contains read mates which span across all five variant positions

#remove intermediate files
rm ERR1309511_es7b.sam
rm 5322.sam
rm 5338.sam
rm 5343.sam
rm 5349.sam
