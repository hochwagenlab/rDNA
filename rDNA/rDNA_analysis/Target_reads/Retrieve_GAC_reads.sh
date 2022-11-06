####                                  Retrieve GAC reads                                  ####
# Original: May 2021; modified: Feb 2022                                                     #
# Daniel Sultanov                                                                            #
# Prerequisites:                                                                             #
# -mapped rDNA reads from ERR1309433 (ERR1309433_rDNA.bam; see "Alignment and variant call") #
# -annotation files for region of interest (.bed; see below how to create)                   #
##############################################################################################

######################################### GAC in 25S ##########################################
#25S rRNA coordinates: 1226-1283
#From CNE analysis I found 2 GAC varinats in ERR1309433 with VF 28%
#Check if the variants are present in the same set of rDNA copies (in the same short reads)
#Note: this script works fine but there are probably more efficient and scalable way to do this
###############################################################################################

#################### Make annotation .bed files ####################
#format: CHROM  START       END
#        S288C  0-based     non-inclusive
#Note: use rDNA coordinates (1-base inclusive):

#For the whole GAC:
#GAC.bed 
#S288C  4603    4661

#For variant positions:
#GAC_25S_1248.bed
#S288C  4638    4639
#GAC_25S_1237.bed
#S288C  4649    4650
####################################################################

##Load modules
module purge
module load samtools/intel/1.3.1

##Retrieve reads
#for the whole GAC
samtools view -h -L GAC.bed ERR1309433_rDNA.bam > ERR1309433_GAC.sam
#for variant positions - do sequentially
samtools view -h -L GAC_25S_1248.bed ERR1309433_GAC.sam > ERR1309433_GAC_1248.sam
samtools view -h -L GAC_25S_1237.bed ERR1309433_GAC_1248.sam > ERR1309433_GAC_var_pos.sam
#Result: ERR1309433_GAC_var_pos.sam contains read mates which span across the two variant positions

##Remove intermediate files
rm ERR1309433_GAC.sam
rm ERR1309433_GAC_1248.sam
