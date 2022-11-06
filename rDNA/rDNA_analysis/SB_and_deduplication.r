##################    Filter Strand bias (SB) and deduplicate data after "Alignment and variant call"    #################
#Daniel Sultanov                                                                                                         #                                                #
#Original: Jun 2020; modified: Feb 2022                                                                                  #
#                                                                                                                        #
#post-process .vcf file with all isolates has been already processed and pre-filtered by "Alignment and variant call.sh: #
#10 < POS < 9100, GC < 0.6, VF>0.004, HRUN < 4                                                                           #
##########################################################################################################################

##read table created after "Alignment and variant call"
isolates <- read.table("~/path/to/filtered/variant/list/isolates.txt",
                       "\t",
                       header = F,
                       fill = T)

#add column names if they are not present
colnames(isolates) <- c("ACCESSION","CHROM","POS","ID","REF","ALT","QUAL","FILTER","DP","VF","SB","DP4","INDEL","HRUN")

##Filter strand bias
#From Guo et al. BMC Genomics 2012, 13:666, they removed top 10% SB values as a cutoff (or 90% percentile) instead of a fixed value.
#top fixed value in the data is 2147483647 (remove it as well). reported SB values are PHRED-scaled
#note: I removed ~ top 12% *entries* with highest SB (top 12% SB value itself is different)

quantile((subset(isolates, SB!=2147483647))$SB, 0.8855)
isolates<-subset(isolates, SB<144)

##Deduplicate
#remove replicates; one isolate (BAM) containted 12 replicates

isolates<-subset(isolates, !(ACCESSION %in% c("ERR1308597",
"ERR1308642",
"ERR1309125",
"ERR1309273",
"ERR1309040",
"ERR1308817",
"ERR1309482",
"ERR1308858",
"ERR1309524",
"ERR1308638",
"ERR1309161",
"ERR1308663",
"ERR1308819",
"ERR1308802",
"ERR1308803",
"ERR1308804",
"ERR1308778",
"ERR1309311",
"ERR1309243",
"ERR1309322",
"ERR1309093",
"ERR1309006",
"ERR1309053",
"ERR1309392",
"ERR1308699",
"ERR1309402",
"ERR1308628",
"ERR1308626",
"ERR1308728",
"ERR1309022",
"ERR1309346"))
