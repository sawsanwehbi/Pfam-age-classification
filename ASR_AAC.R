# set library path to my home directory where packages are installed
myPaths <- .libPaths()
myPaths <- c(myPaths, "/home/u27/sawsanwehbi/R/x86_64-pc-linux-gnu-library/4.0")
myPaths <- c(myPaths[2], myPaths[1])
.libPaths(myPaths)  # add new path

# Load the ape package
library(ape)
library(seqinr)
library(bioseq)
library(stringr)
library(bio3d)

PFAM_AAC <- read.csv('CSV_files/Pfam_data.csv', header=T)
Pfam_ASRprobabilities <- read.csv('Conserved_AAC/pfam_asr_aac_5seq.csv', header = T)
colnames(Pfam_ASRprobabilities)[3:22] <- gsub('p_',"", colnames(Pfam_ASRprobabilities)[3:22])
#PFAM_IDs <- PFAM_AAC$pfamIDs

# Run to get ancestral frequencies of TM sites in TM pfams
ancestral_tm_sites_df <- read.csv( 'ancestral_tm_sites.csv', header = T)
PFAM_IDs <- ancestral_tm_sites_df$TM_Pfams

## Get ASR probability distribution of AA frequencies at NODE 1
pfam_asr_aac_df <- data.frame(matrix(ncol = 21))

for ( i in 1:length(PFAM_IDs)) {
  print(i) 
  Pfam_Asr <- read.table(paste0('Conserved_ASR/', PFAM_IDs[i] ,'_nogap_alignment.fasta.state'), header = T)
  Root_asr <- Pfam_Asr[which(Pfam_Asr$Node == 'Node1'),]
  Root_asr <- Root_asr[as.numeric(str_split_1(ancestral_tm_sites_df$Ancestral_TMsites[i], ',')),]
  max_prob <- sapply(1:nrow(Root_asr), function(x){max(Root_asr[x,4:23])})
  Root_asr <- Root_asr[which(max_prob >  0.4),]
  len_con_asr <- nrow(Root_asr)
  names(len_con_asr) <- 'Conserved_length'
  pfam_asr_prob_aac <- colSums(Root_asr[,4:23])/len_con_asr 
  pfam_asr_seq <- c(pfam_asr_prob_aac , len_con_asr) 
  pfam_asr_aac_df <- rbind(pfam_asr_aac_df,   pfam_asr_seq ) }

pfam_asr_aac_df <- pfam_asr_aac_df[-1,]
colnames(pfam_asr_aac_df) <- names( pfam_asr_seq)
pfam_asr_aac_df <- cbind(PFAM_IDs, pfam_asr_aac_df)
#write.csv(pfam_asr_aac_df, 'pfam_asr_aac_5seq0.4.csv') 
write.csv(pfam_asr_aac_df, 'pfam_asr_aac_5seq0.4_TM.csv') 

