# set library path to my home directory where packages are installed
myPaths <- .libPaths()
myPaths <- c(myPaths, "/home/u27/sawsanwehbi/R/x86_64-pc-linux-gnu-library/4.0")
myPaths <- c(myPaths[2], myPaths[1])
.libPaths(myPaths)  # add new path

# Load the ape package
library(ape)
library(protr)
library(seqinr)
library(bio3d)


#ConAAsites_df <- data.frame(matrix(ncol= 21))
#for (pfam in 1:length(PFAM_IDs)) {
#PFAM_aln <- seqinr::read.fasta(paste0('IQ_Trees/', PFAM_IDs[pfam] ,'_alignment.fasta'),as.string = F)

# getting conserved sites using conservation score from bio3d package
#PFAM_bio3daln <- bio3d::read.fasta(paste0('IQ_Trees/', PFAM_IDs[pfam] ,'_alignment.fasta'))
#conservation_Score <- conserv(PFAM_bio3daln, 'identity', 'blosum62', normalize.matrix = TRUE)
#con_seq <- toupper(gsub('-','',paste0(PFAM_bio3daln$ali[,which(conservation_Score > 0.75)], collapse = '')))
#con_seq <- gsub('X','', con_seq)
#con_seq <- gsub('J','', con_seq)
#con_seq <- gsub('B','', con_seq)
#con_seq <- gsub('U','', con_seq)
#con_seq <- gsub('O','', con_seq)
#con_seq <- gsub('Z','', con_seq)
#ConAAfreq <- extractAAC(con_seq)
#ConLen <- nchar(con_seq)

#FreqMaxAApersite <- vector()
#MaxAApersite <- vector()
#for ( i in 1:length(PFAM_aln[[1]])){
#SiteAA <- sapply(1:length(PFAM_aln), function(x){PFAM_aln[[x]][i]})
#NbMaxAApersite <- table(SiteAA)[which.max(table(SiteAA))]
#FreqMaxAApersite[i] <- NbMaxAApersite/(length(PFAM_aln) - length(which(SiteAA == '-'))) 
#MaxAApersite[i] <- names(table(SiteAA)[which.max(table(SiteAA))])}

# Remove indels from vectors
#NonIndelMaxAApersite <- MaxAApersite[which( MaxAApersite != '-')]
#NonIndelFreqMaxAApersite <- FreqMaxAApersite[which( MaxAApersite != '-')]

#ConAAsites <- NonIndelMaxAApersite[which(NonIndelFreqMaxAApersite 
  #           > median(NonIndelFreqMaxAApersite ))]

## getting AAC for an alignment of the conserved sites
#ConAlnsites <- which(FreqMaxAApersite > quantile(NonIndelFreqMaxAApersite,0.10) & MaxAApersite != '-' )
#ConAlnAA <- unlist(lapply(1:length(PFAM_aln), function(x){PFAM_aln[[x]][ConAlnsites]}))
#ConAlnAA <- ConAlnAA[which( ConAlnAA != '-')]
#ConAlnAA <- ConAlnAA[which( ConAlnAA != 'u')]
#ConAlnAA <- ConAlnAA[which( ConAlnAA != 'o')]
#ConAlnAA <- ConAlnAA[which( ConAlnAA != 'x')]
#ConAlnAA <- ConAlnAA[which( ConAlnAA != 'j')]
#ConAlnAA <- ConAlnAA[which( ConAlnAA != 'b')]
#ConAlnAA <- ConAlnAA[which( ConAlnAA != 'z')]
#ConAAfreq <- extractAAC(toupper(paste0(ConAlnAA, collapse = '')))
#ConLen <- length(ConAlnAA)

#names(ConLen) <- 'Conserved_concatlength'
#ConAAfreq <- c(ConAAfreq, ConLen)
#ConAAsites_df <- rbind(ConAAsites_df , ConAAfreq) }

#colnames(ConAAsites_df) <- c(names(ConAAfreq))
#ConAAsites_df <- ConAAsites_df[-1,]
#ConAAsites_DF <- cbind(PFAM_IDs, ConAAsites_df)
#sum(ConAAsites_DF[2,2:21])
#write.csv(ConAAsites_DF,'Pfam_ConScore0.75AAC.csv', row.names = F)
#write.csv(ConAAsites_DF,'Pfam_ConAln10quantileAAC.csv', row.names = F)

PFAM_AAC <- read.csv('CSV_files/Pfam_data.csv', header=T)
PFAM_IDs <- PFAM_AAC$pfamIDs
conseq_len <- vector()
i <- 33
for ( i in 1:length(PFAM_IDs)) {
print(i) 
Pfam_tree <- read.tree(paste0('EnforceRoot_trees/', PFAM_IDs[i] ,'_enforceroot.treefile')) 
PFAM_bio3daln <- bio3d::read.fasta(paste0('IQ_Trees/', PFAM_IDs[i] ,'_alignment.fasta'))
#PFAM_bio3daln <- bio3d::read.fasta(paste0('Fasta_files/GeneRaxInvalidPFAMs_fixed/', PFAM_IDs[i] ,'_alignment.fasta'))
# Gap frequency in all the alignment, in bacteria only and in archaea only
AllGapfreq_site <- sapply(1:ncol(PFAM_bio3daln$ali), function(x){
  length(which(PFAM_bio3daln$ali[,x] == '-'))/nrow(PFAM_bio3daln$ali)})
Bac_Gapfreq_site <- sapply(1:ncol(PFAM_bio3daln$ali), function(x){
  length(which(PFAM_bio3daln$ali[which(grepl('Bacteria',PFAM_bio3daln$id)),x] == '-'))/
    length(which(grepl('Bacteria',PFAM_bio3daln$id)))})
Bac_nonGapnb_site <- sapply(1:ncol(PFAM_bio3daln$ali), function(x){
  length(which(PFAM_bio3daln$ali[which(grepl('Bacteria',PFAM_bio3daln$id)),x] != '-'))})
Arc_Gapfreq_site <- sapply(1:ncol(PFAM_bio3daln$ali), function(x){
  length(which(PFAM_bio3daln$ali[which(grepl('Archaea',PFAM_bio3daln$id)),x] == '-'))/
    length(which(grepl('Archaea',PFAM_bio3daln$id)))})
Arc_nonGapnb_site <- sapply(1:ncol(PFAM_bio3daln$ali), function(x){
  length(which(PFAM_bio3daln$ali[which(grepl('Archaea',PFAM_bio3daln$id)),x] != '-'))})

  if (PFAM_AAC$ancestor[i] == 'LUCA' | PFAM_AAC$ancestor[i] == 'preLUCA') {
    Pfam_nogap_aln <- PFAM_bio3daln$ali[,which(AllGapfreq_site < 0.5 &
                      #Bac_Gapfreq_site < 0.8 & Arc_Gapfreq_site < 0.8)]}
                      Bac_nonGapnb_site > 4 & Arc_nonGapnb_site > 4)]}
  else {Pfam_nogap_aln <- PFAM_bio3daln$ali[,which(AllGapfreq_site < 0.5)] }

conseq_len[i] <- ncol(Pfam_nogap_aln)

  allgap_seqs <- sapply(1:nrow(Pfam_nogap_aln),function(x){all(Pfam_nogap_aln[x,] == '-')})
  if (any(allgap_seqs)) {
    Pfam_tree_dropped <- drop.tip(Pfam_tree , rownames(Pfam_nogap_aln)[which(allgap_seqs)])
    #Pfam_nogap_aln <- Pfam_nogap_aln[-which(allgap_seqs),]
    bio3d::write.fasta(ids= rownames(Pfam_nogap_aln)[-which(allgap_seqs)], seqs=Pfam_nogap_aln[-which(allgap_seqs),],
    gap= TRUE, file=paste0('Conserved_ASR/',PFAM_IDs[i],'_nogap_alignment.fasta')) 
    Pfam_tree_dropped$node.label <- 1:length(Pfam_tree_dropped$node.label)
    write.tree(Pfam_tree_dropped,paste0('Conserved_ASR/', PFAM_IDs[i] ,'_ASR.treefile') ) }
 else { Pfam_tree$node.label <- 1:length(Pfam_tree$node.label)
 bio3d::write.fasta(ids=PFAM_bio3daln$id, seqs=Pfam_nogap_aln,
      gap= TRUE, file=paste0('Conserved_ASR/',PFAM_IDs[i],'_nogap_alignment.fasta')) 
   write.tree(Pfam_tree, paste0('Conserved_ASR/', PFAM_IDs[i] ,'_ASR.treefile')) }  }

write.csv(conseq_len , 'Conservedseq_ASR_lenmorethan4.csv')

