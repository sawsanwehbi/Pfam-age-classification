## This code classifies Clan ages based on the ages of its Pfam constituents

PFAM_AAC <- read.csv('Pfam_data.csv', header = T)
PFAM_AAC$ancestor[which(PFAM_AAC$avg_transfer > 0.6)]  <- 'unclassifiable'
unique_Clans <- unique(PFAM_AAC$Clans)
Clan_ancestor <- vector()
# 19 CLANS have LACA and LBCA pfams , if we consider them LUCA R=0.37, G is first
for ( i in 1:length(unique_Clans)) {
  if (length(which(grepl('LUCA', PFAM_AAC$ancestor[which(PFAM_AAC$Clans %in% unique_Clans[i])]))) > 1)
  {Clan_ancestor[i] <- 'preLUCA'} 
  else if (length(which(grepl('preLUCA', PFAM_AAC$ancestor[which(PFAM_AAC$Clans %in% unique_Clans[i])]))) > 0)
  { Clan_ancestor[i] <- 'preLUCA'} 
  else if (length(which(grepl('LUCA', PFAM_AAC$ancestor[which(PFAM_AAC$Clans %in% unique_Clans[i])]))) == 1)
  {Clan_ancestor[i] <- 'LUCA'} 
  else if (  'LBCA' %in% PFAM_AAC$ancestor[which(PFAM_AAC$Clans %in% unique_Clans[i])] &
             'LACA' %in% PFAM_AAC$ancestor[which(PFAM_AAC$Clans %in% unique_Clans[i])])
  {Clan_ancestor[i] <- 'LUCA'} 
  else if (  'LBCA' %in% PFAM_AAC$ancestor[which(PFAM_AAC$Clans %in% unique_Clans[i])] |
             'LACA' %in% PFAM_AAC$ancestor[which(PFAM_AAC$Clans %in% unique_Clans[i])]
             & !'unclassifiable' %in% PFAM_AAC$ancestor[which(PFAM_AAC$Clans %in% unique_Clans[i])]) 
  {Clan_ancestor[i] <- 'postLUCA'} 
  else if((length(which(PFAM_AAC$ancestor[which(PFAM_AAC$Clans %in% unique_Clans[i])] %in% c('post-LBCA',Bacterial_supergroups))) > 2) &
          (length(which(PFAM_AAC$ancestor[which(PFAM_AAC$Clans %in% unique_Clans[i])] %in% Archaeal_supergroups)) > 1))
  {Clan_ancestor[i] <- 'unclassifiable' }
  else if((length(which(unique(PFAM_AAC$ancestor[which(PFAM_AAC$Clans %in% unique_Clans[i])]) %in% c('post-LBCA',Bacterial_supergroups))) > 2) |
          (length(which(unique(PFAM_AAC$ancestor[which(PFAM_AAC$Clans %in% unique_Clans[i])]) %in% Archaeal_supergroups)) > 1) 
          & !'unclassifiable' %in% PFAM_AAC$ancestor[which(PFAM_AAC$Clans %in% unique_Clans[i])]) 
  {Clan_ancestor[i] <- 'postLUCA' }
  else if((length(which(unique(PFAM_AAC$ancestor[which(PFAM_AAC$Clans %in% unique_Clans[i])]) %in% c('post-LBCA',Bacterial_supergroups))) == 1) |
          (length(which(unique(PFAM_AAC$ancestor[which(PFAM_AAC$Clans %in% unique_Clans[i])]) %in% Archaeal_supergroups)) == 1) | 
          (length(which(unique(PFAM_AAC$ancestor[which(PFAM_AAC$Clans %in% unique_Clans[i])]) %in% c('post-LBCA',Bacterial_supergroups))) == 2)
          & !'unclassifiable' %in% PFAM_AAC$ancestor[which(PFAM_AAC$Clans %in% unique_Clans[i])]) 
  {Clan_ancestor[i] <- 'modern' }
  else {  Clan_ancestor[i] <- 'unclassifiable'}}

Clan_ancestors_df <- data.frame(unique_Clans,Clan_ancestor)
colnames(Clan_ancestors_df) <- c("Clans", "Clan_ancestor")
write.csv(Clan_ancestors_df, 'Clan_ancestors.csv', row.names = FALSE)
