# Load the ape package
library(ape)
library(protr)
library(seqinr)
library(tidyr)
library(tidyverse)
library(matrixStats)
library(phytools)
library(phangorn)
library(common)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(stringr)
library(gridtext)
library(grid)
library(ggpubr)
library(bio3d)
library(dplyr)
library(weights)
library(philentropy)
library(deming)
library(scales)
library("msa")
library('geoR')
setwd("/PFAMTrees")

      ####---------------------------------------------------------------####
      ###                   Main analysis and figures                    ###
      ####---------------------------------------------------------------####

### Ancestral aa usage analysis ####
Bacterial_supergroups <- c('CPR', 'PVC', 'FCB', 'Terrabacteria', 'Proteobacteria')
Archaeal_supergroups <- c('Asgard', 'TACK', 'DPANN', 'Euryarchaeota')
All_supergroups <- c(Bacterial_supergroups,Archaeal_supergroups)
AA_properties <- read.csv('AminoAcid_properties.csv', header = T)
Moody_pfams <- read.csv('MoodyPfams_probabilities.csv', header = T)
# The following dataframe contains the ancestral amino acid frequencies of conserved regions of Pfams.
#15 Pfams have been reclassified as 'unclassifiable' due to their high transfer rate estimated by generax
Pfam_ConAAC <- read.csv('Pfam_data_ancestralAAC.csv', header = T)

# Run to get the stringent set of pfams confirmed by moody
Moody_pfams <- read.csv('MoodyPfams_probabilities.csv', header = T)
Moody_nonLUCApfams <- Moody_pfams$pfam_id[which(Moody_pfams$greater_than_75 == 0)]
#Pfam_ConAAC <- Pfam_ConAAC[-which(Pfam_ConAAC$ancestor == 'LUCA' & Pfam_ConAAC$pfamIDs %in% Moody_nonLUCApfams |
 #Pfam_ConAAC$ancestor == 'preLUCA' & Pfam_ConAAC$pfamIDs %in% Moody_nonLUCApfams),]

uniqueClans <- unique(Pfam_ConAAC$clans)
Clan_ConAAfreq_Df <- data.frame(uniqueClans)
for (aa in 2:21) {
  Clan_ConAAfreq <- sapply(1:length(uniqueClans), function (i){
    sum(Pfam_ConAAC[,aa][which(Pfam_ConAAC$clans == uniqueClans[i])]*
          Pfam_ConAAC$Conserved_length[which(Pfam_ConAAC$clans == uniqueClans[i])])/
      sum(Pfam_ConAAC$Conserved_length[which(Pfam_ConAAC$clans == uniqueClans[i])])})
  Clan_ConAAfreq_Df <- cbind(Clan_ConAAfreq_Df , Clan_ConAAfreq) }
colnames(Clan_ConAAfreq_Df)[2:21] <- colnames(Pfam_ConAAC)[2:21]
Clan_ConAAfreq_Df$Clan_Conlength <- sapply(1:length(uniqueClans), function (i){
  max(Pfam_ConAAC$Conserved_length[which(Pfam_ConAAC$clans == uniqueClans[i])])})
colnames(Clan_ConAAfreq_Df)[1] <- 'Clans'

Clan_ancestor <- vector() 
for ( i in 1:length(uniqueClans)) {
  if (length(which(grepl('LUCA', Pfam_ConAAC$ancestor[which(Pfam_ConAAC$clans %in% uniqueClans[i])]))) > 1)
  {Clan_ancestor[i] <- 'preLUCA'} 
  else if (length(which(grepl('preLUCA', Pfam_ConAAC$ancestor[which(Pfam_ConAAC$clans %in% uniqueClans[i])]))) > 0)
  { Clan_ancestor[i] <- 'preLUCA'} 
  else if (length(which(grepl('LUCA', Pfam_ConAAC$ancestor[which(Pfam_ConAAC$clans %in% uniqueClans[i])]))) == 1)
  {Clan_ancestor[i] <- 'LUCA'} 
  else if (  'LBCA' %in% Pfam_ConAAC$ancestor[which(Pfam_ConAAC$clans %in% uniqueClans[i])] &
             'LACA' %in% Pfam_ConAAC$ancestor[which(Pfam_ConAAC$clans %in% uniqueClans[i])])
  {Clan_ancestor[i] <- 'LUCA'} 
  else if (  'LBCA' %in% Pfam_ConAAC$ancestor[which(Pfam_ConAAC$clans %in% uniqueClans[i])] |
             'LACA' %in% Pfam_ConAAC$ancestor[which(Pfam_ConAAC$clans %in% uniqueClans[i])]
             & !'unclassifiable' %in% Pfam_ConAAC$ancestor[which(Pfam_ConAAC$Clans %in% uniqueClans[i])]) 
  {Clan_ancestor[i] <- 'postLUCA'} 
  else if((length(which(Pfam_ConAAC$ancestor[which(Pfam_ConAAC$clans %in% uniqueClans[i])] %in% c('post-LBCA',Bacterial_supergroups))) > 2) &
          (length(which(Pfam_ConAAC$ancestor[which(Pfam_ConAAC$clans %in% uniqueClans[i])] %in% Archaeal_supergroups)) > 1))
  {Clan_ancestor[i] <- 'unclassifiable' }
  else if((length(which(unique(Pfam_ConAAC$ancestor[which(Pfam_ConAAC$clans %in% uniqueClans[i])]) %in% c('post-LBCA',Bacterial_supergroups))) > 2) |
          (length(which(unique(Pfam_ConAAC$ancestor[which(Pfam_ConAAC$clans %in% uniqueClans[i])]) %in% Archaeal_supergroups)) > 1) 
          & !'unclassifiable' %in% Pfam_ConAAC$ancestor[which(Pfam_ConAAC$clans %in% uniqueClans[i])]) 
  {Clan_ancestor[i] <- 'postLUCA' }
  else if((length(which(unique(Pfam_ConAAC$ancestor[which(Pfam_ConAAC$clans %in% uniqueClans[i])]) %in% c('post-LBCA',Bacterial_supergroups))) == 1) |
          (length(which(unique(Pfam_ConAAC$ancestor[which(Pfam_ConAAC$clans %in% uniqueClans[i])]) %in% Archaeal_supergroups)) == 1) | 
          (length(which(unique(Pfam_ConAAC$ancestor[which(Pfam_ConAAC$clans %in% uniqueClans[i])]) %in% c('post-LBCA',Bacterial_supergroups))) == 2)
          & !'unclassifiable' %in% Pfam_ConAAC$ancestor[which(Pfam_ConAAC$clans %in% uniqueClans[i])]) 
  {Clan_ancestor[i] <- 'modern' }
  else {  Clan_ancestor[i] <- 'unclassifiable'} }
Clan_ConAAfreq_Df$Clan_ancestor <- Clan_ancestor
#write.csv(Clan_ConAAfreq_Df, 'Clan_data_ancestralAAC.csv', row.names = F)

allLUCA_Clans_ConAAC <- colWeightedMeans(as.matrix(Clan_ConAAfreq_Df[which(grepl('LUCA',Clan_ConAAfreq_Df$Clan_ancestor) & 
                        Clan_ConAAfreq_Df$Clan_ancestor != 'postLUCA'),2:21]),
                        Clan_ConAAfreq_Df$Clan_Conlength[which(grepl('LUCA',Clan_ConAAfreq_Df$Clan_ancestor) &
                        Clan_ConAAfreq_Df$Clan_ancestor != 'postLUCA')] )
LUCA_Clans_ConAAC <- colWeightedMeans(as.matrix(Clan_ConAAfreq_Df[which(Clan_ConAAfreq_Df$Clan_ancestor == 'LUCA'),2:21]),
                      Clan_ConAAfreq_Df$Clan_Conlength[which(Clan_ConAAfreq_Df$Clan_ancestor == 'LUCA')] )
preLUCA_Clans_ConAAC <- colWeightedMeans(as.matrix(Clan_ConAAfreq_Df[which(Clan_ConAAfreq_Df$Clan_ancestor == 'preLUCA'),2:21]),
                       Clan_ConAAfreq_Df$Clan_Conlength[which(Clan_ConAAfreq_Df$Clan_ancestor == 'preLUCA')] )
postLUCA_Clans_ConAAC <- colWeightedMeans(as.matrix(Clan_ConAAfreq_Df[which(Clan_ConAAfreq_Df$Clan_ancestor == 'postLUCA'),2:21]),
                         Clan_ConAAfreq_Df$Clan_Conlength[which(Clan_ConAAfreq_Df$Clan_ancestor == 'postLUCA')] )
modern_Clans_ConAAC <- colWeightedMeans(as.matrix(Clan_ConAAfreq_Df[which(Clan_ConAAfreq_Df$Clan_ancestor == 'modern'),2:21]),
                       Clan_ConAAfreq_Df$Clan_Conlength[which(Clan_ConAAfreq_Df$Clan_ancestor == 'modern')] )

####### SAM-dependent enzymes #######
# COG2890, gene=PRMC ;PF13649 LUCA, PF17827 LBCA
# COG0621, gene=MIAB ;PF01938 LUCA, PF04055 PRELUCA, PF00919  LBCA
# C0G0144, gene=rsmF ;PF01189 LUCA, PF17125 LUCA 
#Pfam_ConAAC[which(Pfam_ConAAC$pfamIDs =='PF17284'),] 
#weiss_LUCA <- read.csv('WeissPfams_LUCA.csv', header = T, comment.char = '#')
#Moody_LUCA <- read.csv('STable_1.csv', header = T, comment.char = '#')
#WeissSAM_COGs <- weiss_LUCA$COG.ID[which(grepl('SAM',weiss_LUCA$cofactora))]

#### standard errors for ratios #####
library(diagis)
allLUCAweightedse <- vector()
LUCAweightedse <- vector()
preLUCAweightedse <- vector()
postLUCAweightedse <- vector()
colnames(Clan_ConAAfreq_Df)
for (colnb in 2:21) {
  allLUCAweightedse[colnb] <- weighted_se(Clan_ConAAfreq_Df[,colnb][which(grepl('LUCA',Clan_ConAAfreq_Df$Clan_ancestor))],
                                          Clan_ConAAfreq_Df$Clan_Conlength[which(grepl('LUCA',Clan_ConAAfreq_Df$Clan_ancestor))])}
for (colnb in 2:21) {
  LUCAweightedse[colnb] <- weighted_se(Clan_ConAAfreq_Df[,colnb][which(Clan_ConAAfreq_Df$Clan_ancestor == 'LUCA')],
                                       Clan_ConAAfreq_Df$Clan_Conlength[which(Clan_ConAAfreq_Df$Clan_ancestor == 'LUCA')])}
for (colnb in 2:21){
  postLUCAweightedse[colnb] <- weighted_se(Clan_ConAAfreq_Df[,colnb][which(Clan_ConAAfreq_Df$Clan_ancestor == 'postLUCA')],
                                           Clan_ConAAfreq_Df$Clan_Conlength[which(Clan_ConAAfreq_Df$Clan_ancestor == 'postLUCA')])}
for (colnb in 2:21){
  preLUCAweightedse[colnb] <- weighted_se(Clan_ConAAfreq_Df[,colnb][which(Clan_ConAAfreq_Df$Clan_ancestor == 'preLUCA')],
                                          Clan_ConAAfreq_Df$Clan_Conlength[which(Clan_ConAAfreq_Df$Clan_ancestor == 'preLUCA')])}

LUCAclanratio_var <-  (LUCAweightedse[-1]^2)/((postLUCA_Clans_ConAAC)^2) +
  (postLUCAweightedse[-1]^2)*((LUCA_Clans_ConAAC)^2)/((postLUCA_Clans_ConAAC)^4)
LUCAclanratio_se <- sqrt(LUCAclanratio_var)
LUCAclanratio_se <- LUCAclanratio_se[match(AA_properties$Letter,names(LUCAclanratio_se))]
preLUCAclanratio_var <-  (preLUCAweightedse[-1]^2)/((postLUCA_Clans_ConAAC)^2) +
  (postLUCAweightedse[-1]^2)*((preLUCA_Clans_ConAAC)^2)/((postLUCA_Clans_ConAAC)^4)
preLUCAclanratio_se <- sqrt(preLUCAclanratio_var)
preLUCAclanratio_se <- preLUCAclanratio_se[match(AA_properties$Letter,names(preLUCAclanratio_se))]

##### Clan conserved AAC analysis ######
ancient_clanConusage <- LUCA_Clans_ConAAC/postLUCA_Clans_ConAAC
ancient_clanConusage <- ancient_clanConusage[match(AA_properties$Letter,names(ancient_clanConusage))]
clan_wt <- 1 / (AA_properties$Filtered_sd_2000^2)
Conclan_model <- lm(AA_properties$Filtered_avg_2000 ~ ancient_clanConusage)
Conclan_wls_model <- lm(AA_properties$Filtered_avg_2000 ~ ancient_clanConusage,
                        weights = clan_wt)
Conclan_wls_confinterval <- broom::augment(Conclan_wls_model , interval="confidence")
summary(Conclan_wls_model)
sort(ancient_clanConusage, decreasing = T)
summary(lm(AA_properties$TrifonovRecalculated ~ ancient_clanConusage,
           weights = clan_wt))
summary(lm(ancient_clanConusage ~ AA_properties$MolecularWeightsDa,
           weights = 1/(LUCAclanratio_se^2)   ))
summary(lm(ancient_clanConusage ~ AA_properties$MolecularWeightsDa + 
             AA_properties$Filtered_avg_2000, weights = 1/(LUCAclanratio_se^2)))

## usage vs Trifonov 2004
summary(lm( ancient_clanConusage ~ AA_properties$Average_rank_Trifonov2004,
           weights = 1/(LUCAclanratio_se^2) ) )
summary(lm(ancient_clanConusage ~ AA_properties$MolecularWeightsDa + 
          AA_properties$Average_rank_Trifonov2004,     weights = 1/(LUCAclanratio_se^2)))


ancientpreLUCA_clanConusage <- preLUCA_Clans_ConAAC/postLUCA_Clans_ConAAC
ancientpreLUCA_clanConusage <- ancientpreLUCA_clanConusage[match(AA_properties$Letter,names(ancientpreLUCA_clanConusage))]
Conclan_preluca_wls_model <- lm(AA_properties$Filtered_avg_2000 ~ ancientpreLUCA_clanConusage,
                                weights = clan_wt)
summary(Conclan_preluca_wls_model)
summary(lm(ancientpreLUCA_clanConusage~ AA_properties$MolecularWeightsDa))
summary(lm(ancientpreLUCA_clanConusage~ AA_properties$MolecularWeightsDa,
           weights = 1/(preLUCAclanratio_se^2) ) )
summary(lm(ancientpreLUCA_clanConusage ~ ancient_clanConusage))

ancientpostLUCA_clanConusage <- postLUCA_Clans_ConAAC/modern_Clans_ConAAC
ancientpostLUCA_clanConusage <- ancientpostLUCA_clanConusage[match(AA_properties$Letter,names(ancientpostLUCA_clanConusage))]
summary(lm(ancientpostLUCA_clanConusage ~ AA_properties$MolecularWeightsDa))

#### Plot conserved ancient aa usage in luca and preluca####
rm(ordervsusage_df )
ordervsusage_df <- cbind(AA_properties$Filtered_avg_2000 ,AA_properties$Filtered_sd_2000,
                         ancient_clanConusage,  
                         ancientpreLUCA_clanConusage , AA_properties$Moosmann_category )
ordervsusage_df <- as.data.frame(ordervsusage_df)
rownames(ordervsusage_df) <- AA_properties$Letter
colnames(ordervsusage_df) <- c('Trifonov_order', 'Trifonov_se',
                               'ancient_clanusage', 'preluca_clanusage','Moosmann_category')
conLUCA_plot <- ggplot(ordervsusage_df, aes(y = as.numeric(Trifonov_order), 
  x = as.numeric(ancient_clanusage), label = rownames(ordervsusage_df))) + 
  ylab('Trifonov (2000) order') +
  xlab('LUCA clan usage') + 
  geom_text(aes(colour = factor(Moosmann_category)), size = 16) + 
  labs(color='Moosmann category') + 
  theme(axis.text=element_text(size=24),axis.title=element_text(size=32,face="bold"),
  legend.text=element_text(size=27), legend.title = element_text(size=16,face="bold")) + 
  scale_fill_brewer(palette="Dark2") + ylim(0,27) + xlim(0.63,1.14) +
  theme(legend.position = 'none') +
  guides(color = guide_legend(override.aes = list(size = 16))) + labs(color=NULL) + 
  geom_line(aes(y = predict(Conclan_wls_model)), linewidth = 1,color = 'black') +
  annotate(geom="text", y=26, x=0.9, label=paste0("Weighted ",paste0('R',supsc('2')), "= 0.37"), color="black", size = 14) +
  annotate(geom="text", y=23, x=0.9, label="p = 0.003", color="black", size = 14) +
  annotate(geom="text", y=26, x=0.63, label="b)", color="black", size = 16) +
  geom_ribbon(aes(ymin=Conclan_wls_confinterval$.lower, ymax=Conclan_wls_confinterval$.upper), colour=NA, alpha=0.3)

#### Molecular weights vs ancient usage plots #####
AA_properties <- cbind(AA_properties , ancient_clanConusage , ancientpreLUCA_clanConusage )
mweights_model <- lm( ancient_clanConusage ~ AA_properties$MolecularWeightsDa,
                      weights = 1/(LUCAclanratio_se^2))
mweights_model_confinterval <- broom::augment(mweights_model, interval="confidence")
prelucamweights_model <- lm( ancientpreLUCA_clanConusage  ~ AA_properties$MolecularWeightsDa,
                             weights = 1/(preLUCAclanratio_se^2))
prelucamweights_model_confinterval <- broom::augment(prelucamweights_model, interval="confidence")
summary(mweights_model)

mweights_plot <- ggplot(AA_properties, aes(y = as.numeric(ancient_clanConusage), 
  x = as.numeric(MolecularWeightsDa), label = Letter)) + 
  xlab('Molecular Weight (Da)') + ylab('LUCA clan usage')  + 
  theme(legend.position="none") + geom_text(size = 16, color='darkblue') + 
  geom_errorbar(aes(ymin=as.numeric(ancient_clanConusage)-LUCAclanratio_se, ymax=as.numeric(ancient_clanConusage)+LUCAclanratio_se)) +
  theme(axis.text=element_text(size=24),axis.title=element_text(size=32,face="bold"),
  legend.text=element_text(size=14), legend.title = element_text(size=16,face="bold")) + 
  geom_line(aes(y = predict(mweights_model )), linewidth = 1,color = 'black') +
  annotate(geom="text", x=140, y=0.71, label=paste0(paste0('Weighted R',supsc('2')), "= 0.48"), color= "black", size = 14) +
  annotate(geom="text", x=140, y=0.64, label="p = 0.0005", color= "black", size = 14) +
  annotate(geom="text", x=70, y=1.2, label="a)", color="black", size = 16) +
  geom_ribbon(aes(ymin=mweights_model_confinterval$.lower, ymax=mweights_model_confinterval$.upper), colour=NA, alpha=0.3)


prelucamweights_plot <- ggplot(AA_properties, aes(y = as.numeric(ancientpreLUCA_clanConusage ), 
  x = as.numeric(MolecularWeightsDa), label = Letter)) + 
  xlab('Molecular Weight (Da)') + ylab('pre-LUCA clan usage')  + 
  theme(legend.position="none") + geom_text(size = 16, color='darkblue') + 
  geom_errorbar(aes(ymin=as.numeric(ancientpreLUCA_clanConusage )-preLUCAclanratio_se, ymax=as.numeric(ancientpreLUCA_clanConusage )+preLUCAclanratio_se)) +
  theme(axis.text=element_text(size=24),axis.title=element_text(size=32,face="bold"),
  legend.text=element_text(size=14), legend.title = element_text(size=16,face="bold")) + 
  geom_line(aes(y = predict(mweights_model )), linewidth = 1,color = 'black') +
  annotate(geom="text", x=140, y=0.71, label=paste0(paste0('Weighted R',supsc('2')), "= 0.33"), color= "black", size = 14) +
  annotate(geom="text", x=140, y=0.64, label="p = 0.007", color= "black", size = 14) +
  annotate(geom="text", x=70, y=1.2, label="c)", color="black", size = 16) +
  geom_ribbon(aes(ymin=mweights_model_confinterval$.lower, ymax=mweights_model_confinterval$.upper), colour=NA, alpha=0.3)

### PreLUCA vs LUCA MODEL 2 regression plot ####
error_ratio <- preLUCAclanratio_se/LUCAclanratio_se
demingfit <- deming(ancientpreLUCA_clanConusage ~ ancient_clanConusage, weights = 1/error_ratio^2,
                    xstd = LUCAclanratio_se, ystd = preLUCAclanratio_se)
# Use a test model and replace the coefficients with those from the york fit
# to get the confidence interval using the predict function
test_model <- lm( ancientpreLUCA_clanConusage ~ ancient_clanConusage )
test_model$coefficients <- demingfit$coefficients
CI_data <- predict(test_model,  interval = "confidence", level = 0.95)
summary(lm( ancientpreLUCA_clanConusage ~ ancient_clanConusage))
prelucavsluca_df <- data.frame( ancientpreLUCA_clanConusage ,ancient_clanConusage)
prelucavsluca_df$Biochem_properties <- c(NA,NA,'Charged','Charged','Aromatic',NA,'Aromatic',NA,'Charged',
                                         NA,NA,NA,NA,NA,'Charged',NA,NA,NA,'Aromatic','Aromatic' )
prelucavsluca_df$Moosmann_category <- AA_properties$Moosmann_category
# WYVFEH are statistically different between luca and preluca
pvalues <- sapply(2:21, function(x) { 
  ttest <- t.test(Clan_ConAAfreq_Df[,x][which(Clan_ConAAfreq_Df$Clan_ancestor == 'LUCA')],
                  Clan_ConAAfreq_Df[,x][which(Clan_ConAAfreq_Df$Clan_ancestor == 'preLUCA')])
  return(ttest$p.value)})
names(pvalues) <- colnames(Pfam_ConAAC)[2:21]

t.test(Clan_ConAAfreq_Df$H[which(Clan_ConAAfreq_Df$Clan_ancestor == 'LUCA')],
       Clan_ConAAfreq_Df$H[which(Clan_ConAAfreq_Df$Clan_ancestor == 'preLUCA')])

preluca_luca_plot <- ggplot(prelucavsluca_df, aes(x = as.numeric(ancient_clanConusage), y = as.numeric( ancientpreLUCA_clanConusage), 
  label = rownames(prelucavsluca_df))) + 
  geom_vline(xintercept=1) + geom_hline(yintercept=1) + xlim(0.63,1.18) + ylim(0.63,1.18) +
  xlab('LUCA clan usage') + theme(legend.key=element_rect(fill="white")) +
  ylab('Pre-LUCA clan usage') + 
  theme(legend.position = 'none') +
  geom_text(aes(colour = factor(Biochem_properties)), size = 16) +    
  labs(color='Biochemical Properties') +  scale_fill_brewer(palette="Dark2") +  
  theme(axis.text=element_text(size=24),axis.title=element_text(size=32,face="bold"),
        legend.text=element_text(size=14), legend.title = element_text(size=16,face="bold")) + 
  geom_abline(aes(slope = demingfit$coefficients[2], intercept =demingfit$coefficients[1]), linewidth = 1,color = 'blue') +
  annotate(geom="text", y=0.7, x=0.87, label=paste0("",paste0('R',supsc('2')), "= 0.51"), color="black", size =14) +
  annotate(geom="text", y=0.67, x=0.87, label="p = 0.0003", color="black", size = 14) +
  #geom_ribbon(aes(ymax=CI_data[,3], ymin=CI_data[,2]), colour=NA, alpha=0.3) +
  annotate(geom="text", x=0.63, y=1.18, label="d)", color="black", size = 16) +
  geom_abline(aes(slope=1, intercept = 0),linewidth = 1,color = 'darkred') + 
  annotate(geom="text", y=ancientpreLUCA_clanConusage[19]+0.012, x=ancient_clanConusage[19], label="**", color="black", size = 16) +
  annotate(geom="text", y=ancientpreLUCA_clanConusage[18]+0.012, x=ancient_clanConusage[18], label="**", color="black", size = 16) +
  annotate(geom="text", y=ancientpreLUCA_clanConusage[5]+0.012, x=ancient_clanConusage[5], label="*", color="black", size = 16) +
  annotate(geom="text", y=ancientpreLUCA_clanConusage[4]+0.012, x=ancient_clanConusage[4], label="*", color="black", size = 16) +
  annotate(geom="text", y=ancientpreLUCA_clanConusage[20]+0.012, x=ancient_clanConusage[20], label="*", color="black", size = 16) +
  annotate(geom="text", y=ancientpreLUCA_clanConusage[7]+0.012, x=ancient_clanConusage[7], label="*", color="black", size = 16) 

# Plot 4-figure panel
grid.arrange(mweights_plot,conLUCA_plot, prelucamweights_plot, preluca_luca_plot)


      ####---------------------------------------------------------------####
      ###             Supplementary analysis and figures                  ###
      ####---------------------------------------------------------------####

### Supplementary tables 1 and 2 ####
usage_df <- signif(data.frame( ancient_clanConusage,  LUCAclanratio_se,
                               ancientpreLUCA_clanConusage , preLUCAclanratio_se),digits = 3)
usage_df <- cbind(names(ancient_clanConusage), usage_df )
colnames(usage_df) <- c('Amino acid', 'LUCA clan usage',	'LUCA clan usage standard error',	
                        'Pre-LUCA clan usage','Pre-LUCA clan usage standard error')
rownames(usage_df) <- NULL
write.csv(usage_df,'usage-df.csv')
plot(tableGrob(usage_df))

LUCA_pathways <- read.csv('LUCA_Moody_pathways.csv', header = T)
colnames(LUCA_pathways )[1:4] <- c('Pathway', 'KEGG ID', 'Enzyme', 'Associated Pfams present in LUCA')
plot(tableGrob(LUCA_pathways[,1:4], theme = ttheme_default(base_size = 13)))

Pfam_ConAAC[which(Pfam_ConAAC$pfamIDs == 'PF06988'),]
#### Transmembrane vs nontransmembrane ####
# replace ancestral amino acid frequencies of transmembrane pfams with
# the ancestral amino acid frequencies of their membrane embedded regions only
Ancestral_TMsites_AAC <- read.csv('pfam_asr_aac_5seq0.4_DEEPTM.csv', header = T)
DeepTmhmm <- readLines('DEEPTmhmm_consensusPfam.out')
tmhmm_out <- DeepTmhmm[which(grepl('Number of predicted TMRs', DeepTmhmm))]
Tmhmm_score <- data.frame(word(tmhmm_out, 2),as.numeric(word(tmhmm_out, -1)))
colnames(Tmhmm_score) <- c('pfamIDs', 'Tmhmm')
Tmhmm_score <- inner_join(Tmhmm_score , Pfam_ConAAC, by = 'pfamIDs')
Tmhmm_score[which(Tmhmm_score$pfamIDs %in% Ancestral_TMsites_AAC$PFAM_IDs),3:23] <- Ancestral_TMsites_AAC[,3:23]

## TM clans usage
Clans_Tmhmm <- vector()
uniqueClans <- unique(Tmhmm_score$clans)
for ( x in 1:length(uniqueClans)) {
  if (all(Tmhmm_score$Tmhmm[which(Tmhmm_score$clans == uniqueClans[x])] == 0)) {
    Clans_Tmhmm[x] <- 0 } 
  else if (all(Tmhmm_score$Tmhmm[which(Tmhmm_score$clans == uniqueClans[x])] > 0)) {
    Clans_Tmhmm[x] <- 1} 
  else { Clans_Tmhmm[x] <- NA} }
Clans_Tmhmm <- data.frame(uniqueClans,Clans_Tmhmm )
Transmembrane_clans <- uniqueClans[which(Clans_Tmhmm$Clans_Tmhmm > 0)] # 798 clans
Nonmembrane_clans <- uniqueClans[which(Clans_Tmhmm$Clans_Tmhmm == 0)] # 3332 clans

TM_uniqueClans <- unique(Tmhmm_score$clans)
TMClan_ConAAfreq_Df <- data.frame(TM_uniqueClans)
for (aa in 3:22) {
  TMClan_ConAAfreq <- sapply(1:length(TM_uniqueClans), function (i){
    sum(Tmhmm_score[,aa][which(Tmhmm_score$clans == TM_uniqueClans[i])]*
          Tmhmm_score$Conserved_length[which(Tmhmm_score$clans == TM_uniqueClans[i])])/
      sum(Tmhmm_score$Conserved_length[which(Tmhmm_score$clans == TM_uniqueClans[i])])})
  TMClan_ConAAfreq_Df <- cbind(TMClan_ConAAfreq_Df , TMClan_ConAAfreq) }
colnames(TMClan_ConAAfreq_Df)[2:21] <- colnames(Tmhmm_score)[3:22]
TMClan_ConAAfreq_Df$Clan_Conlength <- sapply(1:length(TM_uniqueClans), function (i){
  max(Tmhmm_score$Conserved_length[which(Tmhmm_score$clans == TM_uniqueClans[i])])})
colnames(TMClan_ConAAfreq_Df)[1] <- 'Clans'
TMClan_ConAAfreq_Df <-  merge(TMClan_ConAAfreq_Df , Clan_ConAAfreq_Df[,c(1,23)], by = 'Clans')
colnames(Clans_Tmhmm)[1] <- 'Clans'
TMClan_ConAAfreq_Df <-  merge(TMClan_ConAAfreq_Df , Clans_Tmhmm, by = 'Clans')
TMClan_ConAAfreq_Df <- TMClan_ConAAfreq_Df[-which(is.na(TMClan_ConAAfreq_Df$A)),]


Trans_LUCA_clans_AAC <- colWeightedMeans(as.matrix(TMClan_ConAAfreq_Df[which(TMClan_ConAAfreq_Df$Clan_ancestor == 'LUCA' & TMClan_ConAAfreq_Df$Clans_Tmhmm > 0.5 ),2:21]),
                                         TMClan_ConAAfreq_Df$Clan_Conlength[which(TMClan_ConAAfreq_Df$Clan_ancestor == 'LUCA' & TMClan_ConAAfreq_Df$Clans_Tmhmm  > 0.5 )] )
NonTrans_LUCA_clans_AAC <- colWeightedMeans(as.matrix(TMClan_ConAAfreq_Df[which(TMClan_ConAAfreq_Df$Clan_ancestor == 'LUCA' & TMClan_ConAAfreq_Df$Clans_Tmhmm < 0.5 ),2:21]),
                                            TMClan_ConAAfreq_Df$Clan_Conlength[which(TMClan_ConAAfreq_Df$Clan_ancestor == 'LUCA' & TMClan_ConAAfreq_Df$Clans_Tmhmm  < 0.5 )] )
Trans_postLUCA_clans_AAC <- colWeightedMeans(as.matrix(TMClan_ConAAfreq_Df[which(TMClan_ConAAfreq_Df$Clan_ancestor == 'postLUCA' & TMClan_ConAAfreq_Df$Clans_Tmhmm > 0.5 ),2:21]),
                                             TMClan_ConAAfreq_Df$Clan_Conlength[which(TMClan_ConAAfreq_Df$Clan_ancestor == 'postLUCA' & TMClan_ConAAfreq_Df$Clans_Tmhmm  > 0.5 )] )
NonTrans_postLUCA_clans_AAC <- colWeightedMeans(as.matrix(TMClan_ConAAfreq_Df[which(TMClan_ConAAfreq_Df$Clan_ancestor == 'postLUCA' & TMClan_ConAAfreq_Df$Clans_Tmhmm < 0.5 ),2:21]),
                                                TMClan_ConAAfreq_Df$Clan_Conlength[which(TMClan_ConAAfreq_Df$Clan_ancestor == 'postLUCA' & TMClan_ConAAfreq_Df$Clans_Tmhmm  < 0.5 )] )

#### Error bars on x and y axis
TM_LUCAweightedse <- vector()
TM_postLUCAweightedse <- vector()
colnames(TMClan_ConAAfreq_Df)
for (colnb in 2:21) {
  TM_LUCAweightedse[colnb] <- weighted_se(TMClan_ConAAfreq_Df[,colnb][which(TMClan_ConAAfreq_Df$Clan_ancestor == 'LUCA' & TMClan_ConAAfreq_Df$Clans_Tmhmm  > 0.5 )],
                                          TMClan_ConAAfreq_Df$Clan_Conlength[which(TMClan_ConAAfreq_Df$Clan_ancestor == 'LUCA' & TMClan_ConAAfreq_Df$Clans_Tmhmm  > 0.5 )])}
for (colnb in 2:21){
  TM_postLUCAweightedse[colnb] <- weighted_se(TMClan_ConAAfreq_Df[,colnb][which(TMClan_ConAAfreq_Df$Clan_ancestor == 'postLUCA' & TMClan_ConAAfreq_Df$Clans_Tmhmm  > 0.5 )],
                                              TMClan_ConAAfreq_Df$Clan_Conlength[which(TMClan_ConAAfreq_Df$Clan_ancestor == 'postLUCA' & TMClan_ConAAfreq_Df$Clans_Tmhmm  > 0.5 )])}
TM_LUCAclanratio_var <-  (TM_LUCAweightedse[-1]^2)/((Trans_postLUCA_clans_AAC)^2) +
  (TM_postLUCAweightedse[-1]^2)*((Trans_LUCA_clans_AAC)^2)/((Trans_postLUCA_clans_AAC)^4)
TM_LUCAclanratio_se <- sqrt(TM_LUCAclanratio_var)

nonTM_LUCAweightedse <- vector()
nonTM_postLUCAweightedse <- vector()
colnames(TMClan_ConAAfreq_Df)
for (colnb in 2:21) {
  nonTM_LUCAweightedse[colnb] <- weighted_se(TMClan_ConAAfreq_Df[,colnb][which(TMClan_ConAAfreq_Df$Clan_ancestor == 'LUCA' & TMClan_ConAAfreq_Df$Clans_Tmhmm  < 0.5 )],
                                             TMClan_ConAAfreq_Df$Clan_Conlength[which(TMClan_ConAAfreq_Df$Clan_ancestor == 'LUCA' & TMClan_ConAAfreq_Df$Clans_Tmhmm  < 0.5 )])}
for (colnb in 2:21){
  nonTM_postLUCAweightedse[colnb] <- weighted_se(TMClan_ConAAfreq_Df[,colnb][which(TMClan_ConAAfreq_Df$Clan_ancestor == 'postLUCA' & TMClan_ConAAfreq_Df$Clans_Tmhmm  < 0.5 )],
                                                 TMClan_ConAAfreq_Df$Clan_Conlength[which(TMClan_ConAAfreq_Df$Clan_ancestor == 'postLUCA' & TMClan_ConAAfreq_Df$Clans_Tmhmm  < 0.5 )])}
nonTM_LUCAclanratio_var <-  (nonTM_LUCAweightedse[-1]^2)/((NonTrans_postLUCA_clans_AAC)^2) +
  (nonTM_postLUCAweightedse[-1]^2)*((NonTrans_LUCA_clans_AAC)^2)/((NonTrans_postLUCA_clans_AAC)^4)
nonTM_LUCAclanratio_se <- sqrt(nonTM_LUCAclanratio_var)

table(TMClan_ConAAfreq_Df$Clan_ancestor[which(TMClan_ConAAfreq_Df$Clans_Tmhmm < 0.5)])
#63/218 ,  354/998
#45/218 ,  255/998

Trans_LUCAclan_usage <- Trans_LUCA_clans_AAC/Trans_postLUCA_clans_AAC
Nontrans_LUCAclan_usage <- NonTrans_LUCA_clans_AAC/NonTrans_postLUCA_clans_AAC
Trans_nonTransusage_df <- data.frame(Trans_LUCAclan_usage, Nontrans_LUCAclan_usage)
Trans_nonTrans_clan_model <-  lm(Trans_LUCAclan_usage ~ Nontrans_LUCAclan_usage , weights = 1/TM_LUCAclanratio_se )
Trans_nonTrans_clan_model_confinterval <- broom::augment(Trans_nonTrans_clan_model , interval="confidence")
summary(Trans_nonTrans_clan_model)

##Supplementary figure 3
ggplot(Trans_nonTransusage_df, aes(y = as.numeric(Trans_LUCAclan_usage ), 
                                   x = as.numeric( Nontrans_LUCAclan_usage), label = names(Trans_LUCAclan_usage ))) + 
  geom_text(label = names(Trans_LUCAclan_usage ), size = 8) + 
  geom_errorbar(aes(ymin=Trans_LUCAclan_usage-TM_LUCAclanratio_se, ymax=Trans_LUCAclan_usage+TM_LUCAclanratio_se)) +
  geom_errorbar(aes(xmin=Nontrans_LUCAclan_usage-nonTM_LUCAclanratio_se, xmax=Nontrans_LUCAclan_usage+nonTM_LUCAclanratio_se)) +
  ylab('Transmembrane ancestral LUCA clan usage') +
  xlab('Non-transmembrane ancestral LUCA clan usage') + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=24,face="bold"),
        legend.text=element_text(size=14), legend.title = element_text(size=16,face="bold")) + 
  scale_fill_brewer(palette="Dark2") + ylim(0.5,1.9) + xlim(0.5,1.23) +
  #geom_line(aes(y = predict(Trans_nonTrans_clan_model)), linewidth = 1,color = 'black') +
  annotate(geom="text", y=1.2, x=0.7, label=paste0("Weighted ",paste0('R',supsc('2')), "= 0.25"), color="darkblue", size = 8) +
  annotate(geom="text", y=1.14, x=0.7, label="p = 0.02", color="darkblue", size = 8) +
  geom_abline(aes(slope=1, intercept = 0),linewidth = 1,color = 'darkred')  
#geom_ribbon(aes(ymin=Trans_nonTrans_clan_model_confinterval[,2], ymax=Trans_nonTrans_clan_model_confinterval[,3]), colour=NA, alpha=0.3) 

tmhmm_alpha <- DeepTmhmm[which(grepl('TMhelix', DeepTmhmm))]
tmhmm_alpha_pfams <- unique(gsub('\t.*', "", tmhmm_alpha))
tmhmm_beta <- DeepTmhmm[which(grepl('Beta sheet', DeepTmhmm))]
tmhmm_beta_pfams <- unique(gsub('\t.*', "", tmhmm_beta))
table(Pfam_ConAAC$ancestor[which(Pfam_ConAAC$pfamIDs %in% tmhmm_beta_pfams)])
table(Pfam_ConAAC$ancestor[which(Pfam_ConAAC$pfamIDs %in% tmhmm_alpha_pfams)])


trans_chisq_table <- data.frame(c(22,51), c(342,590), row.names = c('post-LUCA', 'modern'))
colnames(trans_chisq_table ) <- c('Beta sheets', 'Alpha helices')
chisq.test(trans_chisq_table)

trans_chisq_table <- data.frame(c(0,22), c(103,342), row.names = c('LUCA', 'post-LUCA'))
colnames(trans_chisq_table ) <- c('Beta sheets', 'Alpha helices')
chisq.test(trans_chisq_table)

trans_chisq_table <- data.frame(c(0,73), c(103,851), row.names = c('LUCA', 'post-LUCA'))
colnames(trans_chisq_table ) <- c('Beta sheets', 'Alpha helices')
chisq.test(trans_chisq_table)


##### Extremophily compositional biases vs ancient usage ####
# prokaryotic GC website https://chrisgaby.github.io/post/prokaryotic-genome-size/index.html
# O2 requirement from https://bacdive.dsmz.de advanced search 
environment_AAC <- read.csv('Environment_AAC.csv', header = T)
environment_AAC$halophily <- gsub('non-halophilic', 'nonHalophilic', environment_AAC$halophily)
prokaryotic_GC <- read.csv('prokaryotes_GC.csv', header = T)
prokaryotic_GC$organisms <- word(prokaryotic_GC$organisms, 1,2)
prokaryotic_GC <- inner_join(prokaryotic_GC ,environment_AAC, by = 'organisms' )
O2_requirement <- read.csv('Bacdive_Oxygen_requirement.csv', header = T)
O2_requirement <- O2_requirement[which(!duplicated(O2_requirement$organisms)),]
environment_AACwithO2 <- inner_join(environment_AAC , O2_requirement, by='organisms')


Trans_ratio <- (colWeightedMeans(as.matrix(Tmhmm_score[which( Tmhmm_score$Tmhmm < 0.5),3:22]),
                                 Tmhmm_score$Conserved_length[which(Tmhmm_score$Tmhmm < 0.5)] )) /
  (colWeightedMeans(as.matrix(Tmhmm_score[which(Tmhmm_score$Tmhmm > 0.5),3:22]),
                    Tmhmm_score$Conserved_length[which(Tmhmm_score$Tmhmm > 0.5)] ))
Trans_ratio <- Trans_ratio[match( names(ancient_clanConusage),names(Trans_ratio))]

GC_ratio <- (colWeightedMeans(as.matrix(prokaryotic_GC[which(prokaryotic_GC$GC. > 60),22:41]),
                              prokaryotic_GC$protein_count[which(prokaryotic_GC$GC. > 60)]))/
  (colWeightedMeans(as.matrix(prokaryotic_GC[which(prokaryotic_GC$GC. < 40),22:41]),
                    prokaryotic_GC$protein_count[which(prokaryotic_GC$GC. < 40)]))

Ph_ratio <-  (colWeightedMeans(as.matrix(environment_AAC[which(environment_AAC$ph_group == 'acidophile'),18:37]),
                               environment_AAC$protein_count[which(environment_AAC$ph_group == 'acidophile')])) /
  (colWeightedMeans(as.matrix(environment_AAC[which(grepl('alkaliphile', environment_AAC$ph_group)),18:37]),
                    environment_AAC$protein_count[which(grepl('alkaliphile', environment_AAC$ph_group))]))

Halo_ratio <- (colWeightedMeans(as.matrix(environment_AAC[which(environment_AAC$halophily == 'nonHalophilic'),18:37]),
                                environment_AAC$protein_count[which(environment_AAC$halophily == 'nonHalophilic')]))/
  (colWeightedMeans(as.matrix(environment_AAC[which(grepl('halophilic|haloalkaliphilic', environment_AAC$halophily)),18:37]),
                    environment_AAC$protein_count[which(grepl('halophilic|haloalkaliphilic', environment_AAC$halophily))]))
Temp_ratio <- (colWeightedMeans(as.matrix(environment_AAC[which(environment_AAC$temp_group == 'mesophile'),18:37]),
                                environment_AAC$protein_count[which(environment_AAC$temp_group == 'mesophile')])) /
  (colWeightedMeans(as.matrix(environment_AAC[which(grepl('thermophile', environment_AAC$temp_group)),18:37]),
                    environment_AAC$protein_count[which(grepl('thermophile', environment_AAC$temp_group))]))

Oxygen_ratio <- (colWeightedMeans(as.matrix(environment_AACwithO2[which( environment_AACwithO2$Oxygen.tolerance == 'aerobe'),18:37]),
                                  environment_AACwithO2$protein_count[which(environment_AACwithO2$Oxygen.tolerance == 'aerobe')])) /
  (colWeightedMeans(as.matrix(environment_AACwithO2[which( environment_AACwithO2$Oxygen.tolerance =='anaerobe'),18:37]),
                    environment_AACwithO2$protein_count[which( environment_AACwithO2$Oxygen.tolerance =='anaerobe')]))

Extremovsusage_df <- cbind(ancient_clanConusage, Ph_ratio, Halo_ratio,Temp_ratio,GC_ratio, Oxygen_ratio, Trans_ratio )
Extremovsusage_df<- as.data.frame(Extremovsusage_df)
rownames(Extremovsusage_df) <- AA_properties$Letter
colnames(Extremovsusage_df) <- c('ancient_clanusage', 'ph', 'salinity', 'temperature','gc', 'oxygen', 'Trans' )
ph_model <- lm(ancient_clanConusage ~ Ph_ratio)
salinity_model <- lm(ancient_clanConusage ~ Halo_ratio)
temp_model <- lm(ancient_clanConusage ~ Temp_ratio)
gc_model <- lm(ancient_clanConusage ~ GC_ratio)
oxygen_model <- lm(ancient_clanConusage ~ Oxygen_ratio)
Trans_model <- lm(ancient_clanConusage ~ Trans_ratio)
summary(Trans_model)

Ph_plot <- ggplot(Extremovsusage_df, aes(y = as.numeric(ancient_clanusage), 
                                         x = as.numeric(Ph_ratio), label = rownames(Extremovsusage_df))) + 
  xlab('Acidophile:Alkaliphile') + ylab(NULL)  + 
  geom_text(size = 8) + theme(legend.position="none") +
  theme(axis.text=element_text(size=14),axis.title=element_text(size=18,face="bold")) + 
  annotate(geom="text", y=0.71, x=0.85, label="p = 0.99", color= "red", size = 8) +
  annotate(geom="text", y=1.2, x=0.63, label="a)", color= "black", size = 10) +
  scale_y_log10(breaks = log_breaks(),limits = c(0.63,1.2)) + scale_x_log10(breaks = log_breaks(), limits = c(0.63,1.2))


salinity_plot <- ggplot(Extremovsusage_df, aes(y = as.numeric(ancient_clanusage), 
                                               x = as.numeric(Halo_ratio), label = rownames(Extremovsusage_df))) + 
  xlab('Non-halophile:Halophile') +  ylab(NULL)  + 
  geom_text(size =8) + theme(legend.position="none") +
  theme(axis.text=element_text(size=14),axis.title=element_text(size=18,face="bold")) + 
  annotate(geom="text", y=0.71, x=0.85, label="p = 0.83", color= "red", size = 8) +
  annotate(geom="text", y=1.2, x=0.63, label="b)", color= "black", size = 10) +
  scale_y_log10(breaks = log_breaks(),limits = c(0.63,1.2)) + scale_x_log10(breaks = log_breaks(), limits = c(0.63,1.2))

temp_plot <- ggplot(Extremovsusage_df, aes(y = as.numeric(ancient_clanusage), 
                                           x = as.numeric(Temp_ratio), label = rownames(Extremovsusage_df))) + 
  xlab('Mesophile:Thermophile') + ylab(NULL)  +
  geom_text(size = 8) + theme(legend.position="none") + 
  theme(axis.text=element_text(size=14),axis.title=element_text(size=18,face="bold")) + 
  annotate(geom="text", y=0.71, x=0.9, label="p = 0.87", color= "red", size = 8) +
  annotate(geom="text", y=1.35, x=0.63, label="c)", color= "black", size = 10) +
  scale_y_log10(breaks = log_breaks(),limits = c(0.63,1.35)) + scale_x_log10(breaks = log_breaks(),limits = c(0.63,1.35))

gc_plot <- ggplot(Extremovsusage_df, aes(y = as.numeric(ancient_clanusage), 
  x = as.numeric(GC_ratio), label = rownames(Extremovsusage_df))) + 
  xlab('High GC:Low GC') + ylab(NULL)  + 
  geom_text(size = 8) + theme(legend.position="none") + 
  theme(axis.text=element_text(size=14),axis.title=element_text(size=18,face="bold")) + 
  annotate(geom="text", y=0.4, x=1.5, label="p = 0.92", color= "red", size = 8) +
  annotate(geom="text", y=7.5, x=0.32, label="e)", color= "black", size = 10) +
  scale_y_log10(breaks = log_breaks(), limits = c(0.31,7.5)) + scale_x_log10(breaks = log_breaks(), limits = c(0.31,7.5))

oxygen_plot <- ggplot(Extremovsusage_df, aes(y = as.numeric(ancient_clanusage), 
  x = as.numeric(Oxygen_ratio), label = rownames(Extremovsusage_df))) + 
  xlab('Aerobe:Anaerobe') + ylab(NULL)  + 
  geom_text(size = 8) + theme(legend.position="none") + 
  theme(axis.text=element_text(size=14),axis.title=element_text(size=18,face="bold")) + 
  annotate(geom="text", y=0.71, x=0.9, label="p = 0.56", color= "red", size = 8) +
  annotate(geom="text", y=1.35, x=0.63, label="d)", color= "black", size = 10) +
  scale_y_log10(breaks = log_breaks(),limits = c(0.63,1.35)) + scale_x_log10(breaks = log_breaks(),limits = c(0.63,1.35))

trans_plot <- ggplot(Extremovsusage_df, aes(y = as.numeric(ancient_clanusage), 
   x = as.numeric(Trans_ratio), label = rownames(Extremovsusage_df))) + 
  xlab('Non-transmembrane:Transmembrane') + ylab(NULL)  + 
  geom_text(size = 8) + theme(legend.position="none") + 
  theme(axis.text=element_text(size=14),axis.title=element_text(size=18,face="bold")) + 
  annotate(geom="text", y=0.4, x=1.5, label="p = 0.91", color= "red", size = 8) +
  annotate(geom="text", y=7.5, x=0.32, label="f)", color= "black", size = 10 ) +
  scale_y_log10(breaks = log_breaks(), limits = c(0.31,7.5)) + scale_x_log10(breaks = log_breaks(), limits = c(0.31,7.5))


## supplementary figure 2
grid.arrange(Ph_plot, salinity_plot, temp_plot, oxygen_plot,gc_plot,trans_plot , nrow =3, ncol =2,
             left = textGrob(expression(bold('Relative ancestral LUCA usage') ),
                             rot = 90,vjust = 1, gp = gpar(cex = 2)))

### Sulfur amino acid analysis ####
environment_AACwithO2$All_environments <- paste(environment_AACwithO2$temp_group ,environment_AACwithO2$sal_group,
                                                environment_AACwithO2$ph_group, environment_AACwithO2$domain,environment_AACwithO2$Oxygen.tolerance,sep='-')
Sulfur_Species <- environment_AACwithO2$organisms[which(grepl('sulf|Sulf|Methano|methano|Thio', environment_AACwithO2$organisms))]
Sulfur_Species <- Sulfur_Species[which(!duplicated(word(Sulfur_Species,1)))]
Nonsulfur_Species <- environment_AACwithO2$organisms[which(!grepl('sulf|Sulf|Methano|methano|Thio', environment_AACwithO2$organisms))]
Nonsulfur_Species <- Nonsulfur_Species[which(!duplicated(word(Nonsulfur_Species,1)))]
environment_AACwithO2_Nonsulfur <- environment_AACwithO2[which(environment_AACwithO2$organisms %in% Nonsulfur_Species),]
environment_AACwithO2_sulfur <- environment_AACwithO2[which(environment_AACwithO2$organisms %in% Sulfur_Species),]


Nonsulfur_Species_df <- environment_AACwithO2[which(environment_AACwithO2$organisms %in% Nonsulfur_Species),]
Nonsulfur_Species_df <- Nonsulfur_Species_df[,c(1,42,19,28)]
Sulfur_Species_df <- environment_AACwithO2[which(environment_AACwithO2$organisms %in% Sulfur_Species),]
Sulfur_Species_df <- Sulfur_Species_df[,c(1,42,19,28)]
sulfur_species_df <- rbind(Nonsulfur_Species_df ,Sulfur_Species_df )
sulfur_species_df$Environment_H2S <- c(rep('Non-H2S-rich', nrow(Nonsulfur_Species_df )),
                                       rep('H2S-rich', nrow(Sulfur_Species_df )))
write.csv(sulfur_species_df, 'sulfur_nonsulfur_species.csv', row.names = F)

## Methionine supplementary figure 4
Nonsulfur_SaafreqM <- as.data.frame(environment_AACwithO2_Nonsulfur  %>% group_by(All_environments) %>% summarise(Nonsulfur_freq = mean(M)))
Sulfur_SaafreqM <- as.data.frame(environment_AACwithO2_sulfur  %>% group_by(All_environments) %>% summarise(Sulfur_freq = mean(M)))
Saafreq_dfM <- merge(Nonsulfur_SaafreqM,Sulfur_SaafreqM  , by='All_environments')
Saafreq_dfM <- Saafreq_dfM[-which(grepl('---', Saafreq_dfM$All_environments)),]
t.test(Saafreq_dfM$Nonsulfur_freq ,Saafreq_dfM$Sulfur_freq ,paired = T)

Mplot_Df <- data.frame(c(Saafreq_dfM$Sulfur_freq,Saafreq_dfM$Nonsulfur_freq), 
                      c(rep(paste0('H', subsc(2),'S-rich'),length(Saafreq_dfM$Sulfur_freq)),
                      rep(paste0('Non-H', subsc(2),'S-rich'),length(Saafreq_dfM$Nonsulfur_fre))),
                      rep(Saafreq_dfM$All_environments,2))
colnames(Mplot_Df) <- c('S_aa_freq', 'Sulfur_environment', 'All_environments')

Methionine_plot <- ggplot(Mplot_Df, aes(y = S_aa_freq, x =Sulfur_environment,color = All_environments)) +
  geom_jitter(width=0) + xlab(NULL) + ylab('Methionine Frequency') + 
  theme(axis.text.y=element_text(size=20),axis.title=element_text(size=30,face="bold"),axis.text.x =element_text(size=28),
        legend.title = element_text(size=26,face="bold"), legend.text = element_text(size=24)) + 
  theme(legend.position = 'none') +
  stat_summary(geom = "point",fun = "mean",col = "black", size = 5,shape = 24,fill = "black") +
  annotate(geom="text", y=0.033, x=1.5, label="Paired t-test p = 0.01", color= "black", size = 12) +
  labs(color='Habitat & Domain', linetype = 'Habitat & Domain') +
  geom_line(aes(x = Sulfur_environment,group = All_environments,
  color = All_environments, linetype = All_environments), linewidth =1.2 ) + scale_shape_manual(values=c(16,rep(16:18, 6))) +
  #scale_linetype_manual(values=c(1,rep(c(1,3,2,4,5,6,7,8,10),each=2))) + 
  scale_linetype_manual(values=c(1:19)) + 
  scale_x_discrete(labels=c(expression("H"[2]*"S-rich"), expression("Non-H"[2]*"S-rich"))) +
  annotate(geom="text", y=0.033, x = 0.5, label="b)", color= "black", size = 14)

## Cysteine supplementary figure 4
Nonsulfur_SaafreqC <- as.data.frame(environment_AACwithO2_Nonsulfur  %>% group_by(All_environments) %>% summarise(Nonsulfur_freq = mean(C)))
Sulfur_SaafreqC <- as.data.frame(environment_AACwithO2_sulfur  %>% group_by(All_environments) %>% summarise(Sulfur_freq = mean(C)))
Saafreq_dfC <- merge(Nonsulfur_SaafreqC,Sulfur_SaafreqC  , by='All_environments')
Saafreq_dfC <- Saafreq_dfC[-which(grepl('---', Saafreq_dfC$All_environments)),]
t.test(Saafreq_dfC$Nonsulfur_freq ,Saafreq_dfC$Sulfur_freq ,paired = T)
Cplot_Df <- data.frame(c(Saafreq_dfC$Sulfur_freq,Saafreq_dfC$Nonsulfur_freq), 
                      c(rep(paste0('H', subsc(2),'S-rich'),length(Saafreq_dfC$Sulfur_freq)),
                        rep(paste0('Non-H', subsc(2),'S-rich'),length(Saafreq_dfC$Nonsulfur_freq))),
                      rep(Saafreq_dfC$All_environments,2))
colnames(Cplot_Df) <- c('S_aa_freq', 'Sulfur_environment', 'All_environments')

## supplementary figure 4
Cysteine_plot <- ggplot(Cplot_Df, aes(y = S_aa_freq, x =Sulfur_environment,color = All_environments)) +
  geom_jitter(width=0) + xlab(NULL) + ylab('Cysteine Frequency') + 
  theme(axis.text.y=element_text(size=20),axis.title=element_text(size=30,face="bold"),axis.text.x =element_text(size=28),
  legend.title = element_text(size=26,face="bold"), legend.text = element_text(size=24)) + 
  theme(legend.position = 'none',legend.title.align=0.5, legend.key.width = unit(4, 'cm')) +
  stat_summary(geom = "point",fun = "mean",col = "black", size = 5,shape = 24,fill = "black") +
  annotate(geom="text", y=0.017, x=1.5, label="Paired t-test p = 0.005", color= "black", size = 12) +
  labs(color='Habitat & Domain', linetype = 'Habitat & Domain') +
  geom_line(aes(x = Sulfur_environment,group = All_environments,
  color = All_environments, linetype = All_environments), linewidth =1.2 ) + scale_shape_manual(values=c(16,rep(16:18, 6))) +
  scale_linetype_manual(values=c(1:19)) +
  scale_x_discrete(labels=c(expression("H"[2]*"S-rich"), expression("Non-H"[2]*"S-rich"))) +
  annotate(geom="text", y=0.017, x = 0.5, label="a)", color= "black", size = 14)

grid.arrange(Cysteine_plot,Methionine_plot, ncol =2)

#####Clustering#####
Rates <- Pfam_ConAAC$avg_clustering_allPhases
phylostrata_age <- vector()
for ( i in 1:nrow(Pfam_ConAAC)) {
  if (grepl('Terrabacteria|PVC|TACK|Proteobacteria|FCB|CPR|DPANN|Euryarchaeota|Asgard',
            Pfam_ConAAC$ancestor[i])) {
    phylostrata_age[i] <- 'modern'
  } else if (grepl('post-LBCA|post-LACA', Pfam_ConAAC$ancestor[i])) {
    phylostrata_age[i] <- 'modern'
  } else if (grepl('LBCA|LACA', Pfam_ConAAC$ancestor[i])) {
    phylostrata_age[i] <- 'postLUCA'
  } else if (grepl('preLUCA', Pfam_ConAAC$ancestor[i])) {
    phylostrata_age[i] <- 'preLUCA'
    # phylostrata_age[i] <- 'LUCA'
  } else if (grepl('LUCA', Pfam_ConAAC$ancestor[i])) {
    phylostrata_age[i] <- 'LUCA'} 
  else { phylostrata_age[i] <- 'unclassifiable'} }
Phylostrata <- phylostrata_age
barplotdata <- data.frame(Phylostrata,Rates)
barplotdata <- barplotdata[-which(barplotdata$Phylostrata == 'unclassifiable'),]

wilcox.test(barplotdata$Rates[which(barplotdata$Phylostrata == 'modern')],
            barplotdata$Rates[which(barplotdata$Phylostrata == 'postLUCA')])

# supplementary figure 1
for (i in 1:nrow(barplotdata))(
  if (barplotdata$Phylostrata[i] == 'preLUCA') {barplotdata$PhylostrataAge[i] <- 1}
  else if (barplotdata$Phylostrata[i] == 'LUCA') {barplotdata$PhylostrataAge[i] <- 2}
  else if (barplotdata$Phylostrata[i] == 'postLUCA') {barplotdata$PhylostrataAge[i] <- 3}
  else if (barplotdata$Phylostrata[i] == 'modern') {barplotdata$PhylostrataAge[i] <- 4})
group_ordered <- with(barplotdata,reorder(Phylostrata,Rates, median))
barplotdata$PhylostrataAge <- factor(barplotdata$PhylostrataAge)

clustering_mean <- barplotdata  %>%   group_by(PhylostrataAge) %>%  
  summarise(average = mean(Rates))  %>% ungroup()
clustering_mean <- clustering_mean[order(clustering_mean$average),]

ggplot(data=barplotdata, aes(x=PhylostrataAge, y=Rates, 
  color=PhylostrataAge)) +   scale_color_brewer(palette="Dark2") +
  geom_violin(trim = T, position = position_dodge(0.9),show.legend = FALSE) +
  theme(axis.text=element_text(size=20),axis.title=element_text(size=25,face="bold")) +
  geom_boxplot(width=0.1, position = position_dodge(0.9), notch = T, aes(middle = mean(Rates)),
   show.legend = FALSE, outlier.shape = NA, coef = 0) +  ylab('Hydrophobic clustering') + xlab(NULL) +
  geom_line(data = clustering_mean,mapping = aes(x = PhylostrataAge, y =average, group=1),
  color="black",linetype=2) + scale_x_discrete(breaks=c("1","2","3",'4'),
  labels=c("pre-LUCA", "LUCA", "post-LUCA", 'modern')) 


