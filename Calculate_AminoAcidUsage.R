# Load the ape package
library(ape)
library(protr)
library(seqinr)
library(broom)
library(ggrepel)
library(stringr)
library(ggplot2)
library(diagis)
library(phytools)
library(dplyr)


setwd("~/Desktop/PFAM Trees")
## Input CSV files
Bacterial_supergroups <- c('CPR', 'PVC', 'FCB', 'Terrabacteria', 'Proteobacteria')
Archaeal_supergroups <- c('Asgard', 'TACK', 'DPANN', 'Euryarchaeota')
All_supergroups <- c(Bacterial_supergroups,Archaeal_supergroups)
Weiss_AAC <- read.csv('AllWeiss_AAC.csv', header = T)
AA_properties <- read.csv('AminoAcid_properties.csv', header = T)
Weiss_clans <- Clans[match( Weiss_AAC$LUCA_Pfams,classified_pfams$PFAM_IDs)]
Weiss_clans <- Weiss_clans[-which(is.na(Weiss_clans))]
classified_pfams <- read.csv('ClassifiedPFAMs.csv', header = T)
PFAM_AAC <- read.csv('Pfam_data.csv', header = T)
Clans_AAC <- read.csv('Clan_data.csv', header = T)

##### ranking phylostrata ####
phylostrata_age <- vector()
for ( i in 1:nrow(classified_pfams)) {
  if (grepl('Terrabacteria|PVC|TACK|Proteobacteria|FCB|CPR|DPANN|Euryarchaeota|Asgard',
            PFAM_AAC$ancestor[i])) {
    phylostrata_age[i] <- 'modern'
  } else if (grepl('post-LBCA|post-LACA', PFAM_AAC$ancestor[i])) {
    phylostrata_age[i] <- 'modern'
  } else if (grepl('LBCA|LACA', PFAM_AAC$ancestor[i])) {
    phylostrata_age[i] <- 'postLUCA'
  } else if (grepl('preLUCA', PFAM_AAC$ancestor[i])) {
    phylostrata_age[i] <- 'preLUCA'
    # phylostrata_age[i] <- 'LUCA'
  } else if (grepl('LUCA', PFAM_AAC$ancestor[i])) {
    phylostrata_age[i] <- 'LUCA'} 
  #else { phylostrata_age[i] <- NA} }
  else { phylostrata_age[i] <- 'unclassifiable'} }
PFAM_AAC <- cbind(PFAM_AAC, phylostrata_age)
group_phylostrata <- PFAM_AAC$ancestor
for ( i in 1:nrow(PFAM_AAC)) {
  if (grepl('LUCA',group_phylostrata[i])) { 
    group_phylostrata[i] <- 'LUCA' } 
  else if (grepl('Inspect',group_phylostrata[i])) {
    group_phylostrata[i] <- 'Inspect' }
  else {
    group_phylostrata[i] <- 'postLUCA' 
  } }
PFAM_AAC <- PFAM_AAC[-which(PFAM_AAC$Clans == "character(0)"),] #obsolete pfams

##### Get Clans AAC ######
# This loop function shows how to classify clans based on their Pfams
# no need to rerun because its already part of the Clan_data.csv file
#unique_Clans <- unique(PFAM_AAC$Clans)
#Clan_ancestor <- vector() 
#for ( i in 1:length(unique_Clans)) {
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
  else {  Clan_ancestor[i] <- 'unclassifiable'}#}

allLUCA_Clans_AAC <- colWeightedMeans(as.matrix(Clans_AAC[which(grepl('LUCA',Clans_AAC$Clan_ancestor)),2:23]),
                                      Clans_AAC$median_ClanLen[which(grepl('LUCA',Clans_AAC$Clan_ancestor))] )
LUCA_Clans_AAC <- colWeightedMeans(as.matrix(Clans_AAC[which(Clans_AAC$Clan_ancestor == 'LUCA'),2:23]),
                                   Clans_AAC$median_ClanLen[which(Clans_AAC$Clan_ancestor == 'LUCA')] )
preLUCA_Clans_AAC <- colWeightedMeans(as.matrix(Clans_AAC[which(Clans_AAC$Clan_ancestor == 'preLUCA'),2:23]),
                                      Clans_AAC$median_ClanLen[which(Clans_AAC$Clan_ancestor == 'preLUCA')] )
postLUCA_Clans_AAC <- colWeightedMeans(as.matrix(Clans_AAC[which(Clans_AAC$Clan_ancestor == 'postLUCA'),2:23]),
                                       Clans_AAC$median_ClanLen[which(Clans_AAC$Clan_ancestor == 'postLUCA')])
modern_Clans_AAC <- colWeightedMeans(as.matrix(Clans_AAC[which(Clans_AAC$Clan_ancestor == 'modern'),2:23]),
                                     Clans_AAC$median_ClanLen[which(Clans_AAC$Clan_ancestor == 'modern')])

ancient_clanusage <- LUCA_Clans_AAC/postLUCA_Clans_AAC 
ancient_clanusage <- ancient_clanusage[match(AA_properties$Letter,names(ancient_clanusage))]
clan_wt <- 1 / (AA_properties$Filtered_sd_2000^2)
clan_model <- lm(AA_properties$Filtered_avg_2000 ~ ancient_clanusage)
clan_wls_model <- lm(AA_properties$Filtered_avg_2000 ~ ancient_clanusage,
                     weights = clan_wt)
summary(clan_wls_model)
clan_wls_confinterval <- broom::augment(clan_wls_model, interval="confidence")
clan_confinterval <- broom::augment(clan_model, interval="confidence")

postLUCAusage <- postLUCA_Clans_AAC/modern_Clans_AAC
postLUCAusage <- postLUCAusage[match(AA_properties$Letter,names(postLUCAusage))]
summary(lm(AA_properties$Filtered_avg_2000 ~ postLUCAusage,weights =clan_wt))

#### standard errors for ratios #####
LUCAweightedse <- vector()
preLUCAweightedse <- vector()
postLUCAweightedse <- vector()
colnames(Clans_AAC)
for (colnb in 1:22) {
  LUCAweightedse[colnb] <- weighted_se(Clans_AAC[,colnb][which(Clans_AAC$Clan_ancestor == 'LUCA')],
                                       Clans_AAC$median_ClanLen[which(Clans_AAC$Clan_ancestor == 'LUCA')])}
for (colnb in 1:22){
  postLUCAweightedse[colnb] <- weighted_se(Clans_AAC[,colnb][which(Clans_AAC$Clan_ancestor == 'postLUCA')],
                                           Clans_AAC$median_ClanLen[which(Clans_AAC$Clan_ancestor == 'postLUCA')])}
for (colnb in 1:22){
  preLUCAweightedse[colnb] <- weighted_se(Clans_AAC[,colnb][which(Clans_AAC$Clan_ancestor == 'preLUCA')],
                                          Clans_AAC$median_ClanLen[which(Clans_AAC$Clan_ancestor == 'preLUCA')])}
LUCAclanratio_var <-  (LUCAweightedse^2)/((postLUCAavgAAC)^2) +
  (postLUCAweightedse^2)*((LUCAavgAAC)^2)/((postLUCAavgAAC)^4)
LUCAclanratio_se <- sqrt(LUCAclanratio_var)
LUCAclanratio_se <- LUCAclanratio_se[match(AA_properties$Letter,names(LUCAclanratio_se))]
preLUCAclanratio_var <-  (preLUCAweightedse^2)/((postLUCAavgAAC)^2) +
  (postLUCAweightedse^2)*((preLUCAavgAAC)^2)/((postLUCAavgAAC)^4)
preLUCAclanratio_se <- sqrt(preLUCAclanratio_var)
preLUCAclanratio_se <- preLUCAclanratio_se[match(AA_properties$Letter,names(preLUCAclanratio_se))]

#### Order plots ####
rm(ordervsusage_df)
ordervsusage_df <- cbind(AA_properties$Filtered_avg_2000 ,AA_properties$Filtered_sd_2000,
                         ancient_clanusage, AA_properties$Moosmann_category)
ordervsusage_df <- as.data.frame(ordervsusage_df)
rownames(ordervsusage_df) <- AA_properties$Letter
colnames(ordervsusage_df) <- c('Trifonov_order', 'Trifonov_se','ancient_clanusage', 'Moosmann_category')

ggplot(ordervsusage_df, aes(y = as.numeric(Trifonov_order), 
                            x = as.numeric(ancient_clanusage ), label = rownames(ordervsusage_df))) + 
  ylab('Averange Rank Order of Recruitment') +
  xlab('Relative amino acid usage in LUCA clans') + 
  geom_text_repel(aes(colour = factor(Moosmann_category)), size = 8, box.padding = 0.1) + 
  labs(color='Moosmann category') + 
  theme(axis.text=element_text(size=15),axis.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=14), legend.title = element_text(size=16,,face="bold")) + 
  scale_fill_brewer(palette="Dark2") +
  theme(legend.position = c(0.15, 0.18)) + xlim(0.77,1.11) +
  geom_line(aes(y = predict(clan_model)), linewidth = 1,color = 'black') +
  annotate(geom="text", y=5, x=0.95, label=paste0("Weighted ",paste0('R',supsc('2')), "= 0.37"), color="black", size = 10) +
  annotate(geom="text", y=4, x=0.95, label="p = 0.002", color="black", size = 10) +
  geom_ribbon(aes(ymin=clan_confinterval$.lower, ymax=clan_confinterval$.upper), colour=NA, alpha=0.3)

### preluca vs luca MODEL 2 regression plot ####
# save pdf as 9.08 and 8.14
error_ratio <- preLUCAclanratio_se/LUCAclanratio_se
demingfit <- deming(preLUCA_clanusage~LUCA_clanusage , weights = 1/error_ratio^2,
                    xstd = LUCAclanratio_se, ystd = preLUCAclanratio_se)
# Use a test model and replace the coefficients with those from the york fit
# to get the confidence interval using the predict function
test_model <- lm(preLUCA_clanusage ~ LUCA_clanusage)
test_model$coefficients <- demingfit$coefficients
CI_data <- predict(test_model,  interval = "confidence", level = 0.95)

ordervsusage_df$Biochem_properties <- c(NA,NA,'charged','charged',NA,NA,NA,NA,'charged',
                NA,'anti-RSS',NA,NA,NA,'charged',NA,NA,NA,'anti-ROS','anti-ROS' )

ggplot(ordervsusage_df, aes(x = as.numeric(LUCAaausage), y = as.numeric(preLUCAaausage), 
                            label = rownames(ordervsusage_df))) + 
  geom_vline(xintercept=1) + geom_hline(yintercept=1) +
  xlab('Relative amino acid usage in LUCA clans') +
  ylab('Relative amino acid usage in preLUCA clans') + 
  geom_text_repel(aes(colour = factor(Biochem_properties)), size = 8, box.padding = 0.05) + 
  labs(color='Biochemical Properties') +   scale_fill_brewer(palette="Dark2") +  
  theme(axis.text=element_text(size=15),axis.title=element_text(size=20,face="bold"),
        legend.text=element_text(size=14), legend.title = element_text(size=16,,face="bold")) + 
  geom_abline(aes(slope = 0.6702535, intercept = 0.3357908), linewidth = 1,color = 'darkblue') +
  annotate(geom="text", y=0.87, x=1.05, label=paste0("",paste0('R',supsc('2')), "= 0.42"), color="black", size = 10) +
  annotate(geom="text", y=0.85, x=1.05, label="p = 0.002", color="black", size = 10) +
  geom_ribbon(aes(ymax=CI_data[,3], ymin=CI_data[,2]), colour=NA, alpha=0.3) +
  geom_abline(aes(slope=1, intercept = 0),linewidth = 1,color = 'darkred')


                             #### Supplementary Analysis  #####
##### Thermophily ######
charged2polar_Clanratio <- (Clans_AAC$D + Clans_AAC$E + Clans_AAC$K + Clans_AAC$R)/
  (Clans_AAC$G + Clans_AAC$H + Clans_AAC$N + Clans_AAC$P +  Clans_AAC$Q + Clans_AAC$S + Clans_AAC$T)

t.test(sqrt(charged2polar_Clanratio[which(Clans_AAC$Clan_ancestor == 'LUCA')]),
       sqrt(charged2polar_Clanratio[which(Clans_AAC$Clan_ancestor == 'preLUCA')]))

t.test(sqrt(charged2polar_Clanratio[which(Clans_AAC$Clan_ancestor == 'preLUCA')]),
       sqrt(charged2polar_Clanratio[which(Clans_AAC$Clan_ancestor == 'postLUCA')]))

hist(sqrt(charged2polar_Clanratio)) # this transform looks good

Rates <- sqrt(charged2polar_Clanratio )
Phylostrata <- Clans_AAC$Clan_ancestor 
barplotdata <- data.frame(Phylostrata,Rates)
barplotdata <- barplotdata[-which(barplotdata$Phylostrata == 'unclassifiable'),]

#####Clustering#####
# Clan df
#Rates <- Clans_AAC$Clans_clustering
#Phylostrata <- Clans_AAC$Clan_ancestor 
# Pfams df
Rates <- PFAM_AAC$avg_clustering_allPhases
Phylostrata <- PFAM_AAC$phylostrata_age
barplotdata <- data.frame(Phylostrata,Rates)
barplotdata <- barplotdata[-which(barplotdata$Phylostrata == 'unclassifiable'),]
wilcox.test(barplotdata$Rates[which(barplotdata$Phylostrata == 'LUCA')],
            barplotdata$Rates[which(barplotdata$Phylostrata == 'modern')])

#### plot supplementary clustering or thermophily figure ####
for (i in 1:nrow(barplotdata))(
  if (barplotdata$Phylostrata[i] == 'preLUCA') {barplotdata$PhylostrataAge[i] <- 1}
  else if (barplotdata$Phylostrata[i] == 'LUCA') {barplotdata$PhylostrataAge[i] <- 2}
  else if (barplotdata$Phylostrata[i] == 'postLUCA') {barplotdata$PhylostrataAge[i] <- 3}
  else if (barplotdata$Phylostrata[i] == 'modern') {barplotdata$PhylostrataAge[i] <- 4})
group_ordered <- with(barplotdata,reorder(Phylostrata,Rates, median))
barplotdata$PhylostrataAge <- factor(barplotdata$PhylostrataAge)
clustering_mean <- barplotdata  %>%  
  group_by(PhylostrataAge) %>%  
  summarize(average = mean(Rates)) %>% 
  ungroup() 
clustering_mean <- clustering_mean[order(clustering_mean$average),]

ggplot(data=barplotdata, aes(x=PhylostrataAge, y=Rates, 
  color=PhylostrataAge)) +   scale_color_brewer(palette="Dark2") +
  geom_violin(trim = T, position = position_dodge(0.9),show.legend = FALSE) +
  theme(axis.text=element_text(size=20),axis.title=element_text(size=25,face="bold")) +
  geom_boxplot(width=0.1, position = position_dodge(0.9), notch = T, aes(middle = mean(Rates)),
  show.legend = FALSE, outlier.shape = NA, coef = 0) +  ylab('Charged to Polar ratio') + xlab(NULL) +
  geom_line(data = clustering_mean,mapping = aes(x = PhylostrataAge, y =average, group=1),
  color="black",linetype=2) + scale_x_discrete(breaks=c("1","2","3",'4'),
  labels=c("pre-LUCA", "LUCA", "post-LUCA", 'modern')) +
  scale_y_continuous(breaks=c(0.5,1,1.5),labels=c(0.25,1,2.25))# back transform tick marks


### SAMPLE 27 PF00133 sequences ####
sample_27 <- c(which(grepl("valine.*TACK|TACK.*valine" , tRNAPFAM_rowsDS$new_seq_names))[1],
               which(grepl("valine.*Euryarchaeota|Euryarchaeota.*valine" , tRNAPFAM_rowsDS$new_seq_names))[1],
               which(grepl("valine.*DPANN|DPANN.*valine" , tRNAPFAM_rowsDS$new_seq_names))[1],
               which(grepl("valine.*Asgard|Asgard.*valine" , tRNAPFAM_rowsDS$new_seq_names))[1],
               which(grepl("valine.*PVC|PVC.*valine" , tRNAPFAM_rowsDS$new_seq_names))[1],
               which(grepl("valine.*FCB|FCB.*valine" , tRNAPFAM_rowsDS$new_seq_names))[1],
               which(grepl("valine.*CPR|CPR.*valine" , tRNAPFAM_rowsDS$new_seq_names))[1],
               which(grepl("valine.*Terrabacteria|Terrabacteria.*valine" , tRNAPFAM_rowsDS$new_seq_names))[1],
               which(grepl("valine.*Proteobacteria|Proteobacteria.*valine" , tRNAPFAM_rowsDS$new_seq_names))[1],
               which(grepl("leucine.*TACK|TACK.*[.]leucine" , tRNAPFAM_rowsDS$new_seq_names))[1],
               which(grepl("leucine.*Euryarchaeota|Euryarchaeota.*[.]leucine" , tRNAPFAM_rowsDS$new_seq_names))[1],
               which(grepl("G001873845_PF00133.DPANN.MAG: leucine--tRNA.Archaea" , tRNAPFAM_rowsDS$new_seq_names))[1],
               which(grepl("leucine.*Asgard|Asgard.*[.]leu" , tRNAPFAM_rowsDS$new_seq_names))[1],
               which(grepl("leucine.*PVC|PVC.*[.]leu" , tRNAPFAM_rowsDS$new_seq_names))[2],
               which(grepl("leucine.*FCB|FCB.*[.]leucine" , tRNAPFAM_rowsDS$new_seq_names))[1],
               which(grepl("leucine.*CPR|CPR.*[.]leucine" , tRNAPFAM_rowsDS$new_seq_names))[1],
               which(grepl("leucine.*Terrabacteria|Terrabacteria.*[.]leucine" , tRNAPFAM_rowsDS$new_seq_names))[1],
               which(grepl("leucine.*Proteobacteria|Proteobacteria.*[.]leucine" , tRNAPFAM_rowsDS$new_seq_names))[1],
               which(grepl("isoleucine.*TACK|TACK.*isoleucine" , tRNAPFAM_rowsDS$new_seq_names))[1],
               which(grepl("isoleucine.*Euryarchaeota|Euryarchaeota.*isoleucine" , tRNAPFAM_rowsDS$new_seq_names))[1],
               which(grepl("DPANN.MAG: isoleucine--tRNA.Archaea" , tRNAPFAM_rowsDS$new_seq_names))[1],
               which(grepl("Asgard.MAG: isoleucyl-tRNA" , tRNAPFAM_rowsDS$new_seq_names)),
               which(grepl("isoleucine.*PVC|PVC.*isoleucine" , tRNAPFAM_rowsDS$new_seq_names))[1],
               which(grepl("isoleucine.*FCB|FCB.*isoleucine" , tRNAPFAM_rowsDS$new_seq_names))[1],
               which(grepl("isoleucine.*CPR|CPR.*isoleucine" , tRNAPFAM_rowsDS$new_seq_names))[1],
               which(grepl("isoleucine.*Terrabacteria|Terrabacteria.*isoleucine" , tRNAPFAM_rowsDS$new_seq_names))[1],
               which(grepl("isoleucine.*Proteobacteria|Proteobacteria.*isoleucine" , tRNAPFAM_rowsDS$new_seq_names))[1])
sample_27 <- sample_27[-which(is.na(sample_27))]
sample_names <- tRNAPFAM_rowsDS$new_seq_names[sample_27]
sample_names <- gsub(' ','', sample_names)
#write.fasta(as.list(tRNAPFAM_rowsDS$pfam_sequence[sample_27]), sample_names ,
        #    'PF00133_27sequences_new.fasta')

