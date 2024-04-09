# set library path to my home directory where packages are installed
myPaths <- .libPaths()
myPaths <- c(myPaths, "/home/u27/sawsanwehbi/R/x86_64-pc-linux-gnu-library/4.0")
myPaths <- c(myPaths[2], myPaths[1])
.libPaths(myPaths)  # add new path
# change accordingly, the steps above are necessary if running on HPC

# Load the needed packages
library(ape)
library(castor)
library(phangorn)
library(data.table)
library(phytools)
library(purrr)
library(stringr)
library(adephylo)
library(stats)
library(dplyr)

# change pathway to your directory of choice
setwd('/groups/masel/sawsanwehbi/Prokaryotic_pfams/')
getwd()

classified_pfams <- read.csv( 'ClassifiedPFAMs.csv', header = T)
PFAM_IDs <- classified_pfams$PFAM_IDs
Bacterial_supergroups <- c('CPR', 'PVC', 'FCB', 'Terrabacteria', 'Proteobacteria')
Archaeal_supergroups <- c('Asgard', 'TACK', 'DPANN', 'Euryarchaeota')
All_supergroups <- c(Bacterial_supergroups,Archaeal_supergroups)
classified_ancestor <- vector()

# Read in the tree from a file in newick format
for ( i in  1:length(PFAM_IDs)) {
  print(i) 
  # enforceroot trees are NQ.pfam estimated (in iqtree), species-tree reconciled (in generax), midpoint rooted (in R) and DTL-optimized (in generax) in this exact order!
  reconciled.treefile <- read.tree(paste0('EnforceRoot_trees/',PFAM_IDs[i],'_enforceroot.treefile'))
  
  ## Loop to identify and prune intradomain hgt taxa, iterated twice to remove nested hgt
  for (iteration in 1:2){
  ## Step 1: getting monophyletic subclades for each supergroup
monophyletic_subclades <- list()
for ( supergroup in 1:length(All_supergroups)) {
  supergroup_tips <- which(grepl(All_supergroups[supergroup],reconciled.treefile$tip.label))
  if (length(supergroup_tips) == 0) {next} else {
    notmono_subclades <- Filter(length, lapply(split(supergroup_tips, 
    cumsum(c(0, diff(supergroup_tips) > 1))), function(x) x[length(x) > 0])) 
    if(length(notmono_subclades) == 0)  {next}  }
  
  # splitting nonmonophyletic subclades into monophyletic
  for (monoclade in 1:length(notmono_subclades)) {
    if ( is.monophyletic(reconciled.treefile, notmono_subclades[[monoclade]])) {
       monosubclades <- list(notmono_subclades[[monoclade]])
       monophyletic_subclades <- c(monophyletic_subclades, monosubclades) }
    else {
      mrca_nonmono_subclade <- reconciled.treefile$edge[which.edge(reconciled.treefile, notmono_subclades[[monoclade]])[1],]
      sep_mono_subclade1 <- Descendants(reconciled.treefile, mrca_nonmono_subclade[2])[[1]][
        Descendants(reconciled.treefile, mrca_nonmono_subclade[2])[[1]] %in% notmono_subclades[[monoclade]]]
      sep_mono_subclade2 <- notmono_subclades[[monoclade]][-which(notmono_subclades[[monoclade]] %in% sep_mono_subclade1)]
      monosubclades <- list(sep_mono_subclade1 ,sep_mono_subclade2)
      monophyletic_subclades <- c(monophyletic_subclades, monosubclades)
    } } }
  

## Step 2: filter out HGT clades based on mixed nodes and node labels from generax
subclades <- monophyletic_subclades 
HGT_tips <- vector()
    
    ## for each subclade
    for ( clade in 1:length(subclades)) {
      if (length(subclades) == 0) {next}
      # get domain of subclade
      subclade_supergroup <- word(gsub('[.]', " ",
      str_extract(reconciled.treefile$tip.label[subclades[[clade]]], '[[.]].*')),2)[1]
      
      # Step 2: get MRCA, parent MRCA and grandparent MRCA nodes for subclade
    if (length(subclades[[clade]]) == 1) 
       #{subclade_mrca <- subclades[[clade]] }
      {next}
      else {subclade_mrca <- getMRCA(reconciled.treefile, subclades[[clade]])}
      if (subclade_mrca == getRoot(reconciled.treefile)) {next} 
      subclade_mrca_parent <- getParent(reconciled.treefile, subclade_mrca)
      if (subclade_mrca_parent == getRoot(reconciled.treefile)) {next} 
      subclade_mrca_grandparent <- getParent(reconciled.treefile, subclade_mrca_parent)
      if (subclade_mrca_grandparent == getRoot(reconciled.treefile)) 
      {next}
        
         # Step 3: get descendant taxa of parent MRCA node,
          mrca_descendants <- Descendants(reconciled.treefile, subclade_mrca_parent)[[1]]
          parent_mrca_descendants <- Descendants(reconciled.treefile, subclade_mrca_grandparent)[[1]]
          sister_clade <-  mrca_descendants[which(!mrca_descendants %in% subclades[[clade]])]
          outer_descendants <- parent_mrca_descendants[which(!parent_mrca_descendants %in% subclades[[clade]])]
          
          # Step 5: check for domain monophyly
          # if domain monophyly is false it means that the subclade is nested within the other domain
          bacterial_monophyly <- all(grepl(c('Bacteria'),reconciled.treefile$tip.label[mrca_descendants]))
          archaeal_monophyly <- all(grepl(c('Archaea'),reconciled.treefile$tip.label[mrca_descendants]))
          domain_monophyly <- any(bacterial_monophyly, archaeal_monophyly)
         
          # get the domain of the outer clade to make sure the subclade is indeed the one transferred
           # i.e its from a different domain to its sister clade AND the outer clade
          #outer_descendants_domains <- unique(word(gsub('[.]', " ",
          #str_extract(reconciled.treefile$tip.label[outer_descendants], '[[.]].*')),3))
          outer_descendants_supergroups <- unique(word(gsub('[.]', " ",
                    str_extract(reconciled.treefile$tip.label[outer_descendants], '[[.]].*')),2))
          
          subclade_in_outergroup <- subclade_supergroup %in% outer_descendants_supergroups
          
          # Step 6: get node label of subclade parent mrca node
          HGT_node.label <- reconciled.treefile$node.label[subclade_mrca_parent - length(reconciled.treefile$tip.label)]
          
          # domain monophyly and accepted edge need to both be false for a subclade to be a
          # confirmed intradomain HGT
             if (domain_monophyly == F & HGT_node.label == 'T' &  subclade_in_outergroup == F ) {
            HGT_tips <- c(HGT_tips, subclades[[clade]])  }
            else {next}  }  

  if (length(HGT_tips) == length(reconciled.treefile$tip.label) |
    length(HGT_tips) == 0 ) {
    reconciled.treefile <- reconciled.treefile } else {
    reconciled.treefile <- drop.tip(reconciled.treefile, HGT_tips) }
}

# midpoint tree 
reconciled.treefile <- midpoint.root(reconciled.treefile)

# iterate 10x while removing most extreme outlier (ie taxa with longest root to tip distance) for robustness
if (length(reconciled.treefile$tip.label) > 20) {
  nbiteration <- 10} else { nbiteration <- 1} # no need to iterate for trees with few taxa

iteration_classified <- vector()
  for (iteration in 1:nbiteration) {
    root2tipdist <- distRoot(reconciled.treefile)
    dropped_tree <- drop.tip(reconciled.treefile , names(which.max(root2tipdist)))
    reconciled.treefile <- midpoint(dropped_tree)

# 2 edges branching directly from root, then the 4 descendants from those root edges
basal_branches <- which(reconciled.treefile$edge[,1] == getRoot(reconciled.treefile))
root_edges <- reconciled.treefile$edge[basal_branches,][,2]
root_branches_d1 <- which(reconciled.treefile$edge[,1] == root_edges[1])
root_branches_d2 <- which(reconciled.treefile$edge[,1] == root_edges[2])
basal_branches <- c(basal_branches ,root_branches_d1,root_branches_d2)

  #### Classification function ####
  ancient_bacterial_clade <- vector()
  for ( a in 1:length(Bacterial_supergroups)) {
    bacterial_edges <- which.edge(reconciled.treefile, 
    reconciled.treefile$tip.label[grepl(Bacterial_supergroups[a],
    reconciled.treefile$tip.label)])
    ancient_bacterial_clade[a] <- any(bacterial_edges %in% basal_branches) &
    !(all(bacterial_edges %in% basal_branches[1:2]))}

  ancient_archaeal_clade <- vector()
  for ( b in 1:length(Archaeal_supergroups )) {
    archaeal_edges <- which.edge(reconciled.treefile, 
    reconciled.treefile$tip.label[grepl(Archaeal_supergroups[b],
    reconciled.treefile$tip.label)])
    ancient_archaeal_clade[b] <- any(archaeal_edges %in% basal_branches) &
    !(all(archaeal_edges %in% basal_branches[1:2]))}

  # classify according to presence of basal clades  
  if (sum(ancient_bacterial_clade) == 5 & sum(ancient_archaeal_clade) == 4){
    Ancestor_reconciled_tree <- 'LUCA'
  } else if (sum(ancient_bacterial_clade) >=3 & sum(ancient_archaeal_clade) >= 2){
    Ancestor_reconciled_tree <- 'LUCA'
  } else if (sum(ancient_bacterial_clade) >=3 & sum(ancient_archaeal_clade) < 2){
    Ancestor_reconciled_tree <- 'LBCA'
  } else if (sum(ancient_bacterial_clade) <3 & sum(ancient_archaeal_clade) >= 2) {
    Ancestor_reconciled_tree <- 'LACA'
  } else if (sum(ancient_bacterial_clade) == 1 & sum(ancient_archaeal_clade) == 0) {
    Ancestor_reconciled_tree <- Bacterial_supergroups[which(ancient_bacterial_clade)]
  } else if (sum(ancient_bacterial_clade) == 1 & sum(ancient_archaeal_clade) == 0) {
    Ancestor_reconciled_tree <- Bacterial_supergroups[which(ancient_bacterial_clade)]
  }  else if (sum(ancient_bacterial_clade) == 2 & sum(ancient_archaeal_clade) == 0) {
    Ancestor_reconciled_tree <- 'post-LBCA'
  } else if (sum(ancient_bacterial_clade) == 0 & sum(ancient_archaeal_clade) == 1) {
    Ancestor_reconciled_tree <- Archaeal_supergroups[which(ancient_archaeal_clade)]
  } else if (sum(ancient_bacterial_clade) == 0 & sum(ancient_archaeal_clade) == 2) {
    Ancestor_reconciled_tree <- 'post-LACA'
  } else {
    Ancestor_reconciled_tree <- 'unclassifiable' } 
  iteration_classified[iteration] <-  Ancestor_reconciled_tree  } 
  #classified_ancestor[i] <- Ancestor_reconciled_tree 
  classified_ancestor <- rbind( classified_ancestor,iteration_classified) } 

iterated_classification_df <- classified_ancestor
rownames(iterated_classification_df) <- PFAM_IDs
write.csv(iterated_classification_df, 'Iteration_classifications')

# Get robust classifications ie > 6 times classifies as the same thing
all_iterations <- read.csv('Iteration_classifications', header = T)
robust_classifications <- vector()
for ( row in 1:nrow(all_iterations)) {
  if ( length(all_iterations[row,] == 1) ) {
  most_abundant <- all_iterations[row,] }
  most_abundant <- names(table(as.character(all_iterations[row,])))[
  which(table(as.character(all_iterations[row,])) > 6)]
  if (length(most_abundant) == 0) { robust_classifications[row] <- 'unclassifiable'}
  else { robust_classifications[row] <- most_abundant}}

  
#depends on the run, name the ancestor vector according to conditions
classified_pfams <- cbind(PFAM_IDs,robust_classifications )
write.csv(classified_pfams, 'ClassifiedPFAMs.csv', row.names = F)

