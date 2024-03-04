# Load the needed package
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

## This script takes Pfams that have been robustly (ie > 6 times) classified as LUCA in the iterative classified 'ClassifyingTrees.R'
# and re-classifies them as preLUCA or LUCA
classified_pfams <- read.csv('ClassifiedPFAMs_spr5_NQ.csv', header = T)
Bacterial_supergroups <- c('CPR', 'PVC', 'FCB', 'Terrabacteria', 'Proteobacteria')
Archaeal_supergroups <- c('Asgard', 'TACK', 'DPANN', 'Euryarchaeota')
All_supergroups <- c(Bacterial_supergroups,Archaeal_supergroups)
#PFAM_IDs <- classified_pfams$PFAM_IDs[which(grepl('LUCA',
                  # classified_pfams$robust_classifications6new))]

PFAM_IDs <- classified_pfams$PFAM_IDs[which(grepl('LUCA',
                   classified_pfams$robust_classifications))]
                  
classified_ancestor <- vector()
Ancestor_reconciled_tree <- vector()


# Read in the tree from a file in newick format
for ( i in 1:length(PFAM_IDs)) {
  print(i) 
  reconciled.treefile <- read.tree(paste0('EnforceRoot_trees/',PFAM_IDs[i],'_enforceroot.treefile'))
  
  for (iteration in 1:2) {
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
    if (subclade_mrca_grandparent == getRoot(reconciled.treefile)) {next}
    
    # Step 3: get descendant taxa of parent MRCA node,
    mrca_descendants <- Descendants(reconciled.treefile, subclade_mrca_parent)[[1]]
    parent_mrca_descendants <- Descendants(reconciled.treefile, subclade_mrca_grandparent)[[1]]
    outer_descendants <- parent_mrca_descendants[which(!parent_mrca_descendants %in% subclades[[clade]])]
  
    # Step 5: check for domain monophyly
    # if domain monophyly is false it means that the subclade is nested within the other domain
    bacterial_monophyly <- all(grepl(c('Bacteria'),reconciled.treefile$tip.label[mrca_descendants]))
    archaeal_monophyly <- all(grepl(c('Archaea'),reconciled.treefile$tip.label[mrca_descendants]))
    domain_monophyly <- any(bacterial_monophyly, archaeal_monophyly)
    
    # get the domain of the outer clade to make sure the subclade is indeed the one transferred
    # i.e its from a different domain to its sister clade AND the outer clade
    outer_descendants_supergroups <- unique(word(gsub('[.]', " ",
    str_extract(reconciled.treefile$tip.label[outer_descendants], '[[.]].*')),2))
    
    subclade_in_outergroup <- subclade_supergroup %in% outer_descendants_supergroups
    
    # Step 6: get node label of subclade parent mrca node
    HGT_node.label <- reconciled.treefile$node.label[subclade_mrca_parent - length(reconciled.treefile$tip.label)]
    
    # domain monophyly and accepted edge need to both be false for a subclade to be a
    # confirmed intradomain HGT
    if (domain_monophyly == F & HGT_node.label == 'T' &  subclade_in_outergroup == F ) {
      HGT_tips <- c(HGT_tips, subclades[[clade]]) 
    } else {next}  }  
  
  if (length(HGT_tips) == length(reconciled.treefile$tip.label) |
      length(HGT_tips) == 0 ) {
    reconciled.treefile <- reconciled.treefile } else {
      reconciled.treefile <- drop.tip(reconciled.treefile, HGT_tips) } }

reconciled.treefile <- midpoint(reconciled.treefile) 

  # 2 edges branching directly from root, then the 4 descendants from those root edges
  basal_branches <- which(reconciled.treefile$edge[,1] == getRoot(reconciled.treefile))
  root_edges <- reconciled.treefile$edge[basal_branches,][,2]

  #### Classification function ####
  Ancestor_basal_tree <- vector()
  
  for (basalclade in 1:2) {
    if (root_edges[basalclade] <= length(reconciled.treefile$tip.label)) {next} 
    basal_clade <- extract.clade(reconciled.treefile, root_edges[basalclade])
    basal_branches_clade <- which(basal_clade$edge[,1] == getRoot(basal_clade))
    root_edges_clade <- basal_clade$edge[basal_branches_clade,][,2]
    if (any(root_edges_clade < length(basal_clade$tip.label))) {next}
    root_branches_d1_clade <- which(basal_clade$edge[,1] == root_edges_clade[1])
    root_branches_d2_clade <- which(basal_clade$edge[,1] == root_edges_clade[2])
    basal_branches_clade <- c(basal_branches_clade ,root_branches_d1_clade,root_branches_d2_clade)

    ancient_bacterial_clade <- vector()
    for ( a in 1:length(Bacterial_supergroups)) {
      bacterial_edges <- which.edge(basal_clade, 
                                    basal_clade$tip.label[grepl(Bacterial_supergroups[a],
                                                                basal_clade$tip.label)]) 
     #ancient_bacterial_clade[a] <- any(bacterial_edges %in% basal_branches_clade) }
     ancient_bacterial_clade[a] <- any(bacterial_edges %in% basal_branches_clade) &
     !(all(bacterial_edges %in% basal_branches_clade[1:2]))} # MRCA is NOT the subclade root
  
    ancient_archaeal_clade <- vector()
    for ( b in 1:length(Archaeal_supergroups )) {
      archaeal_edges <- which.edge(basal_clade, 
                                   basal_clade$tip.label[grepl(Archaeal_supergroups[b],
                                                               basal_clade$tip.label)]) 
      #ancient_archaeal_clade[b] <- any(archaeal_edges %in% basal_branches_clade)}
      ancient_archaeal_clade[b] <- any(archaeal_edges %in% basal_branches_clade) &
      !(all(archaeal_edges %in% basal_branches_clade[1:2]))} # MRCA is NOT the subclade root

    if (sum(ancient_bacterial_clade) >=3 & sum(ancient_archaeal_clade) >= 2){
      Ancestor_basal_tree[basalclade] <- 'LUCA'
    } else { Ancestor_basal_tree[basalclade] <- 'modern'} }
  
  if (all(grepl('LUCA', Ancestor_basal_tree))) {
    Ancestor_reconciled_tree <- 'preLUCA' }
  else {Ancestor_reconciled_tree <- 'LUCA'}
  classified_ancestor[i] <- Ancestor_reconciled_tree }

# save in csv file
preLUCAclassified <- classified_pfams$robust_classifications
preLUCAclassified[which(grepl('LUCA',
        classified_pfams$robust_classifications))] <- classified_ancestor
classified_pfams <- cbind(classified_pfams,preLUCAclassified)
write.csv(classified_pfams,'ClassifiedPFAMs_spr5_NQ.csv' , row.names= F)
