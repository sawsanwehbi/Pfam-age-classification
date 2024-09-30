# Pfam-age-classification
This repository contains the R code used to generate the data in 'Order of amino acid recruitment into the genetic code resolved by Last Universal Common Ancestor’s protein domains'
Input data including reconciled Pfam tree files (Pfam_trees.zip) [1] and prokaryotic Pfam sequences (Prokaryotic_Pfams.zip) [2] can be found on figshare
1. Wehbi, Sawsan (2024). Prokaryotic_Pfams.zip. figshare. Dataset. https://doi.org/10.6084/m9.figshare.25599819.v1
2. Wehbi, Sawsan (2024). Pfam_trees.zip. figshare. Dataset. https://doi.org/10.6084/m9.figshare.25599813.v1

- ClassifiedPFAMs.csv: contains the classifications of 8291 pfams run first through ClassifyingTrees.R then New_LUCAvspreLUCA.R scripts. 9 Pfams are now obsolete so were subsequently removed from the rest of the analysis leaving a total of 8282 pfams. Pfams with an average transfer rate > 0.6 as estimated by GeneRax have been reclassified as 'unclassifiable'
- Pfam_contempAAC.csv and Pfam_data_ancestralAAC.csv contain the Pfams, final classifications, clans, sequence lengths and amino acid frequencies (contemporary and ancestral respectively)
- ASR_AAC.R calculates the pfam ancestral amino acid frequencies from the .state files estimated by IQ-Tree ancestral sequence reconstruction (-asr option)
- AminoAcid_properties.csv contains several amino acid associated metrics complied from various previous publications
- ConservedSites.R Identifies 'conserved' regions in the pfam alignemnts that are ancestrally reconstructed.
- Bacdive_Oxygen_requirement.csv, prokaryotes_GC.csv and Environment_AAC.csv contains information about the oxygen tolerance, GC content and several other environmental conditions of prokaryotes used for the supplementary analysis in supp figures 2 and 4
- DEEPTmhmm_consensusPfam.out , ancestral_DEEPtm_sites.csv and pfam_asr_aac_5seq0.4_DEEPTM.csv correspond to the predicted transmembrane Pfams (using DeepTMHMM), the sites in the pfam alignment that are embdedded in the membrane and their ancestral amino acid frequencies respectively.
- MoodyPfams_probabilities.csv and LUCA_Moody_pathways.csv are files annotating Pfams based on how likely they are to be present in a LUCA full-gene protein classified by Moody et al., 2024 and their corresponding metabolic pathways.
