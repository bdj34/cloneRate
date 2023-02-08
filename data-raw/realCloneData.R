## code to prepare `realCloneData` dataset goes here

# For now, just load in the clone list
#realCloneData <- readRDS("~/rotation_fall2021/ultrametric_trees_with_clades_as_labeled_all_FOUR_papers_12192022.rds")
setwd("~/Downloads")
#system("git clone https://github.com/margaretefabre/Clonal_dynamics")
setwd("Clonal_dynamics/Phylogenies/Files")
system("ls")

#Load PD34493 tree
load("PD41305/trees/tree_ultra")
PD41305_ultra <- tree_SNV_c_ultra
load("PD41276/trees/tree_ultra")
PD41276_ultra <- tree_SNV_c_ultra
load("PD34493/trees/tree_ultra")
PD34493_ultra <- tree_SNV_c_ultra

# Assign ages at sampling (can't find exact ages)
age_vec <- c("PD34493" = 83, "PD41305" = 73, "PD41276" = 79)

#Scale edge length by age to make time-based
PD34493_ultra$edge.length <- PD34493_ultra$edge.length*age_vec["PD34493"]
PD34493_ultra$age <- age_vec[["PD34493"]]

PD41305_ultra$edge.length <- PD41305_ultra$edge.length*age_vec["PD41305"]
PD41305_ultra$age <- age_vec[["PD41305"]]

PD41276_ultra$edge.length <- PD41276_ultra$edge.length*age_vec["PD41276"]
PD41276_ultra$age <- age_vec[["PD41276"]]

fabre_tree_list <- list(PD34493_ultra, PD41305_ultra, PD41276_ultra)
names(fabre_tree_list) <- c("PD34493", "PD41305", "PD41276")

saveRDS(fabre_tree_list, "~/rotation_fall2021/fabre_trees_all.rds")

usethis::use_data(realCloneData, overwrite = TRUE)
