## code to prepare `realCloneData` dataset goes here

# For now, just load in the clone list
realCloneData <- readRDS("~/rotation_fall2021/ultrametric_trees_with_clades_as_labeled_all_FOUR_papers_12192022.rds")
setwd("~/Downloads")
system("git clone https://github.com/margaretefabre/Clonal_dynamics")
setwd("Clonal_dynamics/Phylogenies/Files")
system("ls")


usethis::use_data(realCloneData, overwrite = TRUE)
