## code to prepare `realCloneData` dataset goes here

rm(list = ls())
library(ggtree)
library(phangorn)


# Fit the first 55 mutations to a time from conception to 40/52 weeks (0.769 years)
# Only done for embyonic clones with > 30 coalescence times if mutation tree is available

truncation_muts <- 55
utero_time_yrs <- 40/52
embryonicMutRate <- truncation_muts/utero_time_yrs

################################## FABRE #######################################

# Load mutation tree
load("~/Downloads/Clonal_dynamics/Phylogenies/Files/PD41305/trees/tree")
PD41305 <- tree_SNV_c
PD41305 <- drop.tip(PD41305, tip= "Ancestral")

PD41305_truncated <- truncate_tree(PD41305, dist = 55)
plot(PD41305_truncated, direction = "down")
axisPhylo(side=2, backward = F)

fabre_truncated <- list(PD41305_truncated)

fabre_truncated[[1]]$metadata <- data.frame("sourcePaper" = "fabre", "ID" = "PD41305")

names(fabre_truncated) <- "PD41305"




################################# Mitchell ######################################
# Download and unzip mitchell data from https://data.mendeley.com/datasets/np54zjkvxr/1
mitchell_tree_files <- list.files("~/Downloads/trees_clonal_dynamics_hematopoiesis",
                                  pattern = ".tree", recursive = TRUE, full.names = TRUE
)
mitchell_tree_names <- gsub("/.*", "", gsub(".*output_", "", mitchell_tree_files))
mitchell_trees <- lapply(mitchell_tree_files, ape::read.tree)

names(mitchell_trees) <- mitchell_tree_names

## Now import the sensitivity and adjust for that as Mitchell did
# Use Mitchell functions from github (https://github.com/emily-mitchell/normal_haematopoiesis)
source("data-raw/mitchell_fns.R")

adjusted_mut_tree_list <- list()
mitchell <- list()
for (i in 1:length(mitchell_trees)) {
  sensitivity <- read.table(list.files(paste0(
    "~/rotation_fall2021/trees_clonal_dynamics_hematopoiesis/filtering_output_",
    names(mitchell_trees)[i]
  ), pattern = "sensitivity", full.names = T), header = T)
  adjusted_mut_tree <- get_corrected_tree(mitchell_trees[[i]], details = NULL, sensitivity, include_SNVs = TRUE, include_indels = FALSE, get_edge_from_tree = TRUE)
  adjusted_mut_tree_list <- append(adjusted_mut_tree_list, list(adjusted_mut_tree))
}
names(adjusted_mut_tree_list) <- mitchell_tree_names

# Truncate all
mitchell_truncated <- truncate_tree(adjusted_mut_tree_list, dist = 55)

# Don't use truncation for CB001 and CB002
mitchell_truncated[["CB001"]] <- adjusted_mut_tree_list[["CB001"]]
mitchell_truncated[["CB002"]] <- adjusted_mut_tree_list[["CB002"]]

# Add metadata
for(i in 1:length(mitchell_truncated)){
  mitchell_truncated[[i]]$metadata <- data.frame("sourcePaper" = "mitchell",
                                                 "ID" = names(mitchell_truncated)[i])
}


################################ Williams #########################################
williams_raw <- readRDS("~/Downloads/PDD_TELO.rds") # Downloaded from: https://github.com/NickWilliamsSanger/mpn_phylogenies_and_evolution
PD5163 <- readRDS("~/Downloads/example_treedat.RDS")$tree

PD5163 <- ape::drop.tip(PD5163, "zeros", collapse.singles = F)

williams_truncated <- list(truncate_tree(PD5163, dist = 55))

williams_truncated[[1]]$metadata <- data.frame("sourcePaper" = "williams", "ID" = "PD5163")
names(williams_truncated) <- "PD5163"


##### Combined and include time-based trees by scaling ##############################
embryonic_mutation_trees <- c(mitchell_truncated, fabre_truncated, williams_truncated)

embryonic_time_trees <- list()
for(i in 1:length(embryonic_mutation_trees)){
  ultra <- embryonic_mutation_trees[[i]]
  if(names(embryonic_mutation_trees)[i] %in% c("CB001", "CB002")){
    # Use real clone data
    ultra <- cloneRate::realCloneData$fullTrees[[names(embryonic_mutation_trees)[i]]]
  }else{
    ultra <- embryonic_mutation_trees[[i]]
    ultra$edge.length <- ultra$edge.length/embryonicMutRate
  }
  embryonic_time_trees <- c(embryonic_time_trees, list(ultra))
}
names(embryonic_time_trees) <- names(embryonic_mutation_trees)

# Save as loadable data object for package users
usethis::use_data(embryonic_mutation_trees, overwrite = TRUE, compress = "bzip2")
usethis::use_data(embryonic_time_trees, overwrite = TRUE, compress = "bzip2")
