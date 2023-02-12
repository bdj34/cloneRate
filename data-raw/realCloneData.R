## code to prepare `realCloneData` dataset goes here

rm(list = ls())
library(ggtree)
library(phangorn)

# Set working directory to be the package, if necessary
# setwd("~/package_development/cloneRate/")

# Unfortunately, this requires a lot of manual checking back and forth between
# the published tree figures from the papers and the data. We exclude nested
# expansions which are the result of annotated driver mutations; any subclone
# within a clone with an annotated driver
# that produces more than one tip. These annotations are not
# explicitly included in the data we download, so we have to manually remove them.
# It makes annotating clones cumbersome and slow, but it only has to be done once,
# and we're more confident that it's right.

# Clone git repository into Downloads folder from:
# Fabre et al. 2022 "The longitudinal dynamics and natural history of clonal haematopoiesis"
# https://github.com/margaretefabre/Clonal_dynamics

# Load ultrametric trees
load("~/Downloads/Clonal_dynamics/Phylogenies/Files/PD41305/trees/tree_ultra")
PD41305_ultra <- tree_SNV_c_ultra
load("~/Downloads/Clonal_dynamics/Phylogenies/Files/PD41276/trees/tree_ultra")
PD41276_ultra <- tree_SNV_c_ultra
load("~/Downloads/Clonal_dynamics/Phylogenies/Files/PD34493/trees/tree_ultra")
PD34493_ultra <- tree_SNV_c_ultra

# Assign ages at sampling (can't find exact ages, manually copied from paper)
age_vec <- c("PD34493" = 83, "PD41305" = 73, "PD41276" = 79)

# Scale edge length by age to make time-based in years. Add age to metadata
PD34493_ultra$edge.length <- PD34493_ultra$edge.length * age_vec["PD34493"]
PD34493_ultra$age <- age_vec[["PD34493"]]
PD41305_ultra$edge.length <- PD41305_ultra$edge.length * age_vec["PD41305"]
PD41305_ultra$age <- age_vec[["PD41305"]]
PD41276_ultra$edge.length <- PD41276_ultra$edge.length * age_vec["PD41276"]
PD41276_ultra$age <- age_vec[["PD41276"]]

# Create list with all three main trees and add names as in Fabre paper
fabre <- list(PD34493_ultra, PD41305_ultra, PD41276_ultra)
names(fabre) <- c("PD34493", "PD41305", "PD41276")

# Make a separate list to keep track of the fabre subclones. And a vector naming the subclones
fabre_clones <- list()
fabreCloneNames <- c()

# Edit/remove some of the extraneous info to make data uniform for all papers
for (i in 1:length(fabre)) {
  fabre[[i]]$direction <- NULL
  fabre[[i]]$top <- NULL
  fabre[[i]]$node.label <- NULL
  fabre[[i]]$has_expanded_clades <- F
  fabre[[i]] <- ape::drop.tip(fabre[[i]], "Ancestral", collapse.singles = F)
}

# One by one, annotate clones as they do in the paper, and remove subclones
# First, PD34493
ID <- names(fabre)[1]
tree <- fabre[[ID]]
fabre[[ID]]$has_expanded_clades <- T

### Plot full tree so we can manually identify clones (compare to  Fabre Fig 3B)
ggtree(tree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Get all nodes under 165 (matching SF3B1 clone, in Figure 3B from Fabre)
cloneTips <- c(tree$tip.label[Descendants(tree, 165)[[1]]])
cloneDriver <- "SF3B1:k666n_DelY"
cloneName <- "clone1"
cloneTree <- keep.tip(tree, cloneTips)

# Plot so we can manually identify potential subclones within the clone (none here)
ggtree(cloneTree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# No subclones, add to the list of clones and move on
fabre_clones <- append(fabre_clones, list(cloneTree))
fabreCloneNames <- c(fabreCloneNames, paste0(ID, "_", cloneName, "_", cloneDriver, "_age", round(fabre[[ID]]$age, 2)))

### Plot full tree again so we can manually identify clones (compare again to Fabre Fig. 3B)
ggtree(tree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Node 98 has a clone with a nested subclone from mutation in U2AF1 driver gene
tipsWithNested <- c(tree$tip.label[Descendants(tree, 98)[[1]]])

# We have a nested subclonal expansion from U2AF1. Load annotation from Fabre github
load(paste0(ID, "/trees/details"))
details_U2AF1 <- details3[details3$Gene == "U2AF1", ] # Check data.frame. Shows node 119
# Subtract 1 because we removed outgroup...remove node 118 and all but one descendant
cloneTips <- tipsWithNested[!tipsWithNested %in% c(tree$tip.label[Descendants(tree, 118)[[1]]][-1])]
cloneDriver <- "unknown"
cloneName <- "clone3"
cloneTree <- keep.tip(tree, cloneTips)

# Plot so we can manually identify potential subclones within the clone (none here)
ggtree(cloneTree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# No subclones, add to the list of clones and move on
fabre_clones <- append(fabre_clones, list(cloneTree))
fabreCloneNames <- c(fabreCloneNames, paste0(ID, "_", cloneName, "_", cloneDriver, "_age", round(fabre[[ID]]$age, 2)))

### Plot full tree again so we can manually identify clones (compare again to Fabre Fig. 3B)
ggtree(tree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Get all nodes under 154
cloneTips <- c(tree$tip.label[Descendants(tree, 154)[[1]]])
cloneDriver <- "unknown"
cloneName <- "clone4"
cloneTree <- keep.tip(tree, cloneTips)

# Plot so we can manually identify potential subclones within the clone (none here)
ggtree(cloneTree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# No subclones, add to the list of clones and move on
fabre_clones <- append(fabre_clones, list(cloneTree))
fabreCloneNames <- c(fabreCloneNames, paste0(ID, "_", cloneName, "_", cloneDriver, "_age", round(fabre[[ID]]$age, 2)))


####### PD41305: second individual from Fabre
# Second patient from fabre, PD41305: No known drivers, no nested expansions to worry about
ID <- names(fabre)[2]
tree <- fabre[[ID]]
fabre[[ID]]$has_expanded_clades <- T

### Plot full tree again so we can manually identify clones (compare to Fabre Fig. 3C)
ggtree(tree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Get nodes under node 143
cloneTips <- c(tree$tip.label[Descendants(tree, 143)[[1]]])
cloneDriver <- "unknown"
cloneName <- "clone1"
cloneTree <- keep.tip(tree, cloneTips)

# Plot so we can manually identify potential subclones within the clone (none here)
ggtree(cloneTree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# No subclones, add to the list of clones and move on
fabre_clones <- append(fabre_clones, list(cloneTree))
fabreCloneNames <- c(fabreCloneNames, paste0(ID, "_", cloneName, "_", cloneDriver, "_age", round(fabre[[ID]]$age, 2)))

### Plot full tree again so we can manually identify clones (compare to Fabre Fig. 3C)
ggtree(tree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Tips under node 102 are clone 2
cloneTips <- c(tree$tip.label[Descendants(tree, 102)[[1]]])
cloneDriver <- "unknown"
cloneName <- "clone2"
cloneTree <- keep.tip(tree, cloneTips)

# Plot so we can manually identify potential subclones within the clone (none here)
ggtree(cloneTree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# No subclones, add to the list of clones and move on
fabre_clones <- append(fabre_clones, list(cloneTree))
fabreCloneNames <- c(fabreCloneNames, paste0(ID, "_", cloneName, "_", cloneDriver, "_age", round(fabre[[ID]]$age, 2)))

### Plot full tree again so we can manually identify clones (compare to Fabre Fig. 3C)
ggtree(tree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# tips under node 119 are clone 3
cloneTips <- c(tree$tip.label[Descendants(tree, 119)[[1]]])
cloneDriver <- "unknown"
cloneName <- "clone3"
cloneTree <- keep.tip(tree, cloneTips)

# Plot so we can manually identify potential subclones within the clone (none here)
ggtree(cloneTree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# No subclones, add to the list of clones and move on
fabre_clones <- append(fabre_clones, list(cloneTree))
fabreCloneNames <- c(fabreCloneNames, paste0(ID, "_", cloneName, "_", cloneDriver, "_age", round(fabre[[ID]]$age, 2)))



####### PD41276: third individual from Fabre (also typo'ed as PD42176 in paper)
# 1 driver with small nested expansions annotated
ID <- names(fabre)[3]
tree <- fabre[[ID]]
fabre[[ID]]$has_expanded_clades <- T

### Plot full tree again so we can manually identify clones (compare to Fabre Fig. 3A)
ggtree(tree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Tips under 86 are the clone, but we must remove nested expansions
tipsWithNested <- c(tree$tip.label[Descendants(tree, 86)[[1]]])
# Load the details to find nested drivers
load("~/rotation_fall2021/fabre_trees/details_41276")
details_CBL <- details3[details3$Gene == "CBL", ] # Shows node 136. Subtract 1 because we removed outgroup...135
details_del20q <- details3[details3$variant_ID == "del20q", ] # Shows node 165. Subtract 1 because we removed outgroup...164
cloneTips <- tipsWithNested[!tipsWithNested %in% c(
  tree$tip.label[Descendants(tree, 135)[[1]]][-1],
  tree$tip.label[Descendants(tree, 164)[[1]]][-1]
)]
cloneDriver <- "SF3B1:k666n"
cloneName <- "clone1"
cloneTree <- keep.tip(tree, cloneTips)

# Plot so we can manually identify potential subclones within the clone (already removed here)
ggtree(cloneTree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# No subclones, add to the list of clones and move on
fabre_clones <- append(fabre_clones, list(cloneTree))
fabreCloneNames <- c(fabreCloneNames, paste0(ID, "_", cloneName, "_", cloneDriver, "_age", round(fabre[[ID]]$age, 2)))

# Name the fabre_clones list
names(fabre_clones) <- fabreCloneNames









############# Download Mitchell data
# Download and unzip mitchell data from https://data.mendeley.com/datasets/np54zjkvxr/1
mitchell_tree_files <- list.files("~/Downloads/trees_clonal_dynamics_hematopoiesis",
  pattern = ".tree", recursive = TRUE, full.names = TRUE
)
mitchell_tree_names <- gsub("/.*", "", gsub(".*output_", "", mitchell_tree_files))
mitchell_trees <- lapply(mitchell_tree_files, ape::read.tree)

names(mitchell_trees) <- mitchell_tree_names

# age vec must be in order of mitchell_trees
age_vec <- c(63, 1, 1, 29, 38, 81, 78, 75, 76, 48)

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

  agedf <- data.frame("age" = age_vec[i], "tip.label" = adjusted_mut_tree$tip.label)

  adjusted_mut_tree$agedf <- agedf

  adjusted_mut_tree$edge.length <- round(adjusted_mut_tree$edge.length)

  tree_ultra_length_one <- make.ultrametric.tree(adjusted_mut_tree)
  ultra_tree <- tree_ultra_length_one
  ultra_tree$edge.length <- ultra_tree$edge.length * age_vec[i]
  mitchell <- append(mitchell, list(ultra_tree))
}

names(adjusted_mut_tree_list) <- names(mitchell_trees)
names(mitchell) <- names(mitchell_trees)

for (i in 1:length(mitchell)) {
  mitchell[[i]]$has_expanded_clades <- F
  mitchell[[i]]$age <- mitchell[[i]]$agedf$age[1]
  mitchell[[i]]$agedf <- NULL
}

# Start a list with clones and a vector of names
mitchell_clones <- list()
mitchellCloneNames <- c()

### First patient: KX003: Only one nested expansion in one of the clones we track: TET2 in largest clone
# See Fig. 3C in Mitchell et al.
ID <- "KX003"
tree <- mitchell[[ID]]
mitchell[[ID]]$has_expanded_clades <- T

### Plot full tree again so we can manually identify clones (compare to Mitchell Fig. 3C)
ggtree(tree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Tips under 342 are a clone, but we have to remove nested expansion from TET2
tipsWithNested <- c(tree$tip.label[Descendants(tree, 342)[[1]]])
cloneDriver <- "unknown"
cloneName <- "clone1" # Not labelled in paper but label from left to right from figure
cloneTree <- keep.tip(tree, tipsWithNested)

# Plot to visually inspect and find/remove subclones
ggtree(cloneTree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Remove tips which are descendants of 195 (by visual inspection of subtree compared to annotated tree in fig3C)
cloneTips <- tipsWithNested[!tipsWithNested %in% c(cloneTree$tip.label[Descendants(cloneTree, 195)[[1]]][-1])]
cloneTree <- keep.tip(cloneTree, cloneTips)

# Add to the list of clones and move on
mitchell_clones <- append(mitchell_clones, list(cloneTree))
mitchellCloneNames <- c(mitchellCloneNames, paste0(ID, "_", cloneName, "_", cloneDriver, "_age", round(mitchell[[ID]]$age, 2)))

### Plot full tree again so we can manually identify clones (compare to Mitchell Fig. 3C)
ggtree(tree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

cloneTips <- c(tree$tip.label[Descendants(tree, 453)[[1]]])
cloneDriver <- "unknown"
cloneName <- "clone2"
cloneTree <- keep.tip(tree, cloneTips)

# Plot to check
ggtree(cloneTree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Add to the list of clones and move on
mitchell_clones <- append(mitchell_clones, list(cloneTree))
mitchellCloneNames <- c(mitchellCloneNames, paste0(ID, "_", cloneName, "_", cloneDriver, "_age", round(mitchell[[ID]]$age, 2)))


# Second patient: KX004: many nested expansions
ID <- "KX004"
tree <- mitchell[[ID]]
mitchell[[ID]]$has_expanded_clades <- T

### Plot full tree again so we can manually identify clones (compare to Mitchell Fig. 3B)
ggtree(tree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# First clone is node 462, must remove 2p_CN_LOH
tipsWithNested <- c(tree$tip.label[Descendants(tree, 462)[[1]]])
cloneDriver <- "DNMT3A"
cloneName <- "clone1"
cloneTree <- keep.tip(tree, tipsWithNested)

# Plot to check
ggtree(cloneTree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Remove tips which are descendants of 154 (by visual inspection of subtree compared to paper annotation)
cloneTips <- tipsWithNested[!tipsWithNested %in% c(cloneTree$tip.label[Descendants(cloneTree, 154)[[1]]][-1])]
cloneTree <- keep.tip(cloneTree, cloneTips)

# Add to the list of clones and move on
mitchell_clones <- append(mitchell_clones, list(cloneTree))
mitchellCloneNames <- c(mitchellCloneNames, paste0(ID, "_", cloneName, "_", cloneDriver, "_age", round(mitchell[[ID]]$age, 2)))


### Plot full tree again so we can manually identify clones (compare to Mitchell Fig. 3B)
ggtree(tree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Tips under node 542 are clone, by visual comparison to Fig 3B
cloneTips <- c(tree$tip.label[Descendants(tree, 542)[[1]]]) # No nested
cloneDriver <- "DNMT3A"
cloneName <- "clone2"
cloneTree <- keep.tip(tree, cloneTips)

# Plot to check (no nested)
ggtree(cloneTree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Add to the list of clones and move on
mitchell_clones <- append(mitchell_clones, list(cloneTree))
mitchellCloneNames <- c(mitchellCloneNames, paste0(ID, "_", cloneName, "_", cloneDriver, "_age", round(mitchell[[ID]]$age, 2)))


### Plot full tree again so we can manually identify clones (compare to Mitchell Fig. 3B)
ggtree(tree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

cloneTips <- c(tree$tip.label[Descendants(tree, 580)[[1]]]) # No nested
cloneDriver <- "unknown"
cloneName <- "clone3"
cloneTree <- keep.tip(tree, cloneTips)

# Plot to check (no nested)
ggtree(cloneTree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Add to the list of clones and move on
mitchell_clones <- append(mitchell_clones, list(cloneTree))
mitchellCloneNames <- c(mitchellCloneNames, paste0(ID, "_", cloneName, "_", cloneDriver, "_age", round(mitchell[[ID]]$age, 2)))


### Plot full tree again so we can manually identify clones (compare to Mitchell Fig. 3B)
ggtree(tree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Tips under node 620, by visual inspection
cloneTips <- c(tree$tip.label[Descendants(tree, 620)[[1]]]) # No nested
cloneDriver <- "DNMT3A"
cloneName <- "clone4"
cloneTree <- keep.tip(tree, cloneTips)

# Plot to check (no nested)
ggtree(cloneTree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Add to the list of clones and move on
mitchell_clones <- append(mitchell_clones, list(cloneTree))
mitchellCloneNames <- c(mitchellCloneNames, paste0(ID, "_", cloneName, "_", cloneDriver, "_age", round(mitchell[[ID]]$age, 2)))


### Plot full tree again so we can manually identify clones (compare to Mitchell Fig. 3B)
ggtree(tree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Tips under node 659, by visual inspection
cloneTips <- c(tree$tip.label[Descendants(tree, 659)[[1]]]) # No nested
cloneDriver <- "unknown"
cloneName <- "clone5"
cloneTree <- keep.tip(tree, cloneTips)

# Plot to check (no nested)
ggtree(cloneTree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Add to the list of clones and move on
mitchell_clones <- append(mitchell_clones, list(cloneTree))
mitchellCloneNames <- c(mitchellCloneNames, paste0(ID, "_", cloneName, "_", cloneDriver, "_age", round(mitchell[[ID]]$age, 2)))


### Plot full tree again so we can manually identify clones (compare to Mitchell Fig. 3B)
ggtree(tree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Tips under 699, by visual inspection
cloneTips <- c(tree$tip.label[Descendants(tree, 699)[[1]]]) # No nested
cloneDriver <- "DNMT3A"
cloneName <- "clone"
cloneTree <- keep.tip(tree, cloneTips)

# Plot to check (no nested)
ggtree(cloneTree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Add to the list of clones and move on
mitchell_clones <- append(mitchell_clones, list(cloneTree))
mitchellCloneNames <- c(mitchellCloneNames, paste0(ID, "_", cloneName, "_", cloneDriver, "_age", round(mitchell[[ID]]$age, 2)))


### Plot full tree again so we can manually identify clones (compare to Mitchell Fig. 3B)
ggtree(tree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Tips under 771, by visual inspection
cloneTips <- c(tree$tip.label[Descendants(tree, 771)[[1]]]) # No nested
cloneDriver <- "unknown"
cloneName <- "clone7"
cloneTree <- keep.tip(tree, cloneTips)

# Plot to check (no nested)
ggtree(cloneTree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Add to the list of clones and move on
mitchell_clones <- append(mitchell_clones, list(cloneTree))
mitchellCloneNames <- c(mitchellCloneNames, paste0(ID, "_", cloneName, "_", cloneDriver, "_age", round(mitchell[[ID]]$age, 2)))


### Plot full tree again so we can manually identify clones (compare to Mitchell Fig. 3B)
ggtree(tree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Tips under 867, by visual inspection
cloneTips <- c(tree$tip.label[Descendants(tree, 867)[[1]]]) # No nested
cloneDriver <- "unknown"
cloneName <- "clone8"
cloneTree <- keep.tip(tree, cloneTips)

# Plot to check (no nested)
ggtree(cloneTree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Add to the list of clones and move on
mitchell_clones <- append(mitchell_clones, list(cloneTree))
mitchellCloneNames <- c(mitchellCloneNames, paste0(ID, "_", cloneName, "_", cloneDriver, "_age", round(mitchell[[ID]]$age, 2)))


# KX007 (no expansions with n >= 10 tips, once nested expansions removed)
# Move on to KX008:
ID <- "KX008"
tree <- mitchell[[ID]]
mitchell[[ID]]$has_expanded_clades <- T

### Plot full tree again so we can manually identify clones (compare to Mitchell Fig. 3A)
ggtree(tree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Tips under node 482, by visual inspection
cloneTips <- c(tree$tip.label[Descendants(tree, 482)[[1]]])
cloneDriver <- "unknown"
cloneName <- "clone1"
cloneTree <- keep.tip(tree, cloneTips)

# Plot to check (no nested)
ggtree(cloneTree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Add to the list of clones and move on
mitchell_clones <- append(mitchell_clones, list(cloneTree))
mitchellCloneNames <- c(mitchellCloneNames, paste0(ID, "_", cloneName, "_", cloneDriver, "_age", round(mitchell[[ID]]$age, 2)))


### Plot full tree again so we can manually identify clones (compare to Mitchell Fig. 3A)
ggtree(tree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Tips under node 549, by visual inspection
cloneTips <- c(tree$tip.label[Descendants(tree, 549)[[1]]])
cloneDriver <- "unknown"
cloneName <- "clone2"
cloneTree <- keep.tip(tree, cloneTips)

# Plot to check (no nested)
ggtree(cloneTree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Add to the list of clones and move on
mitchell_clones <- append(mitchell_clones, list(cloneTree))
mitchellCloneNames <- c(mitchellCloneNames, paste0(ID, "_", cloneName, "_", cloneDriver, "_age", round(mitchell[[ID]]$age, 2)))


### Plot full tree again so we can manually identify clones (compare to Mitchell Fig. 3A)
ggtree(tree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Tips under node 598, by visual inspection
cloneTips <- c(tree$tip.label[Descendants(tree, 598)[[1]]])
cloneDriver <- "unknown"
cloneName <- "clone3"
cloneTree <- keep.tip(tree, cloneTips)

# Plot to check (no nested)
ggtree(cloneTree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Add to the list of clones and move on
mitchell_clones <- append(mitchell_clones, list(cloneTree))
mitchellCloneNames <- c(mitchellCloneNames, paste0(ID, "_", cloneName, "_", cloneDriver, "_age", round(mitchell[[ID]]$age, 2)))


### Plot full tree again so we can manually identify clones (compare to Mitchell Fig. 3A)
ggtree(tree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Tips under node 644, by visual inspection
cloneTips <- c(tree$tip.label[Descendants(tree, 644)[[1]]])
cloneDriver <- "unknown"
cloneName <- "clone4"
cloneTree <- keep.tip(tree, cloneTips)

# Plot to check (no nested)
ggtree(cloneTree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Add to the list of clones and move on
mitchell_clones <- append(mitchell_clones, list(cloneTree))
mitchellCloneNames <- c(mitchellCloneNames, paste0(ID, "_", cloneName, "_", cloneDriver, "_age", round(mitchell[[ID]]$age, 2)))


### Plot full tree again so we can manually identify clones (compare to Mitchell Fig. 3A)
ggtree(tree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

cloneTips <- c(tree$tip.label[Descendants(tree, 703)[[1]]])
cloneDriver <- "unknown"
cloneName <- "clone5"
cloneTree <- keep.tip(tree, cloneTips)

# Plot to check (no nested)
ggtree(cloneTree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Add to the list of clones and move on
mitchell_clones <- append(mitchell_clones, list(cloneTree))
mitchellCloneNames <- c(mitchellCloneNames, paste0(ID, "_", cloneName, "_", cloneDriver, "_age", round(mitchell[[ID]]$age, 2)))

# Name the mitchell_clones list
names(mitchell_clones) <- mitchellCloneNames










############# Download Williams data from https://github.com/NickWilliamsSanger/mpn_phylogenies_and_evolution
# Williams
williams_raw <- readRDS("~/Downloads/PDD_TELO.rds") # Downloaded from: https://github.com/NickWilliamsSanger/mpn_phylogenies_and_evolution
williams <- list()
names_vec <- c()

# Split trees by age and remove zeroes outgroup
for (i in 1:length(williams_raw)) {
  # Combined object may have more than one timepoint
  combined_object <- williams_raw[[i]]
  patientName <- names(williams_raw)[i]

  # Remove zeroes group but keep root
  combined_object$ultratree <- drop.tip(combined_object$ultratree, "zeros", collapse.singles = F)

  # Get agedf with driver metadata for all timepoints
  combined_object$agedf <- combined_object$agedf[combined_object$agedf$age_at_sample_pcy > 1, ]
  unique_ages <- sort(unique(combined_object$agedf$age_at_sample_pcy))
  agedf_all <- combined_object$agedf

  for (j in (1:length(unique_ages))) {
    # Subset separate_object to only contain a single timepoint
    separate_object <- combined_object$ultratree
    age <- unique_ages[j]
    tipsDrop <- combined_object$agedf$tip.label[combined_object$agedf$age_at_sample_pcy != age]
    separate_object <- drop.tip(separate_object, tipsDrop, collapse.singles = T)
    separate_object$age <- age

    # Change data format
    separate_object$edge.length <- as.numeric(separate_object$edge.length)
    # has_expanded_clade set to F until we find a clade making it T
    separate_object$has_expanded_clades <- F

    # Cut agedf to only have metadata from a single timepoint
    agedf <- agedf_all[!agedf_all$tip.label %in% tipsDrop, ]

    separate_object$agedf <- NULL

    williams <- append(williams, list(separate_object))
    names_vec <- c(names_vec, paste0(names(williams_raw)[i], "_", j))
  }
}

names(williams) <- names_vec

# Start list of clones and names
williams_clones <- list()
williamsCloneNames <- c()



### Go one by one, starting with PD7271: No nested
ID <- "PD7271_1"
tree <- williams[[ID]]
williams[[ID]]$has_expanded_clades <- T

### Plot full tree again so we can manually identify clones (compare to Williams Fig. 2)
ggtree(tree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Tips under 120 are a clone
cloneTips <- c(tree$tip.label[Descendants(tree, 120)[[1]]])
cloneDriver <- "JAK2_v617f"
cloneName <- "clone1"
cloneTree <- keep.tip(tree, cloneTips)

# Plot to check (no nested)
ggtree(cloneTree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Add to the list of clones and move on
williams_clones <- append(williams_clones, list(cloneTree))
williamsCloneNames <- c(williamsCloneNames, paste0(ID, "_", cloneName, "_", cloneDriver, "_age", round(williams[[ID]]$age, 2)))



# PD5182_1: No nested. Note: Number after underscore indicates the timepoint (some individuals have multiple sampling times)
ID <- "PD5182_1"
tree <- williams[[ID]]
williams[[ID]]$has_expanded_clades <- T

### Plot full tree again so we can manually identify clones (compare to Williams Fig 3)
ggtree(tree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Tips under node 38 are clonal
cloneTips <- c(tree$tip.label[Descendants(tree, 38)[[1]]])
cloneDriver <- "JAK2_v617f_AND_9pUPD"
cloneName <- "clone1"
cloneTree <- keep.tip(tree, cloneTips)

# Plot to check (no nested)
ggtree(cloneTree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Add to the list of clones and move on
williams_clones <- append(williams_clones, list(cloneTree))
williamsCloneNames <- c(williamsCloneNames, paste0(ID, "_", cloneName, "_", cloneDriver, "_age", round(williams[[ID]]$age, 2)))


# PD5182_2: No clades!
ID <- "PD5182_2"
tree <- williams[[ID]]
### Plot full tree again so we can manually identify clones (compare to Williams Fig 3)
ggtree(tree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()
williams[[ID]]$has_expanded_clades <- F


# PD5182_3: No nested
ID <- "PD5182_3"
tree <- williams[[ID]]
williams[[ID]]$has_expanded_clades <- T

### Plot full tree again so we can manually identify clones (compare to Williams Fig 3)
ggtree(tree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Tips under node 102
cloneTips <- c(tree$tip.label[Descendants(tree, 102)[[1]]])
cloneDriver <- "JAK2_v617f_AND_9pUPD_AND_1q+"
cloneName <- "clone2" # Don't match the other clone name because they aren't the same
cloneTree <- keep.tip(tree, cloneTips)

# Plot to check (no nested)
ggtree(cloneTree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Add to the list of clones and move on
williams_clones <- append(williams_clones, list(cloneTree))
williamsCloneNames <- c(williamsCloneNames, paste0(ID, "_", cloneName, "_", cloneDriver, "_age", round(williams[[ID]]$age, 2)))



# PD5847_1: No nested (1 DMNT3A but only 1 doesn't affect coal times)
ID <- "PD5847_1"
tree <- williams[[ID]]
williams[[ID]]$has_expanded_clades <- T

### Plot full tree again so we can manually identify clones (compare to Williams Fig 3)
ggtree(tree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Tips under node 103
cloneTips <- c(tree$tip.label[Descendants(tree, 103)[[1]]])
cloneDriver <- "JAK2:v617f_AND_9pUPD_AND_TET2:p.N281fs*1"
cloneName <- "clone1"
cloneTree <- keep.tip(tree, cloneTips)

# Plot to check (no nested)
ggtree(cloneTree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Add to the list of clones and move on
williams_clones <- append(williams_clones, list(cloneTree))
williamsCloneNames <- c(williamsCloneNames, paste0(ID, "_", cloneName, "_", cloneDriver, "_age", round(williams[[ID]]$age, 2)))



# PD5179_1: No nested (1 annotated (CUX1) but doesn't split)
ID <- "PD5179_1"
tree <- williams[[ID]]
williams[[ID]]$has_expanded_clades <- T

### Plot full tree again so we can manually identify clones (compare to Williams Fig 3)
ggtree(tree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

cloneTips <- c(tree$tip.label[Descendants(tree, 97)[[1]]])
cloneDriver <- "JAK2:v617f_AND_9+_AND_9q-_AND_1q+"
cloneName <- "clone1"
cloneTree <- keep.tip(tree, cloneTips)

# Plot to check (no nested)
ggtree(cloneTree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Add to the list of clones and move on
williams_clones <- append(williams_clones, list(cloneTree))
williamsCloneNames <- c(williamsCloneNames, paste0(ID, "_", cloneName, "_", cloneDriver, "_age", round(williams[[ID]]$age, 2)))


# PD6634_1: No nested
ID <- "PD6634_1"
tree <- williams[[ID]]
williams[[ID]]$has_expanded_clades <- T

### Plot full tree again so we can manually identify clones (compare to Williams Fig 4)
ggtree(tree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Tips under node 39
cloneTips <- c(tree$tip.label[Descendants(tree, 39)[[1]]])
cloneDriver <- "PPM1D:p.Q462*"
cloneName <- "clone1"
cloneTree <- keep.tip(tree, cloneTips)

# Plot to check (no nested)
ggtree(cloneTree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Add to the list of clones and move on
williams_clones <- append(williams_clones, list(cloneTree))
williamsCloneNames <- c(williamsCloneNames, paste0(ID, "_", cloneName, "_", cloneDriver, "_age", round(williams[[ID]]$age, 2)))



# PD9478_1: No nested
ID <- "PD9478_1"
tree <- williams[[ID]]
williams[[ID]]$has_expanded_clades <- T

### Plot full tree again so we can manually identify clones (compare to Williams Fig 3)
ggtree(tree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Tips under node 85 are a clone
cloneTips <- c(tree$tip.label[Descendants(tree, 85)[[1]]])
cloneDriver <- "JAK2:p.F537_K539delinsL_AND_DNMT3A:p.Y908*"
cloneName <- "clone1"
cloneTree <- keep.tip(tree, cloneTips)

# Plot to check (no nested)
ggtree(cloneTree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Add to the list of clones and move on
williams_clones <- append(williams_clones, list(cloneTree))
williamsCloneNames <- c(williamsCloneNames, paste0(ID, "_", cloneName, "_", cloneDriver, "_age", round(williams[[ID]]$age, 2)))



# PD6629_1: One nested (JAK2 within DMNT3A)
ID <- "PD6629_1"
tree <- williams[[ID]]
williams[[ID]]$has_expanded_clades <- T

### Plot full tree again so we can manually identify clones (compare to Williams Fig 3)
ggtree(tree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Tips under node 34 but must remove the JAK2 subclone (node 39)
tipsWithNested <- c(tree$tip.label[Descendants(tree, 34)[[1]]])
cloneTips <- tipsWithNested[!tipsWithNested %in% tree$tip.label[Descendants(tree, 39)[[1]]][-1]]
cloneDriver <- "DNMT3A:p.R882H"
cloneName <- "clone1"
cloneTree <- keep.tip(tree, cloneTips)

# Plot to check (no nested)
ggtree(cloneTree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Add to the list of clones and move on
williams_clones <- append(williams_clones, list(cloneTree))
williamsCloneNames <- c(williamsCloneNames, paste0(ID, "_", cloneName, "_", cloneDriver, "_age", round(williams[[ID]]$age, 2)))


### Plot full tree again so we can manually identify clones (compare to Williams Fig 3)
ggtree(tree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Tips under node 39 are the JAK2 clone
cloneTips <- c(tree$tip.label[Descendants(tree, 39)[[1]]])
# Keep single nested tip with TET2 because coal time unaffected
cloneDriver <- "DNMT3A:p.R882H_AND_JAK2:p.V617F"
cloneName <- "clone2"
cloneTree <- keep.tip(tree, cloneTips)

# Plot to check (no nested)
ggtree(cloneTree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Add to the list of clones and move on
williams_clones <- append(williams_clones, list(cloneTree))
williamsCloneNames <- c(williamsCloneNames, paste0(ID, "_", cloneName, "_", cloneDriver, "_age", round(williams[[ID]]$age, 2)))



# PD6629_2: Nested expansions to be removed
ID <- "PD6629_2"
tree <- williams[[ID]]
williams[[ID]]$has_expanded_clades <- T

### Plot full tree again so we can manually identify clones (compare to Williams Fig 3)
ggtree(tree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Tips under 31 but must remove subclones
tipsWithNested <- c(tree$tip.label[Descendants(tree, 31)[[1]]])
# Remove the nested expansion with TET2 (node 38) and 9pUPD (node 46)
cloneTips <- tipsWithNested[!tipsWithNested %in% c(
  tree$tip.label[Descendants(tree, 38)[[1]]][-1],
  tree$tip.label[Descendants(tree, 46)[[1]]][-1]
)]
cloneDriver <- "DNMT3A:p.R882H_AND_JAK2:p.V617F"
cloneName <- "clone2" # Match naming with PD6629_1 (first timepoint)
cloneTree <- keep.tip(tree, cloneTips)

# Plot to check (no nested)
ggtree(cloneTree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Add to the list of clones and move on
williams_clones <- append(williams_clones, list(cloneTree))
williamsCloneNames <- c(williamsCloneNames, paste0(ID, "_", cloneName, "_", cloneDriver, "_age", round(williams[[ID]]$age, 2)))



# PD6646_1: No nested in this one
ID <- "PD6646_1"
tree <- williams[[ID]]
williams[[ID]]$has_expanded_clades <- T

### Plot full tree again so we can manually identify clones (compare to Williams Fig 3)
ggtree(tree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Tips under node 95 are DNMT3A CBL clone
cloneTips <- c(tree$tip.label[Descendants(tree, 95)[[1]]])
cloneDriver <- "DNMT3A:p.?_AND_CBL:p.C401S"
cloneName <- "clone1"
cloneTree <- keep.tip(tree, cloneTips)

# Plot to check (no nested)
ggtree(cloneTree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Add to the list of clones and move on
williams_clones <- append(williams_clones, list(cloneTree))
williamsCloneNames <- c(williamsCloneNames, paste0(ID, "_", cloneName, "_", cloneDriver, "_age", round(williams[[ID]]$age, 2)))


### Plot full tree again so we can manually identify clones (compare to Williams Fig 3)
ggtree(tree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Tips under node 112 are DNMT3A JAK2 clone
cloneTips <- c(tree$tip.label[Descendants(tree, 112)[[1]]])
cloneDriver <- "DNMT3A:p.?_AND_JAK2:p.V617F"
cloneName <- "clone2"
cloneTree <- keep.tip(tree, cloneTips)

# Plot to check (no nested)
ggtree(cloneTree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Add to the list of clones and move on
williams_clones <- append(williams_clones, list(cloneTree))
williamsCloneNames <- c(williamsCloneNames, paste0(ID, "_", cloneName, "_", cloneDriver, "_age", round(williams[[ID]]$age, 2)))


### Plot full tree again so we can manually identify clones (compare to Williams Fig 3)
ggtree(tree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Tips under 163 are unknown expansion
cloneTips <- c(tree$tip.label[Descendants(tree, 163)[[1]]])
cloneDriver <- "unknown"
cloneName <- "clone3"
cloneTree <- keep.tip(tree, cloneTips)

# Plot to check (no nested)
ggtree(cloneTree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Add to the list of clones and move on
williams_clones <- append(williams_clones, list(cloneTree))
williamsCloneNames <- c(williamsCloneNames, paste0(ID, "_", cloneName, "_", cloneDriver, "_age", round(williams[[ID]]$age, 2)))


# PD6646_2: Remove 9pUPD
ID <- "PD6646_2"
tree <- williams[[ID]]
williams[[ID]]$has_expanded_clades <- T

### Plot full tree again so we can manually identify clones (compare to Williams Fig 3)
ggtree(tree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Under node 29 but must remove nested subclone
tipsWithNested <- c(tree$tip.label[Descendants(tree, 29)[[1]]])
cloneTips <- tipsWithNested[!tipsWithNested %in% c(tree$tip.label[Descendants(tree, 42)[[1]]][-1])]
cloneDriver <- "DNMT3A:p.?_AND_JAK2:p.V617F"
cloneName <- "clone2" # Match clone name from previous timepoint
cloneTree <- keep.tip(tree, cloneTips)

# Plot to check (no nested)
ggtree(cloneTree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Add to the list of clones and move on
williams_clones <- append(williams_clones, list(cloneTree))
williamsCloneNames <- c(williamsCloneNames, paste0(ID, "_", cloneName, "_", cloneDriver, "_age", round(williams[[ID]]$age, 2)))



# PD5117_1: Remove 9pUPD and TET2
ID <- "PD5117_1"
tree <- williams[[ID]]
williams[[ID]]$has_expanded_clades <- T

### Plot full tree again so we can manually identify clones (compare to Williams Fig 2)
ggtree(tree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Tips under node 93, but must remove the nested subclones
tipsWithNested <- c(tree$tip.label[Descendants(tree, 93)[[1]]])
cloneTips <- tipsWithNested[!tipsWithNested %in% c(
  tree$tip.label[Descendants(tree, 134)[[1]]][-1],
  tree$tip.label[Descendants(tree, 124)[[1]]][-1]
)]
cloneDriver <- "JAK2:p.V617F"
cloneName <- "clone1"
cloneTree <- keep.tip(tree, cloneTips)

# Plot to check (no nested)
ggtree(cloneTree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Add to the list of clones and move on
williams_clones <- append(williams_clones, list(cloneTree))
williamsCloneNames <- c(williamsCloneNames, paste0(ID, "_", cloneName, "_", cloneDriver, "_age", round(williams[[ID]]$age, 2)))



# PD4781_1: No nested
ID <- "PD4781_1"
tree <- williams[[ID]]
williams[[ID]]$has_expanded_clades <- T

### Plot full tree again so we can manually identify clones (compare to Williams Fig 3)
ggtree(tree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Tips under node 41 are a clone
cloneTips <- c(tree$tip.label[Descendants(tree, 41)[[1]]])
cloneDriver <- "JAK2:p.V617F_AND_9pUPD_AND_TET2:p.Q1632*_AND_7p-_AND_7q+"
cloneName <- "clone1"
cloneTree <- keep.tip(tree, cloneTips)

# Plot to check (no nested)
ggtree(cloneTree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Add to the list of clones and move on
williams_clones <- append(williams_clones, list(cloneTree))
williamsCloneNames <- c(williamsCloneNames, paste0(ID, "_", cloneName, "_", cloneDriver, "_age", round(williams[[ID]]$age, 2)))



# PD4781_2: No nested
ID <- "PD4781_2"
tree <- williams[[ID]]
williams[[ID]]$has_expanded_clades <- T

### Plot full tree again so we can manually identify clones (compare to Williams Fig 3)
ggtree(tree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

cloneTips <- c(tree$tip.label[Descendants(tree, 24)[[1]]])
cloneDriver <- "JAK2:p.V617F_AND_9pUPD_AND_TET2:p.Q1632*_AND_7p-_AND_7q+"
cloneName <- "clone1" # Match PD4781_1 naming
cloneTree <- keep.tip(tree, cloneTips)

# Plot to check (no nested)
ggtree(cloneTree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Add to the list of clones and move on
williams_clones <- append(williams_clones, list(cloneTree))
williamsCloneNames <- c(williamsCloneNames, paste0(ID, "_", cloneName, "_", cloneDriver, "_age", round(williams[[ID]]$age, 2)))



# PD5147_1: No nested
ID <- "PD5147_1"
tree <- williams[[ID]]
williams[[ID]]$has_expanded_clades <- T

### Plot full tree again so we can manually identify clones (compare to Williams Fig 4)
ggtree(tree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Tips under node 69 are clonal
cloneTips <- c(tree$tip.label[Descendants(tree, 69)[[1]]])
cloneDriver <- "PPM1D:p.T483fs*3_AND_TET2:p.S657fs*42"
cloneName <- "clone1"
cloneTree <- keep.tip(tree, cloneTips)

# Plot to check (no nested)
ggtree(cloneTree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()

# Add to the list of clones and move on
williams_clones <- append(williams_clones, list(cloneTree))
williamsCloneNames <- c(williamsCloneNames, paste0(ID, "_", cloneName, "_", cloneDriver, "_age", round(williams[[ID]]$age, 2)))



# No expansions > 10 tips (PD5163)
ID <- "PD5163_1"
tree <- williams[[ID]]
### Plot full tree again so we can manually identify clones (compare to Williams Fig 2)
ggtree(tree) + geom_text(aes(label = node), vjust = -.3) + layout_dendrogram()
williams[[ID]]$has_expanded_clades <- F


# Name the mitchell_clones list
names(williams_clones) <- williamsCloneNames











############## Van Egeren dataset downloaded from https://gitlab.com/hormozlab/stemcellsim
# Van Egeren data is only from the clone, so we don't subset
library(rtreefit) # Use this to make trees ultrametric, as williams does

# Clone age is 25 years for tree1 (clone only)
tree1 <- ape::read.tree("~/stemcellsim/StemCellSim/ET1_tree.txt")

# Convert to ultrametric, time-based tree
agedf <- data.frame("age" = 25, "tip.label" = tree1$tip.label)
tree1$agedf <- agedf
ultra_tree1 <- fit_tree(tree1, switch_nodes = c(), cores = 6, model = "nb_tree", niter = 2000)$ultratree

ultra_tree1$agedf <- NULL
ultra_tree1$sensitivity <- NULL

# Set patient age (not clone age or tree age)
ultra_tree1$age <- 34
ultra_tree1$diagnosis_age <- 34

ggtree(ultra_tree1) + layout_dendrogram()


# Clone age is 40 years for tree2
tree2 <- ape::read.tree("~/stemcellsim/StemCellSim/ET2_tree.txt")

# Convert to ultrametric, time-based tree
agedf <- data.frame("age" = 40, "tip.label" = tree2$tip.label)
tree2$agedf <- agedf
ultra_tree2 <- fit_tree(tree2, switch_nodes = c(), cores = 6, model = "nb_tree", niter = 2000)$ultratree

ggtree(ultra_tree2) + layout_dendrogram()

ultra_tree2$agedf <- NULL
ultra_tree2$sensitivity <- NULL

# Set patient age (not clone age or tree age)
ultra_tree2$age <- 63
ultra_tree2$diagnosis_age <- 63

# Note: we don't have full trees for Van Egeren, only clones
vanEgeren_clones <- list(ultra_tree1, ultra_tree2)
names(vanEgeren_clones) <- c(
  "vanEgerenET1_clone1_JAK2:p.V617F_age34",
  "vanEgerenET2_clone1_JAK2:p.V617F_age63"
)







# Make combined clone list
names(fabre_clones) <- paste0(names(fabre_clones), "_fabre")
names(mitchell_clones) <- paste0(names(mitchell_clones), "_mitchell")
names(williams_clones) <- paste0(names(williams_clones), "_williams")
names(vanEgeren_clones) <- paste0(names(vanEgeren_clones), "_vanEgeren")

# Merge all clone_lists into one list
clone_list <- c(fabre_clones, mitchell_clones, williams_clones, vanEgeren_clones)
fullTree_list <- c(fabre, mitchell, williams)

# Combine object containing full trees and clone trees as two separate lists
realCloneData <- list(fullTree_list, clone_list)
names(realCloneData) <- c("fullTrees", "cloneTrees")

# Remove "has_expanded_clades" because it doesn't make sense in the expanded clones list
for (i in c(1:length(realCloneData$cloneTrees))) {
  realCloneData$cloneTrees[[i]]$has_expanded_clades <- NULL
}

# Save as loadable data object for package users
usethis::use_data(realCloneData, overwrite = TRUE, compress = "bzip2")
