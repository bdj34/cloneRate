# Functions to adjust mutation trees and convert to ultrametric
# As they do in Mitchell et al. 2022 (https://www.nature.com/articles/s41586-022-04786-y)
# Pulled from github (https://github.com/emily-mitchell/normal_haematopoiesis)

# adjust edge len functions
get_corrected_tree = function(tree, details, sensitivity_df, include_indels = TRUE, include_SNVs = TRUE,get_edge_from_tree=FALSE) {
  tree_c = tree
  tree_c$edge.length = sapply(tree$edge[,2], correct_edge_length, tree = tree, details = details, sensitivity_df = sensitivity_df, include_indels = include_indels, include_SNVs=include_SNVs,get_edge_from_tree=get_edge_from_tree)
  return(tree_c)
}

#Note - only apply the "get_edge_from_tree" option if starting from an SNV only tree. Function will assume that all the existing edge length is SNVs.
correct_edge_length = function(node, tree, details, sensitivity_df, include_indels = TRUE, include_SNVs = TRUE, get_edge_from_tree=FALSE) {
  daughters <- getTips(tree = tree, node = node)
  #correct SNVs on edge, or set to 0 if want an indel only tree
  if(include_SNVs == TRUE) {
    if(get_edge_from_tree) {
      nSNV=tree$edge.length[tree$edge[,2]==node]
    } else {
      nSNV = sum(details$node == node & details$Mut_type == "SNV")
    }
    all_sens_SNVs <- sensitivity_df[sensitivity_df$Sample %in% daughters,"SNV_sensitivity"]
    branch_SNV_sens = 1 - prod(1-all_sens_SNVs)
    new_nSNV = nSNV/branch_SNV_sens
  } else {
    new_nSNV <- 0
  }
  #correct INDELs on edge, or set to 0 if want an SNV only tree
  if(include_indels == TRUE) {
    nINDEL = sum(details$node == node & details$Mut_type == "INDEL")
    all_sens_INDELs <- sensitivity_df[sensitivity_df$Sample %in% daughters,"INDEL_sensitivity"]
    branch_INDEL_sens = 1 - prod(1-all_sens_INDELs)
    new_nINDEL = nINDEL/branch_INDEL_sens
  } else {
    new_nINDEL <- 0
  }
  new_edge_length = new_nSNV + new_nINDEL
  return(new_edge_length)
}

getTips = function(tree,node) {
  require(ape)
  if(node <= length(tree$tip.label)) {
    daughters <- tree$tip.label[node]
  } else {
    daughters <- extract.clade(tree, node = node)$tip.label
  }
  return(daughters)
}

find.distance <- function(tree, from, to) {
  path <- nodepath(tree, from, to)
  res <- 0
  for (j in 2:length(path)) {
    index <- which(tree$edge[,2] == path[j] & tree$edge[,1] == path[j-1], arr.ind = TRUE)
    res <- res + tree$edge.length[index]
  }
  return(res)
}

length.normalise <- function(orig.tree, new.tree, curr.node, remaining.stick) {
  curr.node.children <- unlist(Descendants(orig.tree, curr.node, "children"))

  for (j in curr.node.children) {
    index <- which(orig.tree$edge[,2] == j & orig.tree$edge[,1] == curr.node, arr.ind = TRUE)

    if (j %in% orig.tree$tip.label) {
      new.tree$edge.length[index] <- remaining.stick
    } else {
      curr.node.tips <- unlist(Descendants(orig.tree, j, "tips"))
      curr.dist <- find.distance(orig.tree, curr.node, j)
      if (curr.dist == 0) {curr.dist <- 0.01} # So that no edge lengths are zero
      desc.lengths <- sapply(curr.node.tips, FUN = find.distance, tree = orig.tree, from = curr.node)
      new.tree$edge.length[index] <- remaining.stick * curr.dist / mean(desc.lengths)
      shorter.stick <- remaining.stick - new.tree$edge.length[index]

      # Call as recursive function
      new.tree <- length.normalise(orig.tree, new.tree, j, shorter.stick)
    }
  }
  return(new.tree)
}

make.ultrametric.tree <- function(tree) {
  root.number <- length(tree$tip.label) + 1
  ultra.tree <- length.normalise(tree, tree, root.number, 1)
  return(ultra.tree)
}
