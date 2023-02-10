## code to prepare `exampleTrees` dataset goes here

# Generate empty list and vector naming the list
tree_list <- list()

# Generate 100 trees, appending to list
for (i in c(1:100)) {
  a <- stats::runif(1, min = 1, max = 2)
  b <- a - 1
  cloneYrs <- 20
  n_tips <- 100
  tree <- cloneRate::simTree(a, b, cloneAge = cloneYrs, n = n_tips)
  nu <- stats::runif(1, 10, 20)

  for (i in c(1:length(tree$edge.length))) {
    tree$edge.length[i] <- stats::rpois(1, nu * tree$edge.length[i])
  }

  tree$params <- data.frame(r = 1, "a" = a, "b" = b, "cloneAge" = cloneYrs, "n" = n_tips, "nu" = nu)
  tree_list <- append(tree_list, list(tree))
}

# Rename object and save
exampleMutTrees <- tree_list
usethis::use_data(exampleMutTrees, overwrite = TRUE, compress = "bzip2")
