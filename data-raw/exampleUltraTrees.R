## code to prepare `exampleTrees` dataset goes here
library(cloneRate)

# Generate empty list and vector naming the list
tree_list <- list()
names_vec <- c()

# Generate 100 trees, appending to list
for (i in c(1:100)){
  a <- stats::runif(1, min = 1, max = 2)
  b <- a-1
  cloneYrs <- 20
  n_tips = 100
  tree <- cloneRate::simTree(a, b, cloneAge = cloneYrs, n=n_tips)
  tree$params <- data.frame(r=1, "a" = a, "b" = b, "cloneAge" = cloneYrs, "n" = n_tips)
  tree_list <- append(tree_list, list(tree))

}

# Name the list so we have the info on params of the trees generated
names(tree_list) <- names_vec

# Rename object and save to data-raw (?) and data
exampleUltraTrees <- tree_list
#save(exampleTrees, file = "data-raw/exampleTrees.rda", compress = "bzip2")
usethis::use_data(exampleUltraTrees, overwrite = TRUE, compress = "bzip2")
