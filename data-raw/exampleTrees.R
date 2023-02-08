## code to prepare `exampleTrees` dataset goes here
library(coalRate)

# Generate empty list and vector naming the list
tree_list <- list()
names_vec <- c()

# Generate 100 trees, appending to list
for (i in c(1:100)){
  a <- stats::runif(1, min = 1, max = 2)
  b <- a-1
  tree <- coalRate::simTree(a, b, cloneAge = 20, n=100)
  tree_list <- append(tree_list, list(tree))
  names_vec <- c(names_vec, paste0("r=1_a=", round(a, 5), "_b=", round(b, 5),
                                   "_T=20_n=100")) # T is used synonymously with cloneAge
}

# Name the list so we have the info on params of the trees generated
names(tree_list) <- names_vec

# Rename object and save to data-raw (?) and data
exampleTrees <- tree_list
#save(exampleTrees, file = "data-raw/exampleTrees.rda", compress = "bzip2")
usethis::use_data(exampleTrees, overwrite = TRUE, compress = "bzip2")
