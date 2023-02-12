## code to prepare `exampleTrees` dataset goes here

# Randomly sample 100 birth rates, keeping net growth rate equal to 1
a <- stats::runif(100, min = 1, max = 2)
b <- a-1

# Generate 100 ultrametric trees
exampleUltraTrees <- cloneRate::simUltra(a = a, b = b, cloneAge = 20, n = 100, nTrees = length(a))

# Save
usethis::use_data(exampleUltraTrees, overwrite = TRUE, compress = "bzip2")
