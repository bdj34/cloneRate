## code to prepare `exampleTrees` dataset goes here

# Randomly sample 100 birth rates, keeping net growth rate equal to 1
a <- stats::runif(100, min = 1, max = 2)
b <- a - 1

# Set mutation rate for each of the 100 trees (uniform on 10 to 20 per yr.)
nu <- stats::runif(100, min = 10, max = 20)

# Generate 100 mutation trees, each with 100 tips and each with its own mutation rate
exampleMutTrees <- cloneRate::simMut(a = a, b = b, cloneAge = 20, n = 100, nu = nu, nTrees = length(a))

# Save
usethis::use_data(exampleMutTrees, overwrite = TRUE, compress = "bzip2")
