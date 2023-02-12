## code to prepare `exampleTrees` dataset goes here

# Randomly sample 100 birth rates, keeping net growth rate equal to 1
a <- stats::runif(100, min = 1, max = 2)
b <- a-1

# Generate 100 mutation trees, each with its own mutation rate (uniform on nu_min to nu_max)
exampleMutTrees <- cloneRate::simMut(a=a, b=b, cloneAge=20, n = 100, nu_min = 10, nu_max =20, nTrees = length(a))

# Save
usethis::use_data(exampleMutTrees, overwrite = TRUE, compress = "bzip2")
