test_that("Number of tips works", {
  n_tips <- sample(seq(10, 100), 1)
  expect_length(simUltra(a = 1, b = .5, cloneAge = 20, n = n_tips)$tip.label, n_tips)
  expect_length(simMut(a = 1, b = .5, cloneAge = 20, n = n_tips, nu = 15)$tip.label, n_tips)
})

test_that("Correct number of trees produced", {
  num_trees <- sample(c(2:10), 1)
  trees <- simUltra(a = 1, b = .1, cloneAge = 20, n = 50, nTrees = num_trees)
  expect_equal(length(trees), num_trees)
})

test_that("nTrees must be equal to length of params vectors (if > 1)", {
  num_trees <- 5
  expect_error(trees <- simUltra(a = 1, b = rep(.1, 3), cloneAge = 20, n = 50, nTrees = num_trees), regexp = "Input parameters must be length 1 or length equal to")
  expect_error(trees <- simUltra(a = rep(1, 2), b = .1, cloneAge = 20, n = 50, nTrees = num_trees), regexp = "Input parameters must be length 1 or length equal to")
  expect_error(trees <- simUltra(a = 1, b = .1, cloneAge = rep(20, 6), n = 50, nTrees = num_trees), regexp = "Input parameters must be length 1 or length equal to")
  expect_error(trees <- simUltra(a = 1, b = rep(.1, 3), cloneAge = rep(20, 5), n = rep(50, 10), nTrees = num_trees), regexp = "Input parameters must be length 1 or length equal to")
  expect_error(trees <- simUltra(a = rep(1, 2), b = .1, cloneAge = 5, n = 50, nTrees = 1), regexp = "Input parameters must be length 1 or length equal to")
})

test_that("birth rate can't be equal to or less than death rate", {
  expect_error(trees <- simUltra(a = 1, b = 1, cloneAge = 20, n = 50), regexp = "Birth rate")
  expect_error(trees <- simUltra(a = 1, b = 1.5, cloneAge = 20, n = 50), "Birth rate")
  expect_error(trees <- simUltra(a = rep(1, 5), b = c(rep(.5, 4), 1.5), cloneAge = 20, n = 50, nTrees = 5), "Birth rate")
})

test_that("warning if expected population size is less than n", {
  expect_warning(tree <- simUltra(a = 1, b = .9, cloneAge = 10, n = 10), regexp = "greater than expected population size")
  expect_warning(tree <- simUltra(a = 1, b = .5, cloneAge = 10, n = 200), regexp = "greater than expected population size")
  warning_vec <- capture_warnings(trees <- simUltra(a = 1, b = .5, cloneAge = 10, n = 200, nTrees = 3))
  expect_match(warning_vec, "greater than expected population size")
})

test_that("convert ultrametric to mutation based trees", {
  nYears <- 20
  ultraTrees <- simUltra(a = 1, b = 0, cloneAge = nYears, n = 20, nTrees = 5, addStem = T)
  mutRate <- 100
  mutTrees <- ultra2mut(ultraTrees, nu = mutRate)
  for (i in c(1:length(mutTrees))) {
    expect_equal(mutTrees[[i]]$metadata$nu, mutRate)
    expect_true(max(ape::branching.times(mutTrees[[i]])) < mutRate * 2 * nYears)
    expect_true(max(ape::branching.times(mutTrees[[i]])) > mutRate * .5 * nYears)
  }
})

test_that("utrametric tree has right cloneAge time if addStem = TRUE", {
  nYears <- 20
  ultraTrees <- simUltra(a = 1, b = 0, cloneAge = nYears, n = 20, nTrees = 5, addStem = T)
  for (i in c(1:length(ultraTrees))) {
    expect_true(max(suppressWarnings(ape::branching.times(ultraTrees[[i]]))) == nYears)
  }
})

test_that("n = 3 trees throw error", {
  expect_error(simUltra(a = 1, b = 0, cloneAge = 20, n = 2, nTrees = 1), regexp = "Number of samples")
  expect_error(simUltra(a = 1, b = 0, cloneAge = 20, n = 3, nTrees = 5), regexp = "Number of samples")
  expect_error(simMut(a = 1, b = 0, cloneAge = 20, n = 3, nTrees = 1, nu = 1), regexp = "Number of samples")
})

test_that("precBits too low throws error", {
  expect_error(simUltra(a = 1, b = 0, cloneAge = 20, n = 20, nTrees = 1, precBits = 20),
    regexp = "alpha value is equal to 1 due to insufficient machine precision."
  )
  expect_error(simUltra(a = 5, b = 0, cloneAge = 100, n = 20, nTrees = 1, precBits = 100),
    regexp = "alpha value is equal to 1 due to insufficient machine precision."
  )
  expect_error(simMut(a = 1, b = 0, cloneAge = 20, n = 20, nu = 10, nTrees = 1, precBits = 20),
    regexp = "alpha value is equal to 1 due to insufficient machine precision."
  )
})

test_that("coal_to_tree produces ultrametric phylo object", {
  expect_true(ape::is.ultrametric(coal_to_tree(runif(9, min = 5, max = 20))))
})
