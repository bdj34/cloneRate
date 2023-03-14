test_that("shared mutations method needs mutation rate", {
  mutTree <- simMut(a = 1, b = 0, cloneAge = 20, nu = 10, n = 50)

  # sharedMuts() will find mutation rate if nu is a col in metadata data.frame
  expect_no_error(sharedMuts(mutTree))

  # Remove metadata and sharedMuts() should fail
  mutTree$metadata <- NULL
  expect_error(sharedMuts(mutTree), regexp = "Need to give a mutation rate")

  # Add nu to tree and sharedMuts() should work
  mutTree$nu <- 10
  expect_no_error(sharedMuts(mutTree))

  # Specify mutation rate and sharedMuts() should work without nu in tree or metadata
  mutTree$nu <- NULL
  expect_no_error(sharedMuts(mutTree, nu = 10))
})

test_that("Max. likelihood returns correct number and names", {
  check.df <- maxLikelihood(realCloneData$cloneTrees)
  expect_equal(check.df$cloneName_result, names(realCloneData$cloneTrees))
  expect_length(realCloneData$cloneTrees, nrow(check.df))
})

test_that("Internal lengths returns correct number and names", {
  check.df <- internalLengths(realCloneData$cloneTrees)
  expect_equal(check.df$cloneName_result, names(realCloneData$cloneTrees))
  expect_length(realCloneData$cloneTrees, nrow(check.df))
})

test_that("Moments returns correct number and names", {
  check.df <- moments(realCloneData$cloneTrees)
  expect_equal(check.df$cloneName_result, names(realCloneData$cloneTrees))
  expect_length(realCloneData$cloneTrees, nrow(check.df))
})

test_that("Shared muts. returns correct number and names", {
  check.df <- sharedMuts(exampleMutTrees)
  expect_equal(check.df$cloneName_result, names(exampleMutTrees)) # Should both be NULL
  expect_length(exampleMutTrees, nrow(check.df))

  # Now give exampleMutTrees names and check
  names(exampleMutTrees) <- paste0(c(1:length(exampleMutTrees)))
  check.df <- sharedMuts(exampleMutTrees)
  expect_equal(check.df$cloneName_result, names(exampleMutTrees))
  expect_length(exampleMutTrees, nrow(check.df))
})

test_that("ultrametric fns. return error when given mutation tree(s)", {
  mutTree <- simMut(a = 1, b = 0, cloneAge = 20, nu = 100, n = 10)
  expect_error(internalLengths(mutTree), regexp = "Tree is not ultrametric.")
  expect_error(maxLikelihood(mutTree), regexp = "Tree is not ultrametric.")
  expect_error(moments(mutTree), regexp = "Tree is not ultrametric.")

  mutTrees <- simMut(a = 1, b = 0, cloneAge = 20, nu = 100, n = 10, nTrees = 2)
  expect_error(internalLengths(mutTrees), regexp = "Tree is not ultrametric.")
  expect_error(maxLikelihood(mutTrees), regexp = "Tree is not ultrametric.")
  expect_error(moments(mutTrees), regexp = "Tree is not ultrametric.")
})

test_that("sharedMuts() throws error for ultrametric tree", {
  ultraTree <- simUltra(a = 1, b = 0, cloneAge = 10, n = 10)
  expect_error(sharedMuts(ultraTree, nu = 10), "Tree should be mutation-based")
})

test_that("sharedMuts throws error or warning for unreasonable alpha", {
  mutTree <- simMut(a = 1, b = 0, cloneAge = 20, nu = 100, n = 10)
  expect_error(sharedMuts(mutTree, alpha = 2), regexp = "alpha must be between 0 and 1")
  expect_warning(sharedMuts(mutTree, alpha = .95), regexp = "1-alpha confidence intervals")
})

test_that("Ultra fns. throw error or warning for unreasonable alpha", {
  ultraTree <- simUltra(a = 1, b = 0, cloneAge = 20, n = 10)
  expect_error(internalLengths(ultraTree, alpha = 2), regexp = "alpha must be between 0 and 1")
  expect_warning(maxLikelihood(ultraTree, alpha = .95), regexp = "1-alpha confidence intervals")
  expect_warning(moments(ultraTree, alpha = .95), regexp = "1-alpha confidence intervals")
})

test_that("all growth rate functions throw warning with non-binary trees", {
  nonBinaryTree <- cloneRate::realCloneData$cloneTrees$PD5847_1_clone1
  expect_warning(internalLengths(nonBinaryTree), regexp = "Tree is not binary.")
  expect_warning(maxLikelihood(nonBinaryTree), regexp = "Tree is not binary.")
  expect_warning(moments(nonBinaryTree), regexp = "Tree is not binary.")

  nonBinaryMutTree <- cloneRate::ultra2mut(nonBinaryTree, nu = 10)
  expect_warning(sharedMuts(nonBinaryMutTree), regexp = "Tree is not binary.")
})

test_that("MCMC gices expected output", {
  tree <- simUltra(a = 1, b = 0, cloneAge = 20, n = 20)
  mcmcOut <- birthDeathMCMC(tree)
  expect_true(inherits(mcmcOut, "data.frame"))
  expect_true(mcmcOut$estimate < 2)

  maxRate <- .5
  mcmcOut <- birthDeathMCMC(tree, maxGrowthRate = maxRate)
  expect_true(mcmcOut$estimate < maxRate)
})

test_that("MCMC gives expected error messages", {
  tree2 <- simUltra(a = 10, b = 0, cloneAge = 40, n = 20)
  expect_no_error(birthDeathMCMC(tree2))

  mutTree <- simMut(a = 1, b = 0, cloneAge = 20, n = 20, nu = 10)
  expect_error(birthDeathMCMC(mutTree), regexp = "not ultrametric")

  numTrees <- 2
  listTrees <- simUltra(a = 1, b = 0, cloneAge = 20, n = 20, nTrees = numTrees)
  expect_equal(nrow(birthDeathMCMC(listTrees)), numTrees)
})
