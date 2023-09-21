test_that("tree gives normalized frequency of 1", {
  tree <- simUltra(a = 1, b = 0, cloneAge = 10, n = 40, addStem = TRUE)
  siteFreq.df <- siteFrequency(tree, includeStem = TRUE)
  expect_equal(sum(siteFreq.df$normalizedFreq), 1)

  tree <- simUltra(a = 1, b = 0, cloneAge = 10, n = 10)
  siteFreq.df <- siteFrequency(tree, includeStem = FALSE)
  expect_equal(sum(siteFreq.df$normalizedFreq), 1)
})

test_that("names preserved in input and output list", {
  tree_list <- simMut(a = 1, b = 0, cloneAge = 10, nu = 10, n = 10, nTrees = 3)
  names(tree_list) <- c("A_1", "B_2", "C_3")
  siteFreq.list <- siteFrequency(tree_list)
  expect_true(grepl(names(tree_list)[1], names(siteFreq.list)[1]))
  expect_true(grepl(names(tree_list)[2], names(siteFreq.list)[2]))
  expect_true(grepl(names(tree_list)[3], names(siteFreq.list)[3]))
})
