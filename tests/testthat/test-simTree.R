test_that("Number of tips works", {
  n_tips <- sample(seq(10, 100), 1)
  expect_equal(length(simTree(1, .5, 20, n_tips)$tip.label), n_tips)
})
