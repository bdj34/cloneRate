test_that("shared mutations method needs mutation rate", {

  mutTree <- simMut(a = 1, b = 0, cloneAge = 20, nu = 10, n = 50)

  # sharedMuts() will find mutation rate if nu is a col in metadata data.frame
  expect_no_error(sharedMuts(mutTree))

  # Remove metadata and sharedMuts() should fail
  mutTree$metadata <- NULL
  expect_error(sharedMuts(mutTree), regexp = "Need to give a mutation rate")

  # Specify mutation rate and sharedMuts() should work
  expect_no_error(sharedMuts(mutTree, nu = 10))
})
