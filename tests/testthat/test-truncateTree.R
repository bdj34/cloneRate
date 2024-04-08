test_that("getTipDescendants agrees with phangorn::Descendants", {

  # Test with full ultrametric tree
  realFullTree <- cloneRate::realCloneData$fullTrees[[sample(1:length(cloneRate::realCloneData$fullTrees), size = 1)]]
  node <- sample((ape::Ntip(realFullTree) + 1):ape::Nnode(realFullTree, internal.only = FALSE), size = 1)
  expect_setequal(getTipDescendants(realFullTree, node), phangorn::Descendants(realFullTree, node, type = "tips")[[1]])

  # Test with ultrametric clone tree
  realCloneTree <- cloneRate::realCloneData$cloneTrees[[sample(1:length(cloneRate::realCloneData$cloneTrees), size = 1)]]
  node <- sample((ape::Ntip(realCloneTree) + 1):ape::Nnode(realCloneTree, internal.only = FALSE), size = 1)
  expect_setequal(getTipDescendants(realCloneTree, node), phangorn::Descendants(realCloneTree, node, type = "tips")[[1]])

  # Test with the embryonic mutation-based tree
  embryonicMutTree <- cloneRate::embryonic_mutation_trees[[sample(1:length(cloneRate::embryonic_mutation_trees), size = 1)]]
  node <- sample((ape::Ntip(embryonicMutTree) + 1):ape::Nnode(embryonicMutTree, internal.only = FALSE), size = 1)
  expect_setequal(getTipDescendants(embryonicMutTree, node), phangorn::Descendants(embryonicMutTree, node, type = "tips")[[1]])

  # Test with the embryonic ultrametric tree
  embryonicUltraTree <- cloneRate::embryonic_time_trees[[sample(1:length(cloneRate::embryonic_time_trees), size = 1)]]
  node <- sample((ape::Ntip(embryonicUltraTree) + 1):ape::Nnode(embryonicUltraTree, internal.only = FALSE), size = 1)
  expect_setequal(getTipDescendants(embryonicUltraTree, node), phangorn::Descendants(embryonicUltraTree, node, type = "tips")[[1]])

  # Test with mutation tree from simulated data
  exampleMutTree <- cloneRate::exampleMutTrees[[sample(1:length(cloneRate::exampleMutTrees), size = 1)]]
  node <- sample((ape::Ntip(exampleMutTree) + 1):ape::Nnode(exampleMutTree, internal.only = FALSE), size = 1)
  expect_setequal(getTipDescendants(exampleMutTree, node), phangorn::Descendants(exampleMutTree, node, type = "tips")[[1]])

  # Test with ultrametric tree from simulated data
  exampleUltraTree <- cloneRate::exampleUltraTrees[[sample(1:length(cloneRate::exampleUltraTrees), size = 1)]]
  node <- sample((ape::Ntip(exampleUltraTree) + 1):ape::Nnode(exampleUltraTree, internal.only = FALSE), size = 1)
  expect_setequal(getTipDescendants(exampleUltraTree, node), phangorn::Descendants(exampleUltraTree, node, type = "tips")[[1]])

})


test_that("getImmediateParent agrees with phangorn::Ancestors", {

  # Test with ultrametric  full tree from real data
  realFullTree <- cloneRate::realCloneData$fullTrees[[sample(1:length(cloneRate::realCloneData$fullTrees), size = 1)]]
  node <- sample(1:ape::Nnode(realFullTree, internal.only = FALSE), size = 1)
  expect_setequal(getImmediateParent(realFullTree, node), phangorn::Ancestors(realFullTree, node, type = "parent")[[1]])

  # Test with ultrametric clone tree
  realCloneTree <- cloneRate::realCloneData$cloneTrees[[sample(1:length(cloneRate::realCloneData$cloneTrees), size = 1)]]
  node <- sample(1:ape::Nnode(realCloneTree, internal.only = FALSE), size = 1)
  expect_setequal(getImmediateParent(realCloneTree, node), phangorn::Ancestors(realCloneTree, node, type = "parent")[[1]])

  # Test with the embryonic mutation-based tree
  embryonicMutTree <- cloneRate::embryonic_mutation_trees[[sample(1:length(cloneRate::embryonic_mutation_trees), size = 1)]]
  node <- sample(1:ape::Nnode(embryonicMutTree, internal.only = FALSE), size = 1)
  expect_setequal(getImmediateParent(embryonicMutTree, node), phangorn::Ancestors(embryonicMutTree, node, type = "parent")[[1]])

  # Test with the embryonic ultrametric tree
  embryonicUltraTree <- cloneRate::embryonic_time_trees[[sample(1:length(cloneRate::embryonic_time_trees), size = 1)]]
  node <- sample(1:ape::Nnode(embryonicUltraTree, internal.only = FALSE), size = 1)
  expect_setequal(getImmediateParent(embryonicUltraTree, node), phangorn::Ancestors(embryonicUltraTree, node, type = "parent")[[1]])

  # Test with mutation tree from simulated data
  exampleMutTree <- cloneRate::exampleMutTrees[[sample(1:length(cloneRate::exampleMutTrees), size = 1)]]
  node <- sample(1:ape::Nnode(exampleMutTree, internal.only = FALSE), size = 1)
  expect_setequal(getImmediateParent(exampleMutTree, node), phangorn::Ancestors(exampleMutTree, node, type = "parent")[[1]])

  # Test with ultrametric tree from simulated data
  exampleUltraTree <- cloneRate::exampleUltraTrees[[sample(1:length(cloneRate::exampleUltraTrees), size = 1)]]
  node <- sample(1:ape::Nnode(exampleUltraTree, internal.only = FALSE), size = 1)
  expect_setequal(getImmediateParent(exampleUltraTree, node), phangorn::Ancestors(exampleUltraTree, node, type = "parent")[[1]])

})




test_that("Truncate tree always returns ultrametric tree (mutation or time-based)", {

  truncated <- truncate_tree(cloneRate::realCloneData$fullTrees, dist = runif(n=1, min = 0, max = 20))
  expect_true(all(ape::is.ultrametric.multiPhylo(truncated)))

  truncated <- truncate_tree(cloneRate::realCloneData$cloneTrees, dist = runif(n=1, min = 0, max = 20))
  expect_true(all(ape::is.ultrametric.multiPhylo(truncated)))

  truncated <- truncate_tree(cloneRate::exampleMutTrees, dist = 55)
  expect_true(all(ape::is.ultrametric.multiPhylo(truncated)))

  truncated <- truncate_tree(cloneRate::exampleUltraTrees, dist = 55)
  expect_true(all(ape::is.ultrametric.multiPhylo(truncated)))

  truncated <- truncate_tree(cloneRate::embryonic_mutation_trees, dist = 55)
  expect_true(all(ape::is.ultrametric.multiPhylo(truncated)))

  truncated <- truncate_tree(cloneRate::embryonic_time_trees, dist = .5)
  expect_true(all(ape::is.ultrametric.multiPhylo(truncated)))
})
