test_that("Kmeansimp works", {
  expect_equal(round(Kmeansimp(matrix(c(1:100), 50, 2), 2)
                     $`proportion of betweenSS with totalSS`, 2),  0.75)
  expect_equal(round(Kmeansimp(matrix(c(1:1000), 200, 5), 4, iter = 25)
                     $`proportion of betweenSS with totalSS`, 2),  0.94)
})
