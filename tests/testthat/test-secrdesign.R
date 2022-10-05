## 2022-10-05 start

RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")

###############################################################################
test_that("random number generator stable", {
    set.seed(1235)
    expect_equal(rnorm(1), -0.6979879, tolerance = 1e-6)
})
###############################################################################

library(secrdesign)

## to avoid ASAN/UBSAN errors on CRAN, following advice of Kevin Ushey
## e.g. https://github.com/RcppCore/RcppParallel/issues/169
Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

###############################################################################

## GAminnr()

msk <- make.mask(type = 'rectangular', spacing = 10, nx = 30, ny = 20, buffer = 0)
alltrps <- make.grid(nx = 29, ny = 19, origin = c(10,10), spacing = 10)
set.seed(123)

# 10 generations for demonstration, use more in practice
opt <- GAminnr(msk, alltrps, ntraps = 20, detectpar = list(lambda0 = 0.5, sigma = 20), 
    detectfn = 'HHN', D = 10, noccasions = 5, ngen = 10, verbose = 0)

test_that("Genetic algorithm on track", {
  expect_equal(opt$des$bestobj, -46.01154, tolerance = 1e-4)
})
###############################################################################
