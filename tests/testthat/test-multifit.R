# fit = "multifit": several models fitted to the same simulated data
# 2026-11-05

test_that("multifit output is summarised by model", {
    skip_on_cran()
    trps <- make.grid(5, 5, spacing = 20, detector = "proximity")
    msk <- make.mask(trps, buffer = 60, nx = 16)
    scen <- make.scenarios(noccasions = 4, D = 8, g0 = 0.2, sigma = 20)
    fa <- list(list(
        null = list(model = list(D~1, g0~1, sigma~1)),
        bhv  = list(model = list(D~1, g0~b, sigma~1))))
    sims <- suppressMessages(
        run.scenarios(3, scen, trps, msk, fit = "multifit", fit.args = fa, seed = 1))
    expect_equal(sims$outputtype, "multipredicted")
    expect_equal(names(sims$output[[1]][[1]]), c("null", "bhv"))

    est <- estimateSummary(sims, "D")
    expect_equal(est$model, c("null", "bhv"))
    est1 <- estimateSummary(sims, "D", model = "bhv")
    expect_equal(est1$EST, est$EST[2])
    expect_equal(estimateSummary(sims, "D", model = 2)$EST, est1$EST)
    expect_error(estimateSummary(sims, "D", model = "nomodel"), "not found")
    expect_error(estimateArray(sims), "several models")
    expect_equal(nrow(countSummary(sims)), 1)

    ss <- select.stats(sims, "D", c("estimate", "RB", "RSE", "COV"))
    expect_equal(names(ss), c("null", "bhv"))
    expect_s3_class(ss$bhv, "selectedstatistics")
    ss1 <- select.stats(sims, "D", c("estimate", "RB"), model = "bhv")
    expect_s3_class(ss1, "selectedstatistics")
    expect_equal(mean(ss1$output[[1]][, "estimate"]), est1$EST)   # same replicates as estimateSummary

    # unnamed models with identical formulae are named model1, model2
    fa2 <- list(list(
        list(model = list(D~1, g0~1, sigma~1)),
        list(model = list(D~1, g0~1, sigma~1), details = list(distribution = "binomial"))))
    sims2 <- suppressMessages(
        run.scenarios(2, scen, trps, msk, fit = "multifit", fit.args = fa2, seed = 1))
    expect_equal(names(sims2$output[[1]][[1]]), c("model1", "model2"))
})
