p0 <- c(0.10, 0.15, 0.15, 0.20, 0.20, 0.10, 0.10)

test_that("se_cont_po recovers beta", {
  fit_wh <- se_cont_po(p0, beta = log(1.5), sd_x = 0.5, n = 400, method = "whitehead")
  fit_ol <- se_cont_po(p0, beta = log(1.5), sd_x = 0.5, n = 400, method = "ologit")
  expect_equal(unname(fit_wh["beta"]), log(1.5))
  expect_equal(unname(fit_ol["beta"]), log(1.5), tolerance = 1e-2)
})

test_that("whitehead and ologit SEs agree closely for a small effect", {
  fit_wh <- se_cont_po(p0, beta = log(1.05), sd_x = 0.5, n = 400, method = "whitehead")
  fit_ol <- se_cont_po(p0, beta = log(1.05), sd_x = 0.5, n = 400, method = "ologit")
  expect_equal(unname(fit_wh["SE"]), unname(fit_ol["SE"]), tolerance = 0.02)
})

test_that("R2 > 0 increases the required n", {
  n0 <- n_precision_cont_po(p0, beta = log(1.5), sd_x = 0.5, ratio_UL = 3, R2 = 0, method = "ologit")
  n_adj <- n_precision_cont_po(p0, beta = log(1.5), sd_x = 0.5, ratio_UL = 3, R2 = 0.3, method = "ologit")
  expect_gt(n_adj$n, n0$n)
})

test_that("n_precision_cont_po's implied CI hits the precision target", {
  fit <- n_precision_cont_po(p0, beta = log(1.5), sd_x = 0.5, ratio_UL = 3, method = "whitehead")
  expect_equal(fit$uci / fit$lci, 3, tolerance = 1e-6)
})
