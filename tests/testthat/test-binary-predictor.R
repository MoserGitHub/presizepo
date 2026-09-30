pc <- c(0.15, 0.10, 0.10, 0.10, 0.10, 0.45)

test_that("build_arm_probs requires exactly one of OR/pe/rr", {
  expect_error(build_arm_probs(pc), "exactly one")
  expect_error(build_arm_probs(pc, OR = 1.5, rr = 1.2), "exactly one")
})

test_that("build_arm_probs returns probabilities summing to 1", {
  probs <- build_arm_probs(pc, OR = 1.7)
  expect_equal(sum(probs$p_trt), 1)
  expect_equal(probs$p_ctrl, pc)
})

test_that("build_arm_probs recovers OR = 1 as pc == p_trt", {
  probs <- build_arm_probs(pc, OR = 1)
  expect_equal(probs$p_trt, pc, tolerance = 1e-8)
})

test_that("se_binary_po recovers the input OR", {
  fit_wh <- se_binary_po(pc, OR = 1.7, n = 400, method = "whitehead")
  fit_ol <- se_binary_po(pc, OR = 1.7, n = 400, method = "ologit")
  expect_equal(fit_wh$OR, 1.7, tolerance = 1e-6)
  expect_equal(fit_ol$OR, 1.7, tolerance = 1e-3)
})

test_that("whitehead and ologit SEs agree closely for a small effect", {
  fit_wh <- se_binary_po(pc, OR = 1.1, n = 400, method = "whitehead")
  fit_ol <- se_binary_po(pc, OR = 1.1, n = 400, method = "ologit")
  expect_equal(fit_wh$SE, fit_ol$SE, tolerance = 0.02)
})

test_that("n_precision_binary_po gives a smaller n for a wider precision target", {
  n_narrow <- n_precision_binary_po(pc, OR = 2, ratio_UL = 2, method = "ologit")
  n_wide <- n_precision_binary_po(pc, OR = 2, ratio_UL = 4, method = "ologit")
  expect_gt(n_narrow$n, n_wide$n)
})

test_that("n_precision_binary_po's implied SE hits the precision target", {
  z <- qnorm(0.975)
  fit <- n_precision_binary_po(pc, OR = 2, ratio_UL = 3, method = "whitehead")
  expect_equal(fit$uci_or / fit$lci_or, 3, tolerance = 1e-6)
})
