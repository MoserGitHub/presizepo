### Simulation confirmation
###
### Simulate ordinal data directly from the assumed proportional-odds model,
### fit it with MASS::polr(), and compare the empirical sampling SD of the
### log-OR (and the average model-based SE) against each formula's SE.
### Good agreement supports that the SE formulas above are correct.

library(MASS)

## Simulate ordinal outcomes from a cumulative-logit model with cutpoints
## `zeta` and linear predictor `eta` (cum. logit = zeta - eta, matching
## fit_cumlogit_weighted()'s convention).
sim_ordinal <- function(eta, zeta) {
  levels_n <- length(zeta) + 1
  u <- runif(length(eta))
  zeta_full <- c(-Inf, zeta, Inf)
  cum <- sapply(zeta_full, function(z) plogis(z - eta))
  lvl <- rowSums(u > cum[, -1, drop = FALSE]) + 1
  pmin(lvl, levels_n)
}

## --- Binary predictor: whitehead vs. ologit(SA) vs. ologit(SN) ---
sim_check_binary <- function(pc, OR, n, r = 1, n_boot = 1500, seed = 7351) {
  set.seed(seed)
  probs <- build_arm_probs(pc, OR = OR)
  zeta <- qlogis(cumsum(probs$p_ctrl)[-length(probs$p_ctrl)])
  b_true <- -log(OR)
  n_ctrl <- round(n / (1 + r)); n_trt <- n - n_ctrl
  x <- c(rep(0, n_ctrl), rep(1, n_trt))

  boot_b <- boot_se <- numeric(n_boot)
  for (i in seq_len(n_boot)) {
    y <- factor(sim_ordinal(b_true * x, zeta), levels = seq_len(length(zeta) + 1), ordered = TRUE)
    fit <- tryCatch(polr(y ~ x, Hess = TRUE), error = function(e) NULL)
    if (is.null(fit)) { boot_b[i] <- boot_se[i] <- NA; next }
    boot_b[i] <- coef(fit)[1]
    boot_se[i] <- sqrt(vcov(fit)[1, 1])
  }

  data.frame(
    quantity = c("empirical SD of b_hat", "mean model-based SE",
                 "se_binary(method = whitehead)", "se_binary(method = ologit)"),
    SE = c(sd(boot_b, na.rm = TRUE), mean(boot_se, na.rm = TRUE),
           se_binary_po(pc, OR = OR, n = n, r = r, method = "whitehead")$SE,
           se_binary_po(pc, OR = OR, n = n, r = r, method = "ologit")$SE)
  )
}

sim_check_binary(pc = p, OR = 1.7, n = 400)
sim_check_binary(pc = p, OR = 3, n = 400)  # larger effect: methods should diverge a bit more

## --- Continuous predictor: se_cont(whitehead) vs. se_cont(ologit) ---
sim_check_cont <- function(p0, beta, sd_x, n, n_boot = 1500, seed = 2864) {
  set.seed(seed)
  levels_n <- length(p0)
  zeta <- qlogis(cumsum(p0)[-levels_n])

  boot_beta <- boot_se <- numeric(n_boot)
  for (i in seq_len(n_boot)) {
    x <- rnorm(n, 0, sd_x)
    y <- factor(sim_ordinal(beta * x, zeta), levels = seq_len(levels_n), ordered = TRUE)
    fit <- tryCatch(polr(y ~ x, Hess = TRUE), error = function(e) NULL)
    if (is.null(fit)) { boot_beta[i] <- boot_se[i] <- NA; next }
    boot_beta[i] <- coef(fit)[1]
    boot_se[i] <- sqrt(vcov(fit)[1, 1])
  }

  data.frame(
    quantity = c("empirical SD of beta_hat", "mean model-based SE",
                 "se_cont(method = whitehead)", "se_cont(method = ologit)"),
    SE = c(sd(boot_beta, na.rm = TRUE), mean(boot_se, na.rm = TRUE),
           se_cont(p0, beta, sd_x, n, method = "whitehead")[["SE"]],
           se_cont(p0, beta, sd_x, n, method = "ologit")[["SE"]])
  )
}

sim_check_cont(p0, beta = log(1.5), sd_x = 0.5, n = 400)
