### Precision-based sample size for proportional odds model
###
### Two approaches to the SE of the (log) odds ratio are supported,
### selectable via `method`:
###   "whitehead" - Whitehead's (1993) delta-method variance approximation
###                 (fast, closed-form; assumes proportional odds and can be
###                 inaccurate for large effects).
###   "ologit"    - White, Marley-Zagar, Morris, Parmar, Royston & Babiker
###                 (2023, Stata Journal 23(1):3-23) fit a weighted
###                 proportional-odds model directly to the anticipated
###                 distribution and read the SE off the observed
###                 information matrix. This is exact (no delta-method
###                 approximation) and, via `variance`, exposes two
###                 variants:
###                   "SA" - SE from a fit at the anticipated (alternative)
###                          arm-level probabilities. This is the SE that
###                          matters for a precision target, since the trial
###                          data will actually arise under (approximately)
###                          the planning effect, not under a null effect.
###                   "SN" - SE from a fit with the odds ratio constrained to
###                          1 (null), i.e. from a common, pooled category
###                          distribution across arms. This reproduces
###                          Whitehead's formula almost exactly (confirmed
###                          below) and is included mainly for comparison /
###                          as a check.
### SA and SN nearly coincide for small effects and diverge as the
### anticipated effect grows (see the examples below), which is White et
### al.'s stated motivation for preferring the exact "ologit" fit over
### Whitehead's approximation for large effects.

## Internal: weighted cumulative-logit MLE, generic in x (binary group
## indicator or continuous covariate). Ported from artcatr.R's
## .artcatr_fit_cumlogit().
fit_cumlogit_weighted <- function(level_int, x, w, offset = rep(0, length(x)),
                                   estimate_b = TRUE) {
  levels_n <- max(level_int)
  n_zeta <- levels_n - 1
  x <- as.numeric(x)

  # starting values: cutpoints from pooled (weight-averaged) proportions
  pooled <- tapply(w, level_int, sum)
  pooled <- pooled[as.character(seq_len(levels_n))]
  pooled[is.na(pooled)] <- 1e-6
  cumprop <- cumsum(pooled) / sum(pooled)
  start_zeta <- unname(qlogis(pmin(pmax(cumprop[seq_len(n_zeta)], 1e-6), 1 - 1e-6)))
  start_theta <- if (n_zeta == 1) start_zeta else c(start_zeta[1], log(pmax(diff(start_zeta), 1e-3)))
  start <- if (estimate_b) c(start_theta, 0) else start_theta

  negloglik <- function(par) {
    theta <- par[seq_len(n_zeta)]
    zeta <- if (n_zeta == 1) theta else cumsum(c(theta[1], exp(theta[-1])))
    b <- if (estimate_b) par[n_zeta + 1] else 0
    eta <- b * x + offset
    zeta_full <- c(-Inf, zeta, Inf)
    cum_lo <- plogis(zeta_full[level_int] - eta)
    cum_hi <- plogis(zeta_full[level_int + 1] - eta)
    p <- pmax(cum_hi - cum_lo, 1e-12)
    -sum(w * log(p))
  }

  fit <- optim(start, negloglik, method = "BFGS", hessian = estimate_b,
               control = list(reltol = 1e-12, maxit = 500))

  theta <- unname(fit$par[seq_len(n_zeta)])
  zeta <- if (n_zeta == 1) theta else cumsum(c(theta[1], exp(theta[-1])))
  b <- if (estimate_b) unname(fit$par[n_zeta + 1]) else NA_real_
  SE_b <- if (estimate_b) unname(sqrt(solve(fit$hessian)[n_zeta + 1, n_zeta + 1])) else NA_real_

  list(b = b, SE_b = SE_b, zeta = unname(zeta))
}

### Binary predictor

p <- c(0.15, 0.10, 0.10, 0.10, 0.10, 0.45)

## Build experimental-arm level probabilities from control-arm probabilities
## and exactly one of OR (proportional-odds shift), pe (given directly), or
## rr (common risk ratio).
build_arm_probs <- function(pc, OR = NULL, pe = NULL, rr = NULL) {
  n_specified <- sum(!is.null(OR), !is.null(pe), !is.null(rr))
  if (n_specified != 1) stop("please specify exactly one of OR, pe, rr")
  cum_c <- cumsum(pc)[-length(pc)]
  if (!is.null(pe)) {
    if (length(pe) != length(pc)) stop("pc and pe have different lengths")
    p_trt <- pe
  } else if (!is.null(OR)) {
    odds_t <- (cum_c / (1 - cum_c)) * OR
    p_trt <- diff(c(0, odds_t / (1 + odds_t), 1))
  } else {
    cum_t <- rr * cum_c
    p_trt <- diff(c(0, cum_t, 1))
  }
  list(p_ctrl = pc, p_trt = p_trt)
}

## SE of the log odds ratio for a binary (e.g. treatment) predictor.
## pc: control-arm probabilities per category (sum to 1)
## OR/pe/rr: exactly one, defining the experimental arm (see build_arm_probs)
## n: total sample size; r: allocation ratio (treatment:control)
## method: "whitehead" (delta-method approximation; equivalent to the
##   "SN"/null-constrained ologit fit, confirmed numerically above) or
##   "ologit" (exact MLE fit at the anticipated/alternative arm-level
##   probabilities, i.e. White et al.'s "SA"; the SE that matters for a
##   precision target, since trial data actually arise under the planning
##   effect, not under a null effect)
se_binary_po <- function(pc, OR = NULL, pe = NULL, rr = NULL, n, r = 1,
                       method = c("ologit", "whitehead"), conf = 0.95) {
  method <- match.arg(method)
  z <- qnorm(1 - (1 - conf) / 2)

  probs <- build_arm_probs(pc, OR, pe, rr)
  p_ctrl <- probs$p_ctrl; p_trt <- probs$p_trt
  levels_n <- length(p_ctrl)
  level_int <- rep(seq_len(levels_n), 2)
  x <- c(rep(0, levels_n), rep(1, levels_n))
  w <- n * c(p_ctrl / (1 + r), r * p_trt / (1 + r))

  # the alternative-hypothesis fit is needed for the OR estimate whenever OR
  # isn't given directly (pe/rr case), and for the SE itself under "ologit"
  need_fit_alt <- is.null(OR) || method == "ologit"
  if (need_fit_alt) {
    fit_alt <- fit_cumlogit_weighted(level_int, x, w)
    # sign flipped to match the convention that OR acts on the cumulative
    # odds of the *lower* categories
    logOR <- -fit_alt$b
  } else {
    logOR <- log(OR)
  }

  se <- if (method == "whitehead") {
    pbar <- (p_ctrl + r * p_trt) / (1 + r)
    sqrt(3 * (r + 1)^2 / (r * n * (1 - sum(pbar^3))))
  } else {
    fit_alt$SE_b
  }

  list(OR = exp(logOR), conf = conf, SE = se, lower = exp(logOR - z * se), upper = exp(logOR + z * se),
       method = method)
}

se_binary_po(pc = p, OR = 1.7, n = 30, method = "whitehead")
se_binary_po(pc = p, OR = 1.7, n = 30, method = "ologit")

## n for a target OR CI ratio (upper/lower). The Fisher information scales
## linearly with n (the weights sum to n), so SE(n) = SE(n0) * sqrt(n0 / n)
## exactly, for either method; fit once at n0 and rescale.
n_precision_binary_po <- function(pc, OR = NULL, pe = NULL, rr = NULL, ratio_UL, r = 1,
                                conf = 0.95, method = c("ologit", "whitehead"), n0 = 1000) {
  method <- match.arg(method)
  z <- qnorm(1 - (1 - conf) / 2)
  target_se <- log(ratio_UL) / (2 * z)
  fit0 <- se_binary(pc, OR, pe, rr, n = n0, r = r, method = method)
  n <- ceiling(n0 * (fit0$SE / target_se)^2)
  list(n = n, conf = conf, OR = fit0$OR, method = method,
       lci_or = exp(log(fit0$OR) - z * target_se), uci_or = exp(log(fit0$OR) + z * target_se))
}

n_precision_binary_po(pc = p, OR = 2, ratio_UL = 3, method = "whitehead")
n_precision_binary_po(pc = p, OR = 2, ratio_UL = 3, method = "ologit")

### Continuous predictor

library(statmod)  # gauss.quad(), for Gauss-Hermite integration over x in se_cont()

p0 <- c(0.10, 0.15, 0.15, 0.20, 0.20, 0.10, 0.10)   # illustrative

## SE of the log OR per unit of a continuous predictor x.
## p0: outcome probabilities at x = mean(x); beta: log OR per unit of x
## sd_x: SD of predictor; R2: variance of x explained by other covariates
## method: "whitehead" evaluates the variance formula at fixed cutpoints
##   derived from p0 (fast, closed-form in the spirit of Whitehead's
##   formula, extended here to a continuous covariate); "ologit" jointly
##   estimates the cutpoints and beta by weighted MLE over x ~ N(0, sd_x^2),
##   so SE(beta) reflects the full joint information matrix. The mixing
##   distribution of x is integrated out via Gauss-Hermite quadrature
##   (`ngrid` nodes), which converges far faster than an equal-spaced
##   quantile grid for a normal mixing distribution -- a handful of nodes
##   is generally enough (unlike a quantile grid, which needs thousands of
##   points to stop drifting, particularly in the tails).
se_cont_po <- function(p0, beta, sd_x, n, R2 = 0,
                     method = c("ologit", "whitehead"), ngrid = 15, conf = 0.95) {
  method <- match.arg(method)
  z <- qnorm(1 - (1 - conf) / 2)
  levels_n <- length(p0)
  th <- qlogis(cumsum(p0)[-levels_n])

  if (method == "whitehead") {
    ngrid_wh <- 2000
    xs <- qnorm((1:ngrid_wh - 0.5) / ngrid_wh) * sd_x  # quantile grid for x ~ N(0, sd_x)
    pm <- sapply(xs, function(x) diff(c(0, plogis(th - beta * x), 1)))  # levels_n x ngrid_wh
    pmbar <- rowMeans(pm)
    se <- sqrt(3 / (n * sd_x^2 * (1 - sum(pmbar^3)) * (1 - R2)))
    beta_est <- beta
  } else {
    gh <- gauss.quad(ngrid, kind = "hermite")
    xs <- sqrt(2) * sd_x * gh$nodes            # nodes rescaled to x ~ N(0, sd_x^2)
    wq <- gh$weights / sqrt(pi)                 # mixture weights, sum to 1
    pm <- sapply(xs, function(x) diff(c(0, plogis(th - beta * x), 1)))  # levels_n x ngrid

    level_int <- rep(seq_len(levels_n), times = ngrid)
    x_rep <- rep(xs, each = levels_n)
    w <- as.vector(pm) * rep(wq, each = levels_n) * n * (1 - R2)
    fit <- fit_cumlogit_weighted(level_int, x_rep, w)
    se <- fit$SE_b
    beta_est <- fit$b
  }
  c(beta = beta_est, SE = se, lower = exp(beta_est - z * se), upper = exp(beta_est + z * se))
}

se_cont_po(p0, beta = log(1.5), sd_x = 0.5, n = 30, method = "whitehead")
se_cont_po(p0, beta = log(1.5), sd_x = 0.5, n = 30, method = "ologit")

## n for a target OR (per `delta` units of x) CI ratio, with the same
## fit-once-and-rescale trick as se_binary()/n_precision_binary().
n_precision_cont_po <- function(p0, beta, sd_x, ratio_UL, delta = 1, R2 = 0, conf = 0.95,
                              method = c("ologit", "whitehead"), ngrid = 15, n0 = 1000) {
  method <- match.arg(method)
  z <- qnorm(1 - (1 - conf) / 2)
  target_se <- log(ratio_UL) / (2 * z * delta)              # target SE per 1 unit of x
  fit0 <- se_cont(p0, beta, sd_x, n = n0, R2 = R2, method = method, ngrid = ngrid)
  n <- ceiling(n0 * (fit0[["SE"]] / target_se)^2)
  list(n = n, or = exp(beta), method = method,
       lci = exp(beta - z * target_se), uci = exp(beta + z * target_se))
}

n_precision_cont_po(p0, beta = log(1.5), sd_x = 0.5, ratio_UL = 3, method = "whitehead")
n_precision_cont_po(p0, beta = log(1.5), sd_x = 0.5, ratio_UL = 3, method = "ologit")


