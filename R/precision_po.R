#' Build experimental-arm category probabilities
#'
#' Given control-arm outcome-category probabilities, builds the
#' experimental-arm probabilities implied by exactly one of a proportional
#' (cumulative) odds ratio `OR`, an experimental-arm probability vector `pe`
#' given directly, or a common risk ratio `rr` applied to the cumulative
#' probabilities.
#'
#' @param pc Numeric vector of control-arm category probabilities, summing
#'   to 1, ordered from the lowest to the highest category.
#' @param OR Proportional (cumulative) odds ratio for the experimental vs.
#'   control arm. Exactly one of `OR`, `pe`, `rr` must be supplied.
#' @param pe Experimental-arm category probabilities, given directly (same
#'   length as `pc`).
#' @param rr A common risk ratio applied to the control-arm cumulative
#'   probabilities.
#'
#' @return A list with elements `p_ctrl` and `p_trt`, the control- and
#'   experimental-arm category probability vectors.
#' @noRd
#' @examples
#' build_arm_probs(c(0.15, 0.10, 0.10, 0.10, 0.10, 0.45), OR = 1.7)
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

#' Standard error of the log odds ratio, binary predictor
#'
#' Computes the standard error (and confidence interval) of the log odds
#' ratio for a binary (e.g. treatment) predictor in a proportional-odds
#' model, at a given total sample size `n`.
#'
#' Two approaches to the SE are supported, selectable via `method`:
#' * `"whitehead"` - Whitehead's (1993) delta-method variance approximation
#'   (fast, closed-form; assumes proportional odds and can be inaccurate for
#'   large effects).
#' * `"ologit"` - White, Marley-Zagar, Morris, Parmar, Royston & Babiker
#'   (2023) fit a weighted proportional-odds model directly to the
#'   anticipated distribution and read the SE off the observed information
#'   matrix. This is exact (no delta-method approximation): the SE is
#'   computed from a fit at the anticipated (alternative) arm-level
#'   probabilities, i.e. their "SA" variant, since trial data actually arise
#'   under (approximately) the planning effect, not under a null effect.
#'
#' `"whitehead"` and `"ologit"` nearly coincide for small effects and diverge
#' as the anticipated effect grows, which is White et al.'s stated
#' motivation for preferring the exact `"ologit"` fit for large effects.
#'
#' @param pc Control-arm category probabilities (sum to 1).
#' @param OR,pe,rr Exactly one, defining the experimental arm; see
#'   [build_arm_probs()].
#' @param n Total sample size (both arms).
#' @param r Allocation ratio, treatment:control.
#' @param method `"ologit"` (default; exact MLE) or `"whitehead"`
#'   (delta-method approximation).
#' @param conf Confidence level for the reported interval.
#'
#' @return A list with elements `OR`, `conf`, `SE` (of the log OR), `lower`,
#'   `upper` (CI limits for the OR), and `method`.
#' @export
#'
#' @examples
#' pc <- c(0.15, 0.10, 0.10, 0.10, 0.10, 0.45)
#' se_binary_po(pc, OR = 1.7, n = 30, method = "whitehead")
#' se_binary_po(pc, OR = 1.7, n = 30, method = "ologit")
se_binary_po <- function(p0, OR = NULL, pe = NULL, rr = NULL, n, r = 1,
                         method = c("ologit", "whitehead"), conf = 0.95) {
  method <- match.arg(method)
  z <- stats::qnorm(1 - (1 - conf) / 2)

  probs <- build_arm_probs(p0, OR, pe, rr)
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

#' Sample size for a target OR precision, binary predictor
#'
#' Computes the total sample size `n` needed so that the confidence
#' interval for the odds ratio has a target ratio of upper to lower limit,
#' for a binary predictor in a proportional-odds model. Since the Fisher
#' information scales linearly with `n`, the SE at any `n` is
#' `SE(n0) * sqrt(n0 / n)` exactly (for either method); the fit is done once
#' at `n0` and rescaled analytically.
#'
#' @inheritParams se_binary_po
#' @param ratio_UL Target ratio of the upper to lower confidence limit for
#'   the OR.
#' @param n0 Sample size at which the one-time reference fit is done before
#'   rescaling; the default is large enough that further increases do not
#'   change the result.
#'
#' @return A list with elements `n`, `conf`, `OR`, `method`, `lci_or`,
#'   `uci_or`.
#' @export
#'
#' @examples
#' pc <- c(0.15, 0.10, 0.10, 0.10, 0.10, 0.45)
#' n_precision_binary_po(pc, OR = 2, ratio_UL = 3, method = "whitehead")
#' n_precision_binary_po(pc, OR = 2, ratio_UL = 3, method = "ologit")
n_precision_binary_po <- function(pc, OR = NULL, pe = NULL, rr = NULL, ratio_UL, r = 1,
                                  conf = 0.95, method = c("ologit", "whitehead"), n0 = 1000) {
  method <- match.arg(method)
  z <- stats::qnorm(1 - (1 - conf) / 2)
  target_se <- log(ratio_UL) / (2 * z)
  fit0 <- se_binary_po(pc, OR, pe, rr, n = n0, r = r, method = method)
  n <- ceiling(n0 * (fit0$SE / target_se)^2)
  list(n = n, conf = conf, OR = fit0$OR, method = method,
       lci_or = exp(log(fit0$OR) - z * target_se), uci_or = exp(log(fit0$OR) + z * target_se))
}


#' Resolve the per-unit log OR for a continuous predictor (internal)
#'
#' `beta` (per raw unit of `x`) and `OR_sd` (odds ratio per SD of `x`) are
#' two ways of specifying the same effect; exactly one must be supplied.
#' This converts either into the per-unit log OR used internally by
#' [se_cont_po()] and [n_precision_cont_po()].
#'
#' @param beta Log OR per unit of `x`, or `NULL`.
#' @param OR_sd Odds ratio per SD of `x`, or `NULL`.
#' @param sd_x SD of the predictor `x`.
#'
#' @return The per-unit log OR (a single number).
#' @noRd
resolve_beta_cont <- function(beta, OR_sd, sd_x) {
  n_specified <- sum(!is.null(beta), !is.null(OR_sd))
  if (n_specified != 1) stop("please specify exactly one of beta, OR_sd")
  if (is.null(beta)) beta <- log(OR_sd) / sd_x
  beta
}

#' Standard error of the log OR per unit of a continuous predictor
#'
#' Computes the standard error (and confidence interval) of the log odds
#' ratio per unit of a continuous predictor `x` in a proportional-odds
#' model, at a given total sample size `n`.
#'
#' @param p0 Outcome category probabilities at `x = mean(x)`.
#' @param beta Log OR per unit of `x`. Exactly one of `beta`, `OR_sd` must
#'   be supplied.
#' @param OR_sd Odds ratio per SD of `x` (i.e. the OR comparing `x` one SD
#'   above the mean to the mean), given as an alternative to `beta` so
#'   `sd_x` and the effect size don't need to be reconciled by hand;
#'   internally converted to `beta = log(OR_sd) / sd_x`. Exactly one of
#'   `beta`, `OR_sd` must be supplied.
#' @param sd_x SD of the predictor `x`.
#' @param n Total sample size.
#' @param R2 Proportion of variance of `x` explained by other covariates
#'   (their inclusion reduces the effective information for `beta` by a
#'   factor of `1 - R2`).
#' @param method `"ologit"` (default) jointly estimates the cutpoints and
#'   `beta` by weighted MLE over `x ~ N(0, sd_x^2)`, so `SE(beta)` reflects
#'   the full joint information matrix; the mixing distribution of `x` is
#'   integrated out via Gauss-Hermite quadrature (`ngrid` nodes), which
#'   converges quickly for a normal mixing distribution. `"whitehead"`
#'   evaluates a closed-form variance formula (in the spirit of Whitehead's
#'   1993 formula, extended here to a continuous covariate) at fixed
#'   cutpoints derived from `p0`.
#' @param ngrid Number of Gauss-Hermite quadrature nodes (only used for
#'   `method = "ologit"`).
#' @param conf Confidence level for the reported interval.
#' @param delta Number of units of `x` over which the OR (and its CI) is
#'   expressed; `beta`, `SE`, and the CI limits are all scaled by `delta`
#'   before exponentiating (the linear predictor is `beta * x`, so this
#'   rescaling is exact). Defaults to 1 (per 1 unit of `x`).
#'
#' @return A named numeric vector with elements `beta`, `SE`, `lower`,
#'   `upper` (the last two on the OR scale), all expressed for a change of
#'   `delta` units of `x`.
#' @export
#'
#' @examples
#' p0 <- c(0.10, 0.15, 0.15, 0.20, 0.20, 0.10, 0.10)
#' se_cont_po(p0, beta = log(1.5), sd_x = 0.5, n = 30, method = "whitehead")
#' se_cont_po(p0, beta = log(1.5), sd_x = 0.5, n = 30, method = "ologit")
#' # equivalently, specifying the effect as an OR per SD of x
#' se_cont_po(p0, OR_sd = 1.5^0.5, sd_x = 0.5, n = 30, method = "whitehead")
#' # OR/CI for a 2-unit (here, 2-SD, since sd_x = 1) change in x
#' se_cont_po(p0, beta = log(1.5), sd_x = 1, n = 30, delta = 2, method = "whitehead")
se_cont_po <- function(p0, beta = NULL, sd_x, n, R2 = 0, OR_sd = NULL,
                       method = c("ologit", "whitehead"), ngrid = 15, conf = 0.95,
                       delta = 1) {
  method <- match.arg(method)
  beta <- resolve_beta_cont(beta, OR_sd, sd_x)
  z <- stats::qnorm(1 - (1 - conf) / 2)
  levels_n <- length(p0)
  th <- stats::qlogis(cumsum(p0)[-levels_n])

  if (method == "whitehead") {
    ngrid_wh <- 2000
    xs <- stats::qnorm((1:ngrid_wh - 0.5) / ngrid_wh) * sd_x  # quantile grid for x ~ N(0, sd_x)
    pm <- sapply(xs, function(x) diff(c(0, stats::plogis(th - beta * x), 1)))  # levels_n x ngrid_wh
    pmbar <- rowMeans(pm)
    se <- sqrt(3 / (n * sd_x^2 * (1 - sum(pmbar^3)) * (1 - R2)))
    beta_est <- beta
  } else {
    gh <- statmod::gauss.quad(ngrid, kind = "hermite")
    xs <- sqrt(2) * sd_x * gh$nodes            # nodes rescaled to x ~ N(0, sd_x^2)
    wq <- gh$weights / sqrt(pi)                 # mixture weights, sum to 1
    pm <- sapply(xs, function(x) diff(c(0, stats::plogis(th - beta * x), 1)))  # levels_n x ngrid

    level_int <- rep(seq_len(levels_n), times = ngrid)
    x_rep <- rep(xs, each = levels_n)
    w <- as.vector(pm) * rep(wq, each = levels_n) * n * (1 - R2)
    fit <- fit_cumlogit_weighted(level_int, x_rep, w)
    se <- fit$SE_b
    beta_est <- fit$b
  }

  beta_d <- delta * beta_est
  se_d <- delta * se
  c(beta = beta_d, SE = se_d, lower = exp(beta_d - z * se_d), upper = exp(beta_d + z * se_d))
}

#' Sample size for a target OR precision, continuous predictor
#'
#' Computes the total sample size `n` needed so that the confidence
#' interval for the OR per `delta` units of a continuous predictor `x` has
#' a target ratio of upper to lower limit, using the same
#' fit-once-and-rescale trick as [n_precision_binary_po()].
#'
#' @inheritParams se_cont_po
#' @param ratio_UL Target ratio of the upper to lower confidence limit for
#'   the OR per `delta` units of `x`.
#' @param delta Number of units of `x` over which the OR (and its CI) is
#'   expressed.
#' @param n0 Sample size at which the one-time reference fit is done before
#'   rescaling; the default is large enough that further increases do not
#'   change the result.
#'
#' @return A list with elements `n`, `or`, `method`, `lci`, `uci`. `or`,
#'   `lci`, and `uci` are expressed for a change of `delta` units of `x`
#'   (not a 1-unit change), matching the scale on which `ratio_UL` was
#'   specified.
#' @export
#'
#' @examples
#' p0 <- c(0.10, 0.15, 0.15, 0.20, 0.20, 0.10, 0.10)
#' n_precision_cont_po(p0, beta = log(1.5), sd_x = 0.5, ratio_UL = 3, method = "whitehead")
#' n_precision_cont_po(p0, beta = log(1.5), sd_x = 0.5, ratio_UL = 3, method = "ologit")
#' # equivalently, specifying the effect as an OR per SD of x
#' n_precision_cont_po(p0, OR_sd = 1.5^0.5, sd_x = 0.5, ratio_UL = 3, method = "whitehead")
n_precision_cont_po <- function(p0, beta = NULL, sd_x, ratio_UL, delta = 1, R2 = 0, conf = 0.95,
                                OR_sd = NULL, method = c("ologit", "whitehead"), ngrid = 15,
                                n0 = 1000) {
  method <- match.arg(method)
  beta <- resolve_beta_cont(beta, OR_sd, sd_x)
  z <- stats::qnorm(1 - (1 - conf) / 2)
  target_se <- log(ratio_UL) / (2 * z * delta)              # target SE per 1 unit of x
  fit0 <- se_cont_po(p0, beta = beta, sd_x = sd_x, n = n0, R2 = R2, method = method, ngrid = ngrid)
  n <- ceiling(n0 * (fit0[["SE"]] / target_se)^2)
  list(n = n, or = exp(delta * beta), method = method,
       lci = exp(delta * (beta - z * target_se)), uci = exp(delta * (beta + z * target_se)))
}


#' Weighted cumulative-logit MLE (internal)
#'
#' Fits a proportional-odds (cumulative logit) model to a weighted table of
#' outcome-category counts, generic in `x` (a binary group indicator or a
#' continuous covariate). Used internally as the "exact" (`method = "ologit"`)
#' engine behind [se_binary_po()] and [se_cont_po()]. Ported from an internal
#' `artcatr.R` helper.
#'
#' @param level_int Integer vector giving the outcome category (1-based) for
#'   each row of the weighted "data".
#' @param x Numeric vector, same length as `level_int`: the predictor value
#'   for each row (0/1 for a binary group, or a continuous covariate value).
#' @param w Numeric vector of (non-negative) weights, same length as
#'   `level_int`.
#' @param offset Numeric vector of offsets added to the linear predictor,
#'   same length as `level_int`. Defaults to no offset.
#' @param estimate_b If `FALSE`, fixes the slope `b` at zero and estimates
#'   only the cutpoints (used for the null-constrained fit).
#'
#' @return A list with elements `b` (slope estimate), `SE_b` (its standard
#'   error, `NA` if `estimate_b = FALSE`), and `zeta` (fitted cutpoints).
#' @noRd
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
  start_zeta <- unname(stats::qlogis(pmin(pmax(cumprop[seq_len(n_zeta)], 1e-6), 1 - 1e-6)))
  start_theta <- if (n_zeta == 1) start_zeta else c(start_zeta[1], log(pmax(diff(start_zeta), 1e-3)))
  start <- if (estimate_b) c(start_theta, 0) else start_theta

  negloglik <- function(par) {
    theta <- par[seq_len(n_zeta)]
    zeta <- if (n_zeta == 1) theta else cumsum(c(theta[1], exp(theta[-1])))
    b <- if (estimate_b) par[n_zeta + 1] else 0
    eta <- b * x + offset
    zeta_full <- c(-Inf, zeta, Inf)
    cum_lo <- stats::plogis(zeta_full[level_int] - eta)
    cum_hi <- stats::plogis(zeta_full[level_int + 1] - eta)
    p <- pmax(cum_hi - cum_lo, 1e-12)
    -sum(w * log(p))
  }

  fit <- stats::optim(start, negloglik, method = "BFGS", hessian = estimate_b,
                      control = list(reltol = 1e-12, maxit = 500))

  theta <- unname(fit$par[seq_len(n_zeta)])
  zeta <- if (n_zeta == 1) theta else cumsum(c(theta[1], exp(theta[-1])))
  b <- if (estimate_b) unname(fit$par[n_zeta + 1]) else NA_real_
  SE_b <- if (estimate_b) unname(sqrt(solve(fit$hessian)[n_zeta + 1, n_zeta + 1])) else NA_real_

  list(b = b, SE_b = SE_b, zeta = unname(zeta))
}
