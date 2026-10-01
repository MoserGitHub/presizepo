#' Build experimental-arm category probabilities
#'
#' Given control-arm outcome-category probabilities, builds the
#' experimental-arm probabilities implied by a proportional (cumulative)
#' odds ratio `OR`.
#'
#' @param pc Numeric vector of control-arm category probabilities, summing
#'   to 1, ordered from the lowest to the highest category.
#' @param OR Proportional (cumulative) odds ratio for the experimental vs.
#'   control arm.
#'
#' @return A list with elements `p_ctrl` and `p_trt`, the control- and
#'   experimental-arm category probability vectors.
#' @noRd
#' @examples
#' build_arm_probs(c(0.15, 0.10, 0.10, 0.10, 0.10, 0.45), OR = 1.7)
build_arm_probs <- function(pc, OR) {
  cum_c <- cumsum(pc)[-length(pc)]
  odds_t <- (cum_c / (1 - cum_c)) * OR
  p_trt <- diff(c(0, odds_t / (1 + odds_t), 1))
  list(p_ctrl = pc, p_trt = p_trt)
}

#' Standard error of the log odds ratio, binary predictor
#'
#' Computes the standard error (and confidence interval) of the log odds
#' ratio for a binary (e.g. treatment) predictor in a proportional-odds
#' model, at a given total sample size `n`.
#'
#' Two approaches to the SE are supported, selectable via `method`:
#' * `"whitehead"` - Whitehead's (1993) formula, originally derived for
#'   hypothesis testing with Type I error \eqn{\alpha} and power \eqn{1-\beta}
#'   via \eqn{n = (z_{\alpha/2} + z_\beta)^2 / \text{Information}}. For
#'   precision-based calculations, the same information structure is used to
#'   compute the standard error: \eqn{\text{SE}(\log \text{OR}) = \sqrt{\frac{3(r+1)^2}{rn(1 - \sum_k \bar{p}_k^3)}}}
#'   where \eqn{\bar{p}_k = (p_{c,k} + r \cdot p_{t,k}) / (1 + r)} is the
#'   weighted average probability for category \eqn{k}. The denominator
#'   \eqn{(1 - \sum_k \bar{p}_k^3)} is the information per subject from the
#'   ordinal outcome, which measures discriminatory power and depends only on
#'   the outcome distribution, not on \eqn{\alpha} or \eqn{\beta}.
#' * `"ologit"` (default) - White, Marley-Zagar, Morris, Parmar, Royston & Babiker
#'   (2023) fit a weighted proportional-odds model to the
#'   anticipated distribution and use the SE from the observed information
#'   matrix. This approach avoids the delta-method approximation and directly
#'   accounts for joint estimation of intercepts and slope.
#'
#' @param pc Control-arm category probabilities (sum to 1).
#' @param OR Common odds ratio (experimental vs
#'   control arm).
#' @param n Total sample size (both arms).
#' @param r Allocation ratio, treatment:control.
#' @param method `"ologit"` (default) or `"whitehead"`
#' @param conf Confidence level.
#'
#' @return A list with elements `conf`, `se_log_OR`, `OR`, `lower_OR`,
#'   `upper_OR` (CI limits for the OR), and `method`.
#' @export
#'
#' @examples
#' pc <- c(0.15, 0.10, 0.10, 0.10, 0.10, 0.45)
#' se_binary_po(pc, OR = 1.7, n = 30, method = "whitehead")
#' se_binary_po(pc, OR = 1.7, n = 30, method = "ologit")
se_binary_po <- function(pc, OR, n, r = 1,
                         method = c("ologit", "whitehead"), conf = 0.95) {
  method <- match.arg(method)
  z <- stats::qnorm(1 - (1 - conf) / 2)

  probs <- build_arm_probs(pc, OR)
  p_ctrl <- probs$p_ctrl
  p_trt <- probs$p_trt
  levels_n <- length(p_ctrl)
  level_int <- rep(seq_len(levels_n), 2)
  x <- c(rep(0, levels_n), rep(1, levels_n))
  w <- n * c(p_ctrl / (1 + r), r * p_trt / (1 + r))

  need_fit_alt <- is.null(OR) || method == "ologit"
  if (need_fit_alt) {
    fit_alt <- fit_cumlogit_weighted(level_int, x, w)
    # sign flipped
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

  list(conf = conf, se_log_OR = se, OR = exp(logOR), lower_OR = exp(logOR - z * se), upper_OR = exp(logOR + z * se),
       method = method)
}

#' Sample size for a target common OR precision, binary predictor
#'
#' Computes the total sample size `n` needed so that the confidence
#' interval for the odds ratio has a target ratio of upper to lower limit,
#' for a binary predictor in a proportional-odds model. Since the Fisher
#' information scales linearly with `n`, the SE at any `n` is
#' `SE(n0) * sqrt(n0 / n)`; the fit is done once
#' at `n0` and `n` is rescaled.
#'
#' @param pc Control-arm category probabilities (sum to 1).
#' @param OR Proportional (cumulative) odds ratio for the experimental vs.
#'   control arm.
#' @param ratio_UL Target ratio of the upper to lower confidence limit for
#'   the OR.
#' @param r Allocation ratio, treatment:control.
#' @param conf Confidence level for the reported interval.
#' @param method `"ologit"` (default; exact MLE) or `"whitehead"`
#'   (delta-method approximation).
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
n_precision_binary_po <- function(pc, OR, ratio_UL, r = 1,
                                  conf = 0.95, method = c("ologit", "whitehead"), n0 = 1000) {
  method <- match.arg(method)
  z <- stats::qnorm(1 - (1 - conf) / 2)
  target_se <- log(ratio_UL) / (2 * z)
  fit0 <- se_binary_po(pc, OR, n = n0, r = r, method = method)
  n <- ceiling(n0 * (fit0$se_log_OR / target_se)^2)
  list(n = n, conf = conf, OR = fit0$OR, method = method,
       lci_OR = exp(log(fit0$OR) - z * target_se), uci_OR = exp(log(fit0$OR) + z * target_se))
}


#' Resolve the per-unit log OR for a continuous predictor (internal)
#'
#' Converts either `OR` (log OR per raw unit of `x`) or `OR_sd` (log OR per SD of `x`)
#' into the per-unit log OR used internally by [se_cont_po()] and [n_precision_cont_po()].
#'
#' @param OR Log OR per unit of `x`, or `NULL`.
#' @param OR_sd Log OR per SD of `x`, or `NULL`.
#' @param sd_x SD of the predictor `x`.
#'
#' @return The per-unit log OR (a single number).
#' @noRd
resolve_beta_cont <- function(OR, OR_sd, sd_x) {
  n_specified <- sum(!is.null(OR), !is.null(OR_sd))
  if (n_specified != 1) stop("please specify exactly one of OR, OR_sd")
  if (is.null(OR)) OR <- OR_sd / sd_x
  OR
}

#' Resolve delta (in raw units of x) for a continuous predictor (internal)
#'
#' `delta` (in raw units of `x`) and `delta_sd` (in SDs of `x`) are two
#' ways of specifying the same change in `x`; at most one may be supplied.
#' This converts either into the raw-unit `delta` used internally by
#' [se_cont_po()] and [n_precision_cont_po()], defaulting to 1 (raw unit)
#' if neither is supplied.
#'
#' @param delta Change in `x`, in raw units, or `NULL`.
#' @param delta_sd Change in `x`, in SDs, or `NULL`.
#' @param sd_x SD of the predictor `x`.
#'
#' @return The change in `x` in raw units (a single number).
#' @noRd
resolve_delta_cont <- function(delta, delta_sd, sd_x) {
  n_specified <- sum(!is.null(delta), !is.null(delta_sd))
  if (n_specified > 1) stop("please specify at most one of delta, delta_sd")
  if (!is.null(delta_sd)) return(delta_sd * sd_x)
  if (is.null(delta)) return(1)
  delta
}

#' Standard error of the log OR per unit of a continuous predictor
#'
#' Computes the standard error (and confidence interval) of the log odds
#' ratio per unit of a continuous predictor `x` in a proportional-odds
#' model, at a given total sample size `n`.
#'
#' @param p0 Outcome category probabilities at `x = 0`.
#' @param OR Odds ratio per unit of `x`. Exactly one of `OR`, `OR_sd` must
#'   be supplied.
#' @param OR_sd Odds ratio per SD of `x`, given as an alternative to `OR` so
#'   `sd_x` and the effect size don't need to be reconciled by hand;
#'   internally converted to `log(OR) = log(OR_sd) / sd_x`. Exactly one of
#'   `OR`, `OR_sd` must be supplied.
#' @param sd_x SD of the predictor `x`. `x` is assumed to be centered (mean
#'   0) and Gaussian, `x ~ N(0, sd_x^2)`.
#' @param n Total sample size.
#' @param r Allocation ratio, treatment:control (default 1 for equal allocation).
#'   Affects the effective information available for estimating the slope, just
#'   as in the binary case.
#' @param R2 Proportion of variance of `x` explained by other covariates
#'   (their inclusion reduces the effective information for `beta` by a
#'   factor of `1 - R2`).
#' @param method `"ologit"` (default) jointly estimates the cutpoints and
#'   `beta` by weighted MLE over `x ~ N(0, sd_x^2)`. The algorithm integrates
#'   the anticipated outcome probabilities over the predictor distribution
#'   using Gauss-Hermite quadrature (with `ngrid` nodes), constructs weighted
#'   pseudo-observations, fits a weighted proportional-odds model, and extracts
#'   the standard error from the inverse Hessian matrix. This approach avoids
#'   delta-method approximation and directly accounts for joint estimation of
#'   cutpoints and slope. `"whitehead"` evaluates a closed-form variance
#'   formula (extended from Whitehead's 1993 formula to continuous predictors):
#'   \eqn{\text{SE}(\beta) = \sqrt{\frac{3}{n \sigma_x^2 (1 - \sum_k \bar{p}_k^3) (1 - R^2)}}}
#'   where \eqn{\bar{p}_k} are outcome category probabilities averaged over the
#'   distribution of `x` (computed via numerical integration over a 2000-point
#'   quantile grid). The information term \eqn{(1 - \sum_k \bar{p}_k^3)} reflects
#'   the ordinal precision of the outcome distribution.
#' @param ngrid Number of Gauss-Hermite quadrature nodes (only used for
#'   `method = "ologit"`).
#' @param conf Confidence level for the reported interval.
#' @param delta Change in `x`, in raw units, over which the OR (and its
#'   CI) is expressed; `beta`, `SE`, and the CI limits are all scaled by
#'   `delta` before exponentiating (the linear predictor is `beta * x`, so
#'   this rescaling is exact). At most one of `delta`, `delta_sd` may be
#'   supplied; defaults to 1 (raw unit of `x`) if neither is given. Note
#'   `delta`/`delta_sd` only change what change in `x` the *output* is
#'   expressed over; they do not alter `OR`/`OR_sd` itself, which is
#'   always a per-1-unit/per-1-SD effect (see Details).
#' @param delta_sd Change in `x`, in SDs of `x`, given as an alternative to
#'   `delta` so `sd_x` doesn't need to be multiplied in by hand; internally
#'   converted to `delta = delta_sd * sd_x`. At most one of `delta`,
#'   `delta_sd` may be supplied.
#'
#' @details
#' **Calculation Methods**:
#'
#' * **Whitehead**: Extends Whitehead's (1993) information-based formula to
#'   continuous predictors. The standard error is
#'   \eqn{\text{SE}(\beta) = \sqrt{\frac{3(r+1)^2}{rn \sigma_x^2 (1 - \sum_k \bar{p}_k^3) (1 - R^2)}}}
#'   where \eqn{r} is the allocation ratio (treatment:control), which reduces
#'   effective information as it deviates from 1. The information term
#'   \eqn{(1 - \sum_k \bar{p}_k^3)} is derived from Whitehead's hypothesis-testing
#'   framework but applied here to precision-based calculations. Average outcome
#'   probabilities \eqn{\bar{p}_k} are computed by numerical integration (quantile
#'   grid method with 2000 points) over `x ~ N(0, sd_x^2)`, and \eqn{R^2} adjusts
#'   for confounding by other covariates.
#'
#' * **ologit**: Fits a weighted proportional-odds model to pseudo-observations
#'   created by: (1) evaluating outcome probabilities at Gauss-Hermite
#'   quadrature nodes scaled to `x ~ N(0, sd_x^2)`, (2) weighting by node
#'   weight × sample size × (1 - R2), and (3) extracting the SE of the slope
#'   coefficient from the inverse Hessian. This approach avoids delta-method
#'   approximation and is more accurate for large effect sizes and skewed
#'   distributions.
#'
#' `x` is assumed to be centered (mean 0) and Gaussian, `x ~ N(0, sd_x^2)`;
#' `p0` is interpreted as the outcome-category probabilities at `x = 0`.
#'
#' `OR_sd` and `delta_sd` both reference the SD of `x`, but play different
#' roles and are easy to conflate: `OR_sd` is *always* the OR for a
#' 1-SD change (it fixes `OR`), while `delta_sd` says how many SDs the
#' *reported* OR should span. Because the log OR is linear in `x`, the
#' reported OR for a `delta_sd`-SD change is `OR_sd^delta_sd`, not `OR_sd`
#' itself. For example, `OR_sd = 1.5` with `delta_sd = 2` reports an OR of
#' `1.5^2 = 2.25` (the OR for a 2-SD change), not `1.5` (see final example
#' below). To report an OR of 1.5 *for* a 2-SD change, supply
#' `OR_sd = sqrt(1.5)` together with `delta_sd = 2`.
#'
#' @return A named numeric vector with elements `se_log_OR`, `OR`, `lower_OR`,
#'   `upper_OR` (the last two on the OR scale), all expressed for the change
#'   in `x` given by `delta`/`delta_sd`.
#' @export
#'
#' @examples
#' p0 <- c(0.10, 0.15, 0.15, 0.20, 0.20, 0.10, 0.10)
#' se_cont_po(p0, OR = 1.5, sd_x = 0.5, n = 30, method = "whitehead")
#' se_cont_po(p0, OR = 1.5, sd_x = 0.5, n = 30, method = "ologit")
#' # equivalently, specifying the effect as OR per SD of x
#' se_cont_po(p0, OR_sd = 1.5, sd_x = 0.5, n = 30, method = "whitehead")
#' # OR/CI for a 2-unit (here, 2-SD, since sd_x = 1) change in x
#' se_cont_po(p0, OR = 1.5, sd_x = 1, n = 30, delta = 2, method = "whitehead")
#' # equivalently, specifying the change directly in SDs of x
#' se_cont_po(p0, OR_sd = 1.5, sd_x = 1, n = 30, delta_sd = 2, method = "whitehead")
#'
#' # OR_sd is per 1 SD; delta_sd only rescales the *output*, so with
#' # delta_sd = 2 the reported OR is OR_sd^2, not OR_sd:
#' fit <- se_cont_po(p0, OR_sd = 1.5, sd_x = 1, n = 30, delta_sd = 2, method = "whitehead")
#' fit[["OR"]]  # 2.25 = 1.5^2, the OR for a 2-SD change, not 1.5
#' # to instead report an OR of 1.5 *for* a 2-SD change, use sqrt(1.5):
#' fit2 <- se_cont_po(p0, OR_sd = sqrt(1.5), sd_x = 1, n = 30, delta_sd = 2, method = "whitehead")
#' fit2[["OR"]]  # 1.5, as intended
se_cont_po <- function(p0, OR = NULL, sd_x, n, r = 1, R2 = 0, OR_sd = NULL,
                       method = c("ologit", "whitehead"), ngrid = 15, conf = 0.95,
                       delta = NULL, delta_sd = NULL) {
  method <- match.arg(method)
  # Store input OR/OR_sd before converting to log scale
  OR_input <- OR
  OR_sd_input <- OR_sd
  # Convert OR and OR_sd to log scale
  if (!is.null(OR)) OR <- log(OR)
  if (!is.null(OR_sd)) OR_sd <- log(OR_sd)
  log_beta <- resolve_beta_cont(OR, OR_sd, sd_x)
  delta <- resolve_delta_cont(delta, delta_sd, sd_x)
  z <- stats::qnorm(1 - (1 - conf) / 2)
  levels_n <- length(p0)
  th <- stats::qlogis(cumsum(p0)[-levels_n])

  if (method == "whitehead") {
    ngrid_wh <- 2000
    xs <- stats::qnorm((1:ngrid_wh - 0.5) / ngrid_wh) * sd_x  # quantile grid for x ~ N(0, sd_x)
    pm <- sapply(xs, function(x) diff(c(0, stats::plogis(th - log_beta * x), 1)))  # levels_n x ngrid_wh
    pmbar <- rowMeans(pm)
    se <- sqrt(3 * (r + 1)^2 / (r * n * sd_x^2 * (1 - sum(pmbar^3)) * (1 - R2)))
  } else {
    gh <- statmod::gauss.quad(ngrid, kind = "hermite")
    xs <- sqrt(2) * sd_x * gh$nodes            # nodes rescaled to x ~ N(0, sd_x^2)
    wq <- gh$weights / sqrt(pi)                 # mixture weights, sum to 1
    pm <- sapply(xs, function(x) diff(c(0, stats::plogis(th - log_beta * x), 1)))  # levels_n x ngrid

    level_int <- rep(seq_len(levels_n), times = ngrid)
    x_rep <- rep(xs, each = levels_n)
    w <- as.vector(pm) * rep(wq, each = levels_n) * n * (1 - R2)
    fit <- fit_cumlogit_weighted(level_int, x_rep, w)
    se <- fit$SE_b
  }

  log_beta_d <- delta * log_beta
  se_d <- delta * se
  c(se_log_OR = se_d, OR = exp(log_beta_d), lower_OR = exp(log_beta_d - z * se_d), upper_OR = exp(log_beta_d + z * se_d))
}

#' Sample size for a target common OR precision, continuous predictor
#'
#' Computes the total sample size `n` needed so that the confidence
#' interval for the OR per `delta` units of a continuous predictor `x` has
#' a target ratio of upper to lower limit. It uses a Gauss-Hermite quadrature
#' approximation for the `ologit` method, see details.
#'
#' @inheritParams se_cont_po
#' @param ratio_UL Target ratio of the upper to lower confidence limit for
#'   the OR per `delta`/`delta_sd` change in `x`.
#' @param delta Change in `x`, in raw units, over which the OR (and its
#'   CI) is expressed. At most one of `delta`, `delta_sd` may be supplied;
#'   defaults to 1 (raw unit of `x`) if neither is given.
#' @param n0 Sample size at which the one-time reference fit is done before
#'   rescaling.
#' @return A list with elements `n`, `or`, `method`, `lci`, `uci`. `or`,
#'   `lci`, and `uci` are expressed for the change in `x` given by
#'   `delta`/`delta_sd` (not a 1-raw-unit change), matching the scale on
#'   which `ratio_UL` was specified.
#' @export
#'
#' @examples
#' p0 <- c(0.10, 0.15, 0.15, 0.20, 0.20, 0.10, 0.10)
#' n_precision_cont_po(p0, OR = 1.5, sd_x = 0.5, ratio_UL = 3, method = "whitehead")
#' n_precision_cont_po(p0, OR = 1.5, sd_x = 0.5, ratio_UL = 3, method = "ologit")
#' # equivalently, specifying the effect as OR per SD of x
#' n_precision_cont_po(p0, OR_sd = 1.5, sd_x = 0.5, ratio_UL = 3, method = "whitehead")
#' # OR of 1.5 per 2-SD change in x, specifying the change directly in SDs
#' n_precision_cont_po(p0, OR_sd = 1.5, sd_x = 0.5, delta_sd = 2, ratio_UL = 3, method = "whitehead")
n_precision_cont_po <- function(p0, OR = NULL, sd_x, ratio_UL, delta = NULL, r = 1, R2 = 0, conf = 0.95,
                                OR_sd = NULL, delta_sd = NULL, method = c("ologit", "whitehead"),
                                ngrid = 15, n0 = 1000) {
  method <- match.arg(method)
  # Convert OR and OR_sd to log scale
  if (!is.null(OR)) OR <- log(OR)
  if (!is.null(OR_sd)) OR_sd <- log(OR_sd)
  log_beta <- resolve_beta_cont(OR, OR_sd, sd_x)
  delta <- resolve_delta_cont(delta, delta_sd, sd_x)
  z <- stats::qnorm(1 - (1 - conf) / 2)
  target_se <- log(ratio_UL) / (2 * z * delta)              # target SE per 1 unit of x
  fit0 <- se_cont_po(p0, OR = log_beta, sd_x = sd_x, n = n0, r = r, R2 = R2, method = method, ngrid = ngrid)
  n <- ceiling(n0 * (fit0[["se_log_OR"]] / target_se)^2)
  list(n = n, OR = exp(delta * log_beta), method = method,
       lci_OR = exp(delta * (log_beta - z * target_se)), uci_OR = exp(delta * (log_beta + z * target_se)))
}


#' Weighted cumulative-logit MLE (internal)
#'
#' Fits a proportional-odds (cumulative logit) model to a weighted table of
#' outcome-category counts, generic in `x` (a binary group indicator or a
#' continuous covariate). Used from artcatr.r
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
