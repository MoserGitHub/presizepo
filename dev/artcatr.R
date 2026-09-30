# artcatr: sample size / power for a two-arm trial with an ordinal
# outcome
#
# R implementation of White, Marley-Zagar, Morris, Parmar, Royston & Babiker
# (2023) "artcat: Sample-size calculation for an ordered categorical outcome",
# https://doi.org/10.1177/1536867X231161934

#' Build level probabilities (internal)
#'
#' @description Internal function for calculatin level probabilities 
#' for control and experimental arms based on input parameters
#' @param pc Level probabilities for control group (cumulative or
#'   non-cumulative, depending on `cumulative`), including the right-most
#'   level (probabilities must sum to 1, or cumulative probabilities must
#'   reach 1)
#' @param cumulative Logical; if `TRUE`, `pc` (and `pe`, if supplied) are
#'   cumulative probabilities, otherwise they are per-level probabilities
#' @param pe Level probabilities for experimental group, or `NULL` to derive
#'   them from `or` or `rr`
#' @param or Anticipated odds ratio (experimental vs control), or `NULL`
#' @param rr Anticipated risk ratio (experimental vs control), or `NULL`
#' @param margin Non-inferiority/equivalence margin, expressed as an odds
#'   ratio applied to the control probabilities
#' @returns A list with components:
#'   p1: Non-cumulative level probabilities for the control group
#'   p2: Non-cumulative level probabilities for the experimental group
#'   p3: Non-cumulative level probabilities for the margin-shifted
#'     control group
#'   levels: Number of levels remaining after dropping any levels with
#'     zero probability in both arms
#' @noRd

.artcatr_build_probs <- function(pc, cumulative, pe, or, rr, margin) {
  p1sum <- if (cumulative) pc else cumsum(pc)

  n_specified <- sum(!is.null(pe), !is.null(or), !is.null(rr))
  if (n_specified == 0) {
    if (margin != 1) {
      message("Note: assuming anticipated odds ratio = 1")
      or <- 1
    } else {
      stop("please specify one of pe, or, rr")
    }
  } else if (n_specified > 1) {
    stop("please don't specify more than one of pe, or, rr")
  }

  if (!is.null(pe)) {
    p2sum <- if (cumulative) pe else cumsum(pe)
    if (length(p2sum) != length(p1sum)) stop("pc and pe have different lengths")
  } else if (!is.null(or)) {
    p2sum <- or * p1sum / (1 - p1sum + or * p1sum)
  } else if (!is.null(rr)) {
    p2sum <- rr * p1sum
    if (isTRUE(all.equal(p1sum[length(p1sum)], 1))) p2sum[length(p2sum)] <- 1
    if (any(p2sum > 1 + 1e-8)) stop("rr option implies pe sums to more than 1")
  }

  if (any(p1sum > 1 + 1e-8)) stop("probabilities in pc sum to more than 1")
  if (any(p2sum > 1 + 1e-8)) stop("probabilities in pe sum to more than 1")

  p3sum <- margin * p1sum / (1 - p1sum + margin * p1sum)

  if (!isTRUE(all.equal(p1sum[length(p1sum)], 1))) {
    stop("pc must include the right-most level's probability (cumulative probabilities must reach 1)")
  }
  if (!isTRUE(all.equal(p2sum[length(p2sum)], 1))) {
    stop("pe must include the right-most level's probability (cumulative probabilities must reach 1)")
  }

  p1 <- diff(c(0, p1sum))
  p2 <- diff(c(0, p2sum))
  p3 <- diff(c(0, p3sum))

  drop <- (p1 == 0) & (p2 == 0)
  if (any(drop)) {
    message(sprintf(
      "Warning: %d level(s) dropped due to zero probability in both arms",
      sum(drop)
    ))
    p1 <- p1[!drop]; p2 <- p2[!drop]; p3 <- p3[!drop]
  }

  if (any(p1 < -1e-8) || any(p2 < -1e-8)) {
    stop(if (cumulative) {
      "decreasing cumulative probabilities found in pc/pe"
    } else {
      "negative probabilities found in pc/pe"
    })
  }
  p1 <- pmax(p1, 0); p2 <- pmax(p2, 0); p3 <- pmax(p3, 0)

  list(p1 = p1, p2 = p2, p3 = p3, levels = length(p1))
}

#' Fit a weighted cumulative-logit model (internal)
#'
#' @description Weighted maximum-likelihood fit of a cumulative-logit
#' model for an ordered outcome with a binary covariate `x`, allowing
#' arbitrary case weights and an optional fixed offset.
#' @param level_int Integer vector of ordered outcome levels (coded
#'   `1:levels_n`) for each observation
#' @param x Binary covariate value (here treatment group indicator)
#' @param w Numeric vector of case weights for each observation
#' @param offset Numeric vector of fixed offsets added to the linear
#'   predictor for each observation (needed for non-inferiority margin)
#' @param estimate_b Logical; if `TRUE`, the coefficient `b` on `x` is
#'   estimated, otherwise it is fixed at zero and only the cutpoints are
#'   estimated
#' @returns A list with components:
#'   b: Estimated coefficient on `x` (the log odds ratio), or `NA` if
#'     `estimate_b` is `FALSE`
#'   SE_b: Standard error of `b`, or `NA` if `estimate_b` is `FALSE`
#'   zeta: Numeric vector of estimated cutpoints (length
#'     `levels_n - 1`)
#' @examples The following example shows that
#'  .artcatr_fit_cumlogit() reproduces the coefficient, SE, and cutpoints
#' from a weighted proportional-odds model fit via MASS::polr(),
#'  Small discrepancy because of different optimizers.
#' 
#' levels_n <- 4
#' level_int <- rep(seq_len(levels_n), 2)
#' x <- c(rep(0, levels_n), rep(1, levels_n))
#' p1 <- c(0.3, 0.3, 0.2, 0.2)
#' p2 <- c(0.15, 0.25, 0.3, 0.3)
#' w <- c(p1 / 2, p2 / 2)
#'
#' fit <- artcatr:::.artcatr_fit_cumlogit(level_int, x, w)
#' fit$b     # 0.7510
#' fit$SE_b  # 3.6513
#' fit$zeta  # -0.8948, 0.3871, 1.5188
#'
#' df <- data.frame(
#'   y = factor(level_int, levels = seq_len(levels_n), ordered = TRUE),
#'   x = x
#' )
#' fit_polr <- MASS::polr(y ~ x, data = df, weights = w, Hess = TRUE)
#' coef(fit_polr)             # 0.7505
#' sqrt(vcov(fit_polr))[1, 1] # 3.6513
#' fit_polr$zeta              # -0.8953, 0.3867, 1.5184
#' @noRd
.artcatr_fit_cumlogit <- function(level_int, x, w, offset = rep(0, length(x)),
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

#' Fit the ologit null/alternative log odds ratio (internal)
#'
#' @description Fit the ologit-based null/alternative log-OR and its two SEs
#' according to "ologit" method of White et al (2023, sec. 2.2): fit a
#' proportional-odds model to the anticipated joint distribution to get the
#' alternative-hypothesis log OR and SE (SA); then fit a null-constrained
#' model (OR fixed at margin) to get null-consistent category probabilities,
#' and refit to get the null SE (SN).
#'
#' @param p1 Numeric vector of non-cumulative anticipated level probabilities
#'   for the control arm
#' @param p2 Numeric vector of non-cumulative anticipated level probabilities
#'   for the experimental arm
#' @param A Numeric scalar allocation ratio (control : experimental)
#' @param margin Numeric scalar odds ratio defining the null hypothesis
#'   (`margin = 1` for a superiority trial)
#' @returns A list with components:
#'   logor: Null-hypothesis log odds ratio, i.e. `b_exptrt + log(margin)`
#'   SA: Standard error of the alternative-hypothesis log odds ratio
#'   SN: Standard error of the log odds ratio under the null-constrained
#'     (margin) model
#'   orcalc: Alternative-hypothesis odds ratio implied by `p1` and `p2`
#' @examples Example that .artcatr_ologit_fit reproduces results from MASS::polr
#' 
#' The alternative-hypothesis fit (logor/SA, with margin = 1) reproduces
#' the coefficient and SE from a weighted proportional-odds model fit
#' via MASS::polr() on the same weighted pseudo-data
#' p1 <- c(0.3, 0.3, 0.2, 0.2)
#' p2 <- c(0.15, 0.25, 0.3, 0.3)
#' fit <- artcatr:::.artcatr_ologit_fit(p1, p2, A = 1, margin = 1)
#' fit$logor # 0.751
#' fit$SA    # 3.651
#'
#' levels_n <- length(p1)
#' level_int <- rep(seq_len(levels_n), 2)
#' exptrt <- c(rep(0, levels_n), rep(1, levels_n))
#' w <- c(p1 / 2, p2 / 2)
#' df <- data.frame(
#'   y = factor(level_int, levels = seq_len(levels_n), ordered = TRUE),
#'   x = exptrt
#' )
#' fit_polr <- MASS::polr(y ~ x, data = df, weights = w, Hess = TRUE)
#' coef(fit_polr)               # 0.751
#' sqrt(vcov(fit_polr))[1, 1]   # 3.651
#' @noRd
.artcatr_ologit_fit <- function(p1, p2, A, margin) {
  levels_n <- length(p1)
  level_int <- rep(seq_len(levels_n), 2)
  exptrt <- c(rep(0, levels_n), rep(1, levels_n))
  w <- c(p1 * A / (A + 1), p2 * 1 / (A + 1))

  # Alternative-hypothesis fit: gives logOR and its SE (SA)
  fit_alt <- .artcatr_fit_cumlogit(level_int, exptrt, w, estimate_b = TRUE)
  
  b_exptrt <- fit_alt$b
  SA <- fit_alt$SE_b
  orcalc <- exp(-b_exptrt)
  logor <- b_exptrt + log(margin)

  # Null-constrained fit — get null-consistent category probabilities
  # The treatment effect is fixed at the null value (-log(margin),
  # via a fixed offset rather than an estimated coefficient;
  # estimate_b = FALSE,
  # and only the cutpoints (zeta) are re-estimated.
  # Those null-consistent cutpoints are then used
  # to predict what the category probabilities in
  # each arm would be under the null hypothesis (OR = margin)

  offset <- -log(margin) * exptrt
  fit_null <- .artcatr_fit_cumlogit(level_int, exptrt, w, offset = offset, estimate_b = FALSE)
  zeta <- fit_null$zeta
  cum_probs <- function(eta) diff(c(0, plogis(zeta - eta), 1))
  pred <- rbind(cum_probs(0), cum_probs(-log(margin))) # row 1: control, row 2: experimental
  poff <- c(pred[1, ] * A / (A + 1), pred[2, ] * 1 / (A + 1))

  # Null SE (SN): refit the treatment-effect model against null-consistent weights
  fit_nullvar <- .artcatr_fit_cumlogit(level_int, exptrt, poff, estimate_b = TRUE)
  SN <- fit_nullvar$SE_b

  list(logor = logor, SA = SA, SN = SN, orcalc = orcalc)
}

#' Sample size / power for a randomized trial with an ordered categorical
#' outcome, analysed by the proportional-odds model
#'
#' Implements the methods of White et al (2023), including their new ologit
#' method (which allows a non-proportional-odds anticipated effect via `pe`,
#' is more accurate for large effects, and supports non-inferiority /
#' substantial-superiority trials via `margin`), and Whitehead's (1993)
#' original method for comparison.
#'
#' @param pc Numeric vector. Anticipated probabilities (or cumulative
#'   probabilities, if `cumulative = TRUE`) in the control arm, one per
#'   outcome level, including the right-most level (probabilities must sum
#'   to 1, or cumulative probabilities must reach 1).
#' @param cumulative Logical. Are `pc`/`pe` cumulative probabilities?
#' @param pe Numeric vector. Anticipated (cumulative) probabilities in the
#'   experimental arm, specified like `pc`. Exactly one of `pe`, `or`, `rr`
#'   must be given (or none, if `margin != 1`, in which case `or = 1` is
#'   assumed).
#' @param or Numeric scalar. Anticipated common odds ratio applied to the
#'   control-arm cumulative probabilities (assumes proportional odds).
#'   An odds ratio below 1 shifts probability towards the right-most level.
#' @param rr Numeric scalar. Anticipated common risk ratio (applied to every
#'   level except the right-most).
#' @param margin Numeric scalar odds ratio defining the null hypothesis.
#'   `margin = 1` (default) gives a superiority trial; otherwise a
#'   non-inferiority or substantial-superiority trial (see `favourable`).
#' @param favourable Logical or `NULL`. Is the left-most outcome level the
#'   most favourable outcome? If `NULL` (default), inferred from the
#'   anticipated average odds ratio vs `margin`.
#' @param power Numeric scalar in (0, 1). Power to achieve; sample size is
#'   calculated. Default `0.8` if neither `power` nor `n` given.
#' @param n Numeric scalar. Total sample size; power is calculated.
#' @param aratio Length-2 numeric vector, allocation ratio control:experimental
#'   (e.g. `c(1, 2)` means 2 experimental participants per 1 control).
#' @param alpha Numeric scalar significance level. Default `0.05`.
#' @param onesided Logical. Is `alpha` one-sided? Default `FALSE` (two-sided).
#' @param method One of `"ologit"` (default, the new method) or `"whitehead"`.
#' @param variance For `method = "ologit"`, one of `"NA"` (default, recommended
#'   by White et al 2023), `"NN"` (equivalent to Whitehead's method), or
#'   `"AA"`. `"NA"` uses the null-hypothesis SE for the alpha term and the
#'   alternative-hypothesis SE for the beta (power) term; `"NN"`/`"AA"` use
#'   one SE for both.
#' @param quiet Logical. Suppress the printed summary?
#'
#' @return A list with the computed sample size or power, per-group
#'   sample sizes, the anticipated probability table, and other design
#'   parameters (mirroring the Stata command's stored results).
#'
#' @references White IR, Marley-Zagar E, Morris TP, Parmar MKB, Royston P,
#'   Babiker AG (2023). artcat: Sample-size calculation for an ordered
#'   categorical outcome. Stata Journal 23(1):3-23.
#'   https://doi.org/10.1177/1536867X231161934
#'
#'   Whitehead J (1993). Sample size calculations for ordered categorical
#'   data. Statistics in Medicine 12:2257-2271.
#'
#' @examples
#' # FLU-IVIG trial (Davey et al 2019), reproducing White et al (2023) sec 4:
#' # default (new) ologit method
#' artcatr(pc = c(.018, .036, .156, .141, .39, .259), or = 1 / 1.77, power = .8,
#'        favourable = FALSE)
#'
#' # Whitehead method, for comparison
#' artcatr(pc = c(.018, .036, .156, .141, .39, .259), or = 1 / 1.77, power = .8,
#'        favourable = FALSE, method = "whitehead")
#'
#' # Non-inferiority trial (margin expressed as an odds ratio)
#' artcatr(pc = c(.010, .021, .099, .103, .384, .383), or = 1, margin = 1.33,
#'        power = .8, favourable = FALSE)
#' @export
artcatr <- function(pc,
                    cumulative = FALSE,
                    pe = NULL, or = NULL, rr = NULL,
                    margin = 1,
                    favourable = NULL,
                    power = NULL, n = NULL,
                    aratio = c(1, 1),
                    alpha = 0.05,
                    onesided = FALSE,
                    method = c("ologit", "whitehead"),
                    variance = c("NA", "NN", "AA"),
                    quiet = FALSE) {

  method <- match.arg(method)
  variance <- match.arg(variance)

  if (!xor(is.null(power), is.null(n))) {
    if (is.null(power) && is.null(n)) power <- 0.8
    else stop("please don't specify both power and n")
  }
  if (!is.null(power) && (power <= 0 || power >= 1)) stop("power must be between 0 and 1")
  if (!is.null(n) && n <= 0) stop("n must be greater than 0")
  if (alpha <= 0 || alpha >= 1) stop("alpha must be between 0 and 1")
  if (margin <= 0) stop("margin must be expressed as an odds ratio greater than 0")
  if (method == "whitehead" && is.null(or)) stop("Whitehead method requires or")
  if (method == "whitehead" && margin != 1) {
    stop("the Whitehead method is not available for non-inferiority trials")
  }

  # capture how the experimental arm was specified, for display purposes
  # (mirrors the "assume or=1" default applied inside .artcatr_build_probs)
  or_display <- or
  pe_display <- pe
  rr_display <- rr
  if (is.null(pe) && is.null(or) && is.null(rr) && margin != 1) or_display <- 1

  probs <- .artcatr_build_probs(pc, cumulative, pe, or, rr, margin)
  p1 <- probs$p1; p2 <- probs$p2; p3 <- probs$p3; levels_n <- probs$levels

  A <- aratio[1] / aratio[2]
  pbar <- (A * p1 + p2) / (A + 1)

  za <- qnorm(1 - alpha / (if (onesided) 1 else 2))
  if (!is.null(power)) zb <- qnorm(power)

  fit <- .artcatr_ologit_fit(p1, p2, A, margin)
  orcalc <- if (!is.null(or)) or else fit$orcalc

  # infer favourable/unfavourable direction if not supplied
  direction_inferred <- is.null(favourable)
  if (is.null(favourable)) favourable <- (orcalc > margin)
  if (isTRUE(all.equal(orcalc, margin))) stop("or = margin makes a trial impossible")

  trialtype <- if (margin == 1) {
    "superiority"
  } else if ((margin < 1 && favourable) || (margin > 1 && !favourable)) {
    "non-inferiority"
  } else {
    "substantial-superiority"
  }

  result <- list(
    n = NA_real_, nC = NA_real_, nE = NA_real_, power = NA_real_,
    n_whitehead = NA_real_, n_ologit_NN = NA_real_, n_ologit_NA = NA_real_, n_ologit_AA = NA_real_,
    power_whitehead = NA_real_, power_ologit_NN = NA_real_, power_ologit_NA = NA_real_, power_ologit_AA = NA_real_
  )

  round_group_sizes <- function(n_raw) {
    nE_raw <- unname(n_raw) / (1 + A)
    nC_raw <- nE_raw * A
    nE <- ceiling(nE_raw)
    nC <- ceiling(nC_raw)
    
    list(n = unname(nC + nE), nC = unname(nC), nE = unname(nE))
  }

  if (method == "whitehead") {
    sumpbar3 <- sum(pbar^3)
    if (!is.null(n)) { # n given -> compute power
      zb_calc <- sqrt(n * (log(or)^2) * (1 - sumpbar3) / 12) - za
      result$power_whitehead <- pnorm(zb_calc)
      result$power <- result$power_whitehead
    } else { # power given -> compute n
      n_raw <- 3 * (A + 1)^2 * (za + zb)^2 / (A * (log(or)^2) * (1 - sumpbar3))
      result$n_whitehead <- n_raw
      g <- round_group_sizes(n_raw)
      
      result$n <- g$n; result$nC <- g$nC; result$nE <- g$nE
    }
  } else { # ologit
    S <- c(N = fit$SN, A = fit$SA)
    for (m1 in c("N", "A")) for (m2 in c("N", "A")) {
      if (m1 == "A" && m2 == "N") next
      tag <- paste0("ologit_", m1, m2)
      if (!is.null(n)) { # n given -> compute power
        zb_calc <- (sqrt(n) * abs(fit$logor) - za * S[[m1]]) / S[[m2]]
        result[[paste0("power_", tag)]] <- pnorm(zb_calc)
      } else { # power given -> compute n
        n_raw <- (za * S[[m1]] + zb * S[[m2]])^2 / fit$logor^2
        result[[paste0("n_", tag)]] <- n_raw
      }
    }
    chosen_tag <- paste0("ologit_", variance)
    if (!is.null(n)) {
      result$power <- result[[paste0("power_", chosen_tag)]]
    } else {
      n_raw <- result[[paste0("n_", chosen_tag)]]
      g <- round_group_sizes(n_raw)
      # g <- n_raw
      result$n <- g$n; result$nC <- g$nC; result$nE <- g$nE
    }
  }

  result$probtable <- data.frame(level = seq_len(levels_n), pc = p1, pe = p2,
                                  pe_null = if (margin != 1) p3 else NA_real_)
  result$logor <- fit$logor
  result$orcalc <- orcalc
  result$favourable <- favourable
  result$direction_inferred <- direction_inferred
  result$trialtype <- trialtype
  result$method <- method
  result$variance <- variance
  result$alpha <- alpha
  result$onesided <- onesided
  result$aratio <- aratio
  result$margin <- margin
  result$cumulative <- cumulative
  result$pc_input <- pc
  result$pe_input <- pe_display
  result$or_input <- or_display
  result$rr_input <- rr_display
  result$power_designed <- if (is.null(n)) power else NA_real_
  result$n_designed <- if (!is.null(n)) n else NA_real_
  result$lefttext <- if (favourable) "most favourable" else "least favourable"
  result$righttext <- if (favourable) "least favourable" else "most favourable"
  result$ltgt <- if (favourable) ">" else "<"

  class(result) <- "artcatr"
  if (!quiet) print(result)

  invisible(result)
}


#' Print an `artcatr` object
#'
#' @param x An object returned by artcatr().
#' @param digits Number of significant digits to round probabilities to.
#' @param probtable Logical; include the table of anticipated probabilities?
#' @param ... Not used.
#' @return `x`, invisibly.
#' @export
print.artcatr <- function(x, digits = 3, probtable = TRUE, ...) {
  fmt <- function(v) format(signif(v, digits + 1), trim = TRUE)
  lab <- function(s) sprintf("%-38s", s)
  hline <- function() cat(strrep("-", 70), "\n", sep = "")
  
  cat("artcatr: sample size for an ordered categorical outcome (proportional-odds model)\n
      based on the Stata implementation from White et al. (2023)")
  hline()
  cat(lab("Type of trial"), x$trialtype, "\n", sep = "")
  cat(lab("Favourable/unfavourable outcome"),
      if (x$favourable) "favourable" else "unfavourable",
      if (x$direction_inferred) "  (inferred)" else "", "\n", sep = "")
  
  if (x$trialtype == "superiority") {
    cat(lab("Null hypothesis"), "odds ratio = 1\n", sep = "")
    cat(lab("Superiority region"), sprintf("odds ratio %s 1\n", x$ltgt), sep = "")
  } else if (x$trialtype == "non-inferiority") {
    cat(lab("Null hypothesis"), sprintf("odds ratio = %s\n", fmt(x$margin)), sep = "")
    cat(lab("Non-inferiority region"), sprintf("odds ratio %s %s\n", x$ltgt, fmt(x$margin)), sep = "")
  } else {
    cat(lab("Null hypothesis"), sprintf("odds ratio = %s\n", fmt(x$margin)), sep = "")
    cat(lab("Substantial-superiority region"), sprintf("odds ratio %s %s\n", x$ltgt, fmt(x$margin)), sep = "")
  }
  
  cat(lab("Allocation ratio C:E"), paste(x$aratio, collapse = ":"), "\n", sep = "")
  
  pc_text <- paste(fmt(x$pc_input), collapse = " ")
  cat(lab("Anticipated probabilities, control"), pc_text,
      if (x$cumulative) " (cumulative)" else "", "\n", sep = "")
  pe_text <- if (!is.null(x$pe_input)) {
    paste(fmt(x$pe_input), collapse = " ")
  } else if (!is.null(x$or_input)) {
    paste0("given by odds ratio = ", fmt(x$or_input))
  } else {
    paste0("given by risk ratio = ", fmt(x$rr_input))
  }
  cat(lab("                      experimental"), pe_text,
      if (x$cumulative && !is.null(x$pe_input)) " (cumulative)" else "", "\n", sep = "")
  if (is.null(x$or_input)) {
    cat(lab("Anticipated average odds ratio"), fmt(x$orcalc), "\n", sep = "")
  }
  
  if (probtable) {
    cat("\nTable of anticipated probabilities\n")
    has_null <- x$margin != 1
    cat(sprintf("%-4s%-20s%8s%8s%s\n", "", "", "C", "E", if (has_null) sprintf("%10s", "E null") else ""))
    pt <- x$probtable
    levels_n <- nrow(pt)
    for (i in seq_len(levels_n)) {
      levlab <- if (i == 1) x$lefttext else if (i == levels_n) x$righttext else ""
      cat(sprintf("%-4d%-20s%8.3f%8.3f%s\n", i, levlab, pt$pc[i], pt$pe[i],
                  if (has_null) sprintf("%10.3f", pt$pe_null[i]) else ""))
    }
  }
  
  cat("\n")
  cat(lab("Alpha"), sprintf("%s (%s)\n", fmt(x$alpha), if (x$onesided) "one-sided" else "two-sided"), sep = "")
  if (!is.na(x$power_designed)) cat(lab("Power (designed)"), fmt(x$power_designed), "\n", sep = "")
  if (!is.na(x$n_designed)) cat(lab("Total sample size (designed)"), x$n_designed, "\n", sep = "")
  cat(lab("Method"),
      if (x$method == "whitehead") "Whitehead" else sprintf("ologit (variance %s)", x$variance),
      "\n\n", sep = "")
  
  if (!is.na(x$n_designed)) {
    cat(lab("Power (calculated)"), fmt(x$power), "\n", sep = "")
  } else {
    cat(lab("Total sample size (calculated)"), x$n, "\n", sep = "")
    cat(lab("Sample size per group (calculated)"), sprintf("%d %d", x$nC, x$nE), "\n", sep = "")
  }
  hline()
  invisible(x)
}

#' Tidy summary table of anticipated probabilities from an `artcatr` object
#'
#' @param object An object returned by artcatr().
#' @param digits Number of significant digits to round probabilities to.
#' @param ... Not used.
#' @return A data frame with columns `level`, `description` (favourable/
#'   unfavourable labels for the extreme levels), `control`, `experimental`,
#'   and (for non-inferiority/substantial-superiority designs) `experimental_null`.
#' @export
summary.artcatr <- function(object, digits = 3, ...) {
  pt <- object$probtable
  levels_n <- nrow(pt)
  description <- rep("", levels_n)
  description[1] <- object$lefttext
  description[levels_n] <- object$righttext

  out <- data.frame(
    level = pt$level,
    description = description,
    control = round(pt$pc, digits),
    experimental = round(pt$pe, digits)
  )
  if (object$margin != 1) out$experimental_null <- round(pt$pe_null, digits)

  attr(out, "design") <- list(
    trialtype = object$trialtype, favourable = object$favourable,
    margin = object$margin, orcalc = object$orcalc, method = object$method,
    variance = object$variance, alpha = object$alpha, onesided = object$onesided,
    aratio = object$aratio, n = object$n, nC = object$nC, nE = object$nE,
    power = object$power, n_designed = object$n_designed, power_designed = object$power_designed
  )
  class(out) <- c("summary.artcatr", "data.frame")
  out
}

#' Print a `summary.artcatr` object
#'
#' @param x An object returned by [summary.artcatr()].
#' @param ... Not used.
#' @return `x`, invisibly.
#' @export
print.summary.artcatr <- function(x, ...) {
  d <- attr(x, "design")
  cat("Design:", d$trialtype,
      sprintf("(%s outcome; method = %s%s)\n",
              if (d$favourable) "favourable" else "unfavourable",
              if (d$method == "whitehead") "Whitehead" else "ologit",
              if (d$method != "whitehead") sprintf(", variance = %s", d$variance) else ""))
  if (is.na(d$n_designed)) {
    # power was given as the design target; n was solved for
    cat(sprintf("Sample size (calculated): %d total (C = %d, E = %d), power (designed) = %.3f\n",
                d$n, d$nC, d$nE, d$power_designed))
  } else {
    # n was given as the design target; power was solved for
    cat(sprintf("Sample size (designed): %d total, power (calculated) = %.3f\n",
                d$n_designed, d$power))
  }
  cat("\n")
  print.data.frame(x, row.names = FALSE)
  invisible(x)
}
