## rivals.R -- the competing tests of the rival studies, taken from the EDGE battery library so that both papers use
## the same validated implementations (simulations/_battery_tests.R: line-by-line ports of ebrahim.gof 2.8.0, checked
## against the package to 1e-8). The grouped Farrington statistic uses ef_exact.R with random tie-breaking.
SIMDIR <- "./code/lib"
source(file.path(SIMDIR, "_battery_tests.R"))
source("./code/ef_exact.R")

## every test on one fitted data set; y, X (with intercept) as used in the fit
rival_tests <- function(y, X, giviti = TRUE) {
  d   <- data.frame(y = y, X[, -1, drop = FALSE])
  f   <- as.formula(paste("y ~", paste(colnames(X)[-1], collapse = " + ")))
  fq  <- bt_fit(list(d = d, f = f))
  if (is.null(fq) || !isTRUE(fq$fit$converged)) return(NULL)
  o   <- sample.int(length(y))                          # random tie-breaking for the rank groups
  ef  <- ef_test(y[o], fq$p_raw[o], fq$X[o, , drop = FALSE], 10)
  ef5 <- ef_test(y[o], fq$p_raw[o], fq$X[o, , drop = FALSE], 5)
  ef20 <- ef_test(y[o], fq$p_raw[o], fq$X[o, , drop = FALSE], 20)
  stk <- bt_stukel(fq)
  ## Stukel-W (Addendum 3): the added directions built from the linear predictor winsorized at its 2.5% and 97.5%
  ## sample quantiles, so one extreme record cannot carry the score
  fw <- fq; qe <- stats::quantile(fq$eta, c(.025, .975), names = FALSE); fw$eta <- pmin(pmax(fq$eta, qe[1]), qe[2])
  jw <- bt_stukel_joint_stat(fw)
  stkW <- if (is.finite(jw$chi)) stats::pchisq(jw$chi, jw$df, lower.tail = FALSE) else NA_real_
  c(EF = unname(ef["p_EF_chisq"]), HL = unname(ef["p_HL_chisq"]),
    EF_stat = unname(ef["EF"]), HL_stat = unname(ef["HL"]), C = unname(ef["C"]),
    EF5 = unname(ef5["p_EF_chisq"]), HL5 = unname(ef5["p_HL_chisq"]),
    EF20 = unname(ef20["p_EF_chisq"]), HL20 = unname(ef20["p_HL_chisq"]),
    PH = bt_ph(fq$y, fq$ph, 10), Tsiatis = bt_tsiatis(fq, 10), Xie = bt_xie(fq),
    Stukel = unname(stk["Stk.joint"]), StukelW = stkW, Cubic = bt_cubic(fq$y, fq$eta),
    GiViTI = if (giviti) bt_giviti(fq$y, fq$p_raw, 0.95) else NA_real_)
}
