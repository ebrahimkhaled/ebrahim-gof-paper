## realdata_screen.R -- every candidate data set in the installed packages, the same analysis for each, all
## reported: the linear logit model a textbook uses for it, the same model with the complementary log-log
## link (AIC difference: positive = cloglog fits better, the signature of an asymmetric link), and the
## HL, EF and Stukel tests on the logit fit at ten groups. Nothing is dropped from the output.
source("./code/ef_exact.R")
suppressMessages({ library(aplore3); library(MASS) })
expand <- function(d, dead, total) {                     # binomial counts -> one row per subject
  i <- rep(seq_len(nrow(d)), d[[total]])
  y <- unlist(lapply(seq_len(nrow(d)), function(k) rep(c(1, 0), c(d[[dead]][k], d[[total]][k] - d[[dead]][k]))))
  cbind(d[i, , drop = FALSE], y = y)
}
data("fbeetle", package = "VGAM"); data("budworm", package = "doBy")
cand <- list(
  beetle   = list(d = expand(fbeetle, "dead", "n"), f = y ~ logdose,
                  src = "Bliss (1935) flour-beetle mortality, VGAM::fbeetle"),
  budworm  = list(d = transform(expand(budworm, "ndead", "ntotal"), ldose = log2(dose)), f = y ~ sex + ldose,
                  src = "Collett (2003) tobacco budworm, doBy::budworm"),
  icu      = list(d = transform(icu, y = as.integer(sta == "Died")), f = y ~ age + can + cpr + inf + sys + type + loc,
                  src = "Hosmer et al. (2013) ICU, aplore3::icu"),
  glow500  = list(d = transform(glow500, y = as.integer(fracture == "Yes")),
                  f = y ~ age + weight + priorfrac + momfrac + armassist + raterisk,
                  src = "Hosmer et al. (2013) GLOW, aplore3::glow500"),
  burn1000 = list(d = transform(burn1000, y = as.integer(death == "Dead")), f = y ~ age + tbsa + race + inh_inj + flame,
                  src = "Hosmer et al. (2013) burn injury, aplore3::burn1000"),
  lowbwt   = list(d = transform(lowbwt, y = as.integer(low == "< 2500 g")), f = y ~ age + lwt + race + smoke + ptl + ht + ui,
                  src = "Hosmer et al. (2013) low birth weight, aplore3::lowbwt"),
  pima     = list(d = transform(rbind(Pima.tr, Pima.te), y = as.integer(type == "Yes")),
                  f = y ~ npreg + glu + bp + skin + bmi + ped + age, src = "Pima diabetes, MASS::Pima.tr + Pima.te")
)
res <- lapply(names(cand), function(nm) {
  c0 <- cand[[nm]]
  fl <- glm(c0$f, binomial("logit"), data = c0$d)
  fc <- glm(c0$f, binomial("cloglog"), data = c0$d)
  X  <- model.matrix(fl); y <- fl$y; p <- fitted(fl)
  t  <- ef_test(y, p, X); s <- stukel_joint(y, p, X, fl$linear.predictors)
  data.frame(data = nm, n = length(y), events = sum(y), AIC_logit_minus_cloglog = round(AIC(fl) - AIC(fc), 2),
             HL = round(t["HL"], 2), C = round(t["C"], 2), EF = round(t["EF"], 2),
             p_HL = signif(t["p_HL_chisq"], 3), p_EF = signif(t["p_EF_chisq"], 3),
             p_EF_limit_sumlam = round(t["sum_lam"], 2), p_Stukel = signif(s["p_Stukel"], 3), source = c0$src,
             row.names = NULL)
})
R <- do.call(rbind, res)
options(width = 250); print(R)
write.csv(R, "./results/realdata_screen.csv", row.names = FALSE)
