##-------------------------------------------------------------------------------
## This script implement the Chapter Bounded Count Regression model of
## the Minicourse
## Bazán, J. L. (2026). Regression and Latent Variable Models for Bounded Outcomes.
## XLVIII Jornadas Nacionales de Estadística 2026. Universidad de Valparaíso, 
## Octuber 20-23, 2026.
##--------------------------------------------------------------------------------


##1. It installs required packages only when they are not already available.

cran_packages <- c(
  "knitr",
  "gamlss",
  "gamlss.dist",
  "gamlss.data",
  "numDeriv"
)

if (dir.exists("Rlib")) {
  .libPaths(c(normalizePath("Rlib"), .libPaths()))
}

install_if_missing <- function(packages) {
  missing_packages <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing_packages) > 0L) {
    install.packages(missing_packages, repos = "https://cloud.r-project.org")
  }
}

install_if_missing(cran_packages)

if (!requireNamespace("RobustDiscrete", quietly = TRUE)) {
  if (!dir.exists("Rlib")) {
    dir.create("Rlib", recursive = TRUE)
  }

  install.packages(
    "https://github.com/cbayesr/RobustDiscrete_installer/raw/main/RobustDiscrete_0.2.0.tar.gz",
    repos = NULL,
    type = "source",
    lib = "Rlib"
  )

  .libPaths(c(normalizePath("Rlib"), .libPaths()))
}


## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  message = FALSE,
  warning = FALSE
)




## ----load-gamlss--------------------------------------------------------------
library(gamlss)
library(gamlss.dist)


## ----diagnostic-helpers-------------------------------------------------------
compare_models <- function(fit_bi, fit_bb, n) {
  tab <- data.frame(
    model = c("Binomial", "Beta-binomial"),
    AIC = c(AIC(fit_bi), AIC(fit_bb)),
    BIC = c(GAIC(fit_bi, k = log(n)), GAIC(fit_bb, k = log(n)))
  )
  tab$AIC <- round(tab$AIC, 2)
  tab$BIC <- round(tab$BIC, 2)
  tab
}

residual_summary <- function(fit) {
  r <- residuals(fit)
  c(
    mean = mean(r),
    sd = sd(r),
    min = min(r),
    median = median(r),
    max = max(r)
  )
}

plot_residual_diagnostics <- function(fit_bi, fit_bb, main_label) {
  old_par <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  on.exit(par(old_par))

  plot(fitted(fit_bi), residuals(fit_bi),
       xlab = "Fitted values", ylab = "Normalized quantile residuals",
       main = paste(main_label, "- BI"))
  abline(h = 0, col = 2, lty = 2)
  qqnorm(residuals(fit_bi), main = "Normal Q-Q: BI")
  qqline(residuals(fit_bi), col = 2)

  plot(fitted(fit_bb), residuals(fit_bb),
       xlab = "Fitted values", ylab = "Normalized quantile residuals",
       main = paste(main_label, "- BB"))
  abline(h = 0, col = 2, lty = 2)
  qqnorm(residuals(fit_bb), main = "Normal Q-Q: BB")
  qqline(residuals(fit_bb), col = 2)
}

compare_three_models <- function(fit_bi, fit_bb, fit_b2b, n) {
  tab <- data.frame(
    model = c("Binomial", "Beta-binomial", "B2B"),
    AIC = c(AIC(fit_bi), AIC(fit_bb), AIC(fit_b2b)),
    BIC = c(
      GAIC(fit_bi, k = log(n)),
      GAIC(fit_bb, k = log(n)),
      GAIC(fit_b2b, k = log(n))
    )
  )
  tab$AIC <- round(tab$AIC, 2)
  tab$BIC <- round(tab$BIC, 2)
  tab
}

plot_residual_diagnostics_three <- function(fit_bi, fit_bb, fit_b2b, main_label) {
  old_par <- par(mfrow = c(3, 2), mar = c(4, 4, 3, 1))
  on.exit(par(old_par))

  fits <- list(BI = fit_bi, BB = fit_bb, B2B = fit_b2b)
  for (model_name in names(fits)) {
    fit <- fits[[model_name]]
    plot(fitted(fit), residuals(fit),
         xlab = "Fitted values", ylab = "Normalized quantile residuals",
         main = paste(main_label, "-", model_name))
    abline(h = 0, col = 2, lty = 2)
    qqnorm(residuals(fit), main = paste("Normal Q-Q:", model_name))
    qqline(residuals(fit), col = 2)
  }
}


## 3. ----simulated-binomial-regression--------------------------------------------
set.seed(2026)
n <- 400
m <- 20
x <- runif(n, -1.5, 1.5)

eta <- -0.4 + 1.1 * x
mu <- plogis(eta)
y <- rBI(n, bd = m, mu = mu)

dat_bin <- data.frame(
  y = y,
  failures = m - y,
  x = x
)

fit_bin_BI <- gamlss(
  cbind(y, failures) ~ x,
  family = BI,
  data = dat_bin,
  trace = FALSE
)

fit_bin_BB <- gamlss(
  cbind(y, failures) ~ x,
  sigma.formula = ~ 1,
  family = BB,
  data = dat_bin,
  trace = FALSE
)

fit_bin_BI
summary(fit_bin_BI)

fit_bin_BB
summary(fit_bin_BB)

coef(fit_bin_BI, what = "mu")
coef(fit_bin_BB, what = "mu")

sigma_hat_bin <- exp(coef(fit_bin_BB, what = "sigma"))
delta_hat_bin <- sigma_hat_bin / (1 + sigma_hat_bin)
delta_hat_bin

compare_models(fit_bin_BI, fit_bin_BB, n)

deviance(fit_bin_BI)
deviance(fit_bin_BB)

round(rbind(
  Binomial = residual_summary(fit_bin_BI),
  `Beta-binomial` = residual_summary(fit_bin_BB)
), 3)


## ----binomial-data-residual-diagnostics, fig.width=8, fig.height=6, fig.cap="Residual diagnostics for the binomial and beta-binomial fits to data simulated from the binomial regression model."----
pdf("Figure2_binomial_residuals.pdf", width = 8, height = 6)
plot_residual_diagnostics(fit_bin_BI, fit_bin_BB, "Binomial data")
dev.off()


## 4. ----simulated-beta-binomial-regression---------------------------------------
set.seed(2027)
n <- 400
m <- 20
x <- runif(n, -1.5, 1.5)

eta <- -0.4 + 1.1 * x
mu <- plogis(eta)
delta <- 0.30
sigma <- delta / (1 - delta)

y <- rBB(n, bd = m, mu = mu, sigma = sigma)

dat_bb <- data.frame(
  y = y,
  failures = m - y,
  x = x
)

fit_bb_BI <- gamlss(
  cbind(y, failures) ~ x,
  family = BI,
  data = dat_bb,
  trace = FALSE
)

fit_bb_BB <- gamlss(
  cbind(y, failures) ~ x,
  sigma.formula = ~ 1,
  family = BB,
  data = dat_bb,
  trace = FALSE
)

fit_bb_BI
summary(fit_bb_BI)

fit_bb_BB
summary(fit_bb_BB)

coef(fit_bb_BI, what = "mu")
coef(fit_bb_BB, what = "mu")

sigma_hat <- exp(coef(fit_bb_BB, what = "sigma"))
delta_hat <- sigma_hat / (1 + sigma_hat)
delta_hat

compare_models(fit_bb_BI, fit_bb_BB, n)

deviance(fit_bb_BI)
deviance(fit_bb_BB)

round(rbind(
  Binomial = residual_summary(fit_bb_BI),
  `Beta-binomial` = residual_summary(fit_bb_BB)
), 3)


## ----beta-binomial-data-residual-diagnostics, fig.width=8, fig.height=6, fig.cap="Residual diagnostics for the binomial and beta-binomial fits to data simulated from the beta-binomial regression model."----
pdf("Figure3_bb_residuals.pdf", width = 8, height = 6)
plot_residual_diagnostics(fit_bb_BI, fit_bb_BB, "Beta-binomial data")
dev.off()




## 6.----simulated-b2b-regression-------------------------------------------------
if (dir.exists("Rlib")) {
  .libPaths(c("Rlib", .libPaths()))
}

library(numDeriv)
library(RobustDiscrete)

rb2b <- function(n, bd, mu, delta, nu) {
  shape1 <- delta * (1 - nu) / nu
  shape2 <- (1 - delta) * (1 - nu) / nu
  w <- rbeta(n, shape1 = shape1, shape2 = shape2)
  sigma_w <- w / (1 - w)
  rBB(n, bd = bd, mu = mu, sigma = sigma_w)
}

set.seed(2028)
n <- 200
m <- 20
x <- runif(n, -1.5, 1.5)

eta <- -0.4 + 1.1 * x
mu <- plogis(eta)
delta <- 0.30
nu <- 0.60
y <- rb2b(n, bd = m, mu = mu, delta = delta, nu = nu)

dat_b2b <- data.frame(
  y = y,
  failures = m - y,
  x = x
)

control_b2b <- gamlss.control(n.cyc = 100, trace = FALSE)

fit_b2b_BI <- gamlss(
  cbind(y, failures) ~ x,
  family = BI,
  data = dat_b2b,
  control = control_b2b
)

fit_b2b_BB <- gamlss(
  cbind(y, failures) ~ x,
  sigma.formula = ~ 1,
  family = BB,
  data = dat_b2b,
  control = control_b2b
)

fit_b2b_B2B <- gamlss(
  cbind(y, failures) ~ x,
  sigma.formula = ~ 1,
  nu.formula = ~ 1,
  family = B2B,
  data = dat_b2b,
  mu.start = fit_b2b_BB$mu.fv,
  sigma.start = fit_b2b_BB$sigma.fv / (1 + fit_b2b_BB$sigma.fv),
  nu.start = rep(0.5, n),
  control = control_b2b
)

fit_b2b_BI
summary(fit_b2b_BI)

fit_b2b_BB
summary(fit_b2b_BB)

fit_b2b_B2B
summary(fit_b2b_B2B)

coef(fit_b2b_BI, what = "mu")
coef(fit_b2b_BB, what = "mu")
coef(fit_b2b_B2B, what = "mu")

delta_hat_b2b <- plogis(coef(fit_b2b_B2B, what = "sigma"))
nu_hat_b2b <- plogis(coef(fit_b2b_B2B, what = "nu"))
c(delta_hat = delta_hat_b2b, nu_hat = nu_hat_b2b)

compare_three_models(fit_b2b_BI, fit_b2b_BB, fit_b2b_B2B, n)

deviance(fit_b2b_BI)
deviance(fit_b2b_BB)
deviance(fit_b2b_B2B)

round(rbind(
  Binomial = residual_summary(fit_b2b_BI),
  `Beta-binomial` = residual_summary(fit_b2b_BB),
  B2B = residual_summary(fit_b2b_B2B)
), 3)


## ----b2b-data-residual-diagnostics, fig.width=8, fig.height=9, fig.cap="Residual diagnostics for the binomial, beta-binomial, and B2B fits to data simulated from the B2B regression model."----
pdf("Figure4_b2b_residuals.pdf", width = 8, height = 9)
plot_residual_diagnostics_three(fit_b2b_BI, fit_b2b_BB, fit_b2b_B2B, "B2B data")
dev.off()

pdf("Figure4a_b2b_residuals_bi_bb.pdf", width = 8, height = 6)
plot_residual_diagnostics(fit_b2b_BI, fit_b2b_BB, "B2B data")
dev.off()

pdf("Figure4b_b2b_residuals_b2b_only.pdf", width = 8, height = 3.5)
old_par <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
plot(fitted(fit_b2b_B2B), residuals(fit_b2b_B2B),
     xlab = "Fitted values", ylab = "Normalized quantile residuals",
     main = "B2B data - B2B")
abline(h = 0, col = 2, lty = 2)
qqnorm(residuals(fit_b2b_B2B), main = "Normal Q-Q: B2B")
qqline(residuals(fit_b2b_B2B), col = 2)
par(old_par)
dev.off()


## 7. ----repository-b2b-example-code----------------------------------
 if (dir.exists("Rlib")) {
   .libPaths(c("Rlib", .libPaths()))
 }
 
# library(gamlss)
 library(numDeriv)
 library(RobustDiscrete)
# 
# data(aep, package = "gamlss.data")
# 
 M2 <- gamlss(
   y ~ ward + loglos + year,
   sigma.formula = ~ year + ward,
   family = BB,
   data = aep
 )
 
 M5 <- gamlss(
   y ~ ward + loglos + year,
   sigma.formula = ~ year + ward,
   nu.formula = ~ year,
   family = B2B,
   data = aep,
   mu.start = M2$mu.fv,
   sigma.start = M2$sigma.fv / (1 + M2$sigma.fv)
 )
 
 summary(M5)
 AIC(M2, M5)
 GAIC(M2, M5, k = log(nrow(aep)))
