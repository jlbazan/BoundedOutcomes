##-------------------------------------------------------------------------------
## This script implements the Chapter Bounded Continuous Regression model of
## the Minicourse
## Bazán, J. L. (2026). Regression and Latent Variable Models for Bounded Outcomes.
## XLVIII Jornadas Nacionales de Estadística 2026. Universidad de Valparaíso,
## October 20-23, 2026.
##--------------------------------------------------------------------------------


## 1. It installs required packages only when they are not already available.

cran_packages <- c(
  "knitr",
  "gamlss",
  "gamlss.dist",
  "ggplot2"
)

install_if_missing <- function(packages) {
  missing_packages <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing_packages) > 0L) {
    install.packages(missing_packages, repos = "https://cloud.r-project.org")
  }
}

install_if_missing(cran_packages)


## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  message = FALSE,
  warning = FALSE
)

#' # Practical Estimation Using GAMLSS
#' 
#' 
## ----packages-----------------------------------------------------------------
# Load the modeling functions, distributions, and plotting tools.
library(gamlss)
library(gamlss.dist)
library(ggplot2)

#' 
## ----generate-families, results='hide'----------------------------------------
# Generate logit-transformed JSBt families from GAMLSS baselines.
gen.Family(family = "GU", type = "logit")
gen.Family(family = "RG", type = "logit")
gen.Family(family = "exGAUS", type = "logit")
gen.Family(family = "SN1", type = "logit")
gen.Family(family = "SN2", type = "logit")
gen.Family(family = "PE", type = "logit")
gen.Family(family = "LO", type = "logit")
gen.Family(family = "NO", type = "logit")
gen.Family(family = "TF", type = "logit")

#' 
## ----simulate-data------------------------------------------------------------
# Use the same covariate values in both simulations.
set.seed(660)

n <- 500
x <- rnorm(n)

# Beta-generated case: x enters the conditional mean.
beta_data <- data.frame(Scenario = "Beta-generated", x = x)
beta_data$eta <- -0.40 + 1.20 * beta_data$x
beta_data$true_location <- plogis(beta_data$eta)
beta_data$sigma <- 0.20
beta_data$nu <- 0
beta_data$y <- rBE(
  n = n,
  mu = beta_data$true_location,
  sigma = beta_data$sigma
)

# Logistic-LSN1-generated case: x enters the real-line location.
lsn1_data <- data.frame(Scenario = "Logistic-LSN1-generated", x = x)
lsn1_data$eta <- -0.40 + 1.20 * lsn1_data$x
lsn1_data$true_location <- plogis(lsn1_data$eta)
lsn1_data$sigma <- 0.45
lsn1_data$nu <- 5
set.seed(650)
lsn1_data$y <- rlogitSN1(
  n = n,
  mu = lsn1_data$eta,
  sigma = lsn1_data$sigma,
  nu = lsn1_data$nu
)

# Store both datasets in a list so the same fitting code can be reused.
simulations <- list(
  Beta = beta_data,
  LSN1 = lsn1_data
)

lapply(simulations, function(dat) summary(dat$y))

#' 
## ----exploratory-plots, fig.width=9, fig.height=4-----------------------------
# Combine the two simulations only for graphical display.
sim_plot_data <- do.call(rbind, simulations)
mean_curve <- sim_plot_data[order(sim_plot_data$Scenario, sim_plot_data$x), ]

p1 <- ggplot(sim_plot_data, aes(x = y)) +
  geom_histogram(bins = 30, color = "black", fill = "grey80") +
  facet_wrap(~ Scenario) +
  theme_classic() +
  labs(x = "Simulated bounded response", y = "Count")

p2 <- ggplot(sim_plot_data, aes(x = x, y = y)) +
  geom_point(alpha = 0.35) +
  geom_line(
    data = mean_curve,
    aes(y = true_location),
    color = "red",
    linewidth = 1
  ) +
  facet_wrap(~ Scenario) +
  theme_classic() +
  labs(x = "Covariate x", y = "Simulated bounded response")

print(p1)
print(p2)

#' 
#' ## Fitting JSBt Bounded Regression Models
#' 
#' For each simulated dataset, all fitted models use the same predictor in
#' the location or mean regression and keep the remaining parameters fixed.
#' For the beta model, `mu` is the conditional mean of the response. For
#' the JSBt models, the linear predictor is assigned to the real-line
#' location parameter of the baseline distribution and is transformed back
#' to the unit interval by the inverse-logit map.
#' 
## ----fit-models---------------------------------------------------------------
ctrl <- gamlss.control(n.cyc = 100, trace = FALSE)

# Fit the beta benchmark and all JSBt candidates to one dataset.
fit_regression_models <- function(dat) {
  list(
    Beta = gamlss(
      y ~ x,
      sigma.formula = ~ 1,
      data = dat,
      family = BE,
      control = ctrl
    ),
    LGU = gamlss(
      y ~ x,
      sigma.formula = ~ 1,
      data = dat,
      family = logitGU,
      control = ctrl
    ),
    LRG = gamlss(
      y ~ x,
      sigma.formula = ~ 1,
      data = dat,
      family = logitRG,
      control = ctrl
    ),
    LLo = gamlss(
      y ~ x,
      sigma.formula = ~ 1,
      data = dat,
      family = logitLO,
      control = ctrl
    ),
    LN = gamlss(
      y ~ x,
      sigma.formula = ~ 1,
      data = dat,
      family = logitNO,
      control = ctrl
    ),
    LEG = gamlss(
      y ~ x,
      sigma.formula = ~ 1,
      nu.formula = ~ 1,
      data = dat,
      family = logitexGAUS,
      control = ctrl
    ),
    LPE = gamlss(
      y ~ x,
      sigma.formula = ~ 1,
      nu.formula = ~ 1,
      data = dat,
      family = logitPE,
      control = ctrl
    ),
    LSN1 = gamlss(
      y ~ x,
      sigma.formula = ~ 1,
      nu.formula = ~ 1,
      data = dat,
      family = logitSN1,
      control = ctrl
    ),
    LSN2 = gamlss(
      y ~ x,
      sigma.formula = ~ 1,
      nu.formula = ~ 1,
      data = dat,
      family = logitSN2,
      control = ctrl
    ),
    LT = gamlss(
      y ~ x,
      sigma.formula = ~ 1,
      nu.formula = ~ 1,
      data = dat,
      family = logitTF,
      control = ctrl
    )
  )
}

# Apply the same model set to the beta and Logistic-LSN1 simulations.
fits <- lapply(simulations, fit_regression_models)

#' 
#' ## Model Comparison
#' 
#' The following table compares the fitted beta and JSBt bounded regression
#' models for both simulated datasets using AIC and BIC. Smaller values
#' indicate a better trade-off between fit and complexity.
#' 
## ----model-comparison---------------------------------------------------------
# Build one AIC/BIC table containing both simulation scenarios.
criteria <- do.call(
  rbind,
  lapply(names(fits), function(scenario_name) {
    model_list <- fits[[scenario_name]]
    data.frame(
      Scenario = scenario_name,
      Model = names(model_list),
      AIC = sapply(model_list, AIC),
      BIC = sapply(model_list, BIC),
      row.names = NULL
    )
  })
)

criteria <- criteria[order(criteria$Scenario, criteria$AIC), ]

# Select the best model in each simulation case by minimum AIC.
best_names <- sapply(names(fits), function(scenario_name) {
  scenario_criteria <- criteria[criteria$Scenario == scenario_name, ]
  scenario_criteria$Model[which.min(scenario_criteria$AIC)]
})

# Keep the selected fitted objects for diagnostics and plots.
best_models <- lapply(names(fits), function(scenario_name) {
  fits[[scenario_name]][[best_names[scenario_name]]]
})
names(best_models) <- names(fits)

criteria

#' 
#' ## Residual Diagnostics
#' 
#' Randomized quantile residuals are standard diagnostics in GAMLSS models.
#' For a continuous response, if the fitted model is adequate, the
#' residuals should be approximately standard normal.
#' 
## ----residual-diagnostics, fig.width=9, fig.height=4--------------------------
# Randomized quantile residuals should be approximately normal.
residual_list <- lapply(best_models, residuals)

par(mfrow = c(length(best_models), 2))
for (scenario_name in names(best_models)) {
  rq <- residual_list[[scenario_name]]
  best_name <- best_names[scenario_name]

  hist(
    rq,
    main = paste("Residuals:", scenario_name, "-", best_name),
    xlab = "Randomized quantile residuals",
    col = "grey80",
    border = "black"
  )
  qqnorm(rq, main = paste("Q-Q plot:", scenario_name, "-", best_name))
  qqline(rq, lwd = 2)
}
par(mfrow = c(1, 1))

data.frame(
  Scenario = names(best_models),
  SelectedModel = best_names[names(best_models)],
  ShapiroPValue = sapply(residual_list, function(rq) shapiro.test(rq)$p.value),
  row.names = NULL
)

#' 
#' ## Fitted Location on the Response Scale
#' 
#' For the beta regression model, the fitted `mu` parameter is already the
#' conditional mean on the unit interval. For a selected JSBt model, the
#' fitted location parameter is estimated on the real-line baseline scale,
#' so the inverse-logit transformation gives a location index on the
#' response scale. The plot is produced for the selected model in each
#' simulation case.
#' 
## ----fitted-location, fig.width=7, fig.height=4-------------------------------
# Put the selected fitted location on the response scale.
fitted_location_data <- do.call(
  rbind,
  lapply(names(best_models), function(scenario_name) {
    dat <- simulations[[scenario_name]]
    best_name <- best_names[scenario_name]
    best_model <- best_models[[scenario_name]]
    dat$SelectedModel <- best_name
    dat$mu_hat <- fitted(best_model, what = "mu")
    dat$location_hat <- if (best_name == "Beta") {
      dat$mu_hat
    } else {
      plogis(dat$mu_hat)
    }
    dat
  })
)

plot_data <- fitted_location_data[
  order(fitted_location_data$Scenario, fitted_location_data$x),
]

ggplot(fitted_location_data, aes(x = x, y = y)) +
  geom_point(alpha = 0.35) +
  geom_line(
    data = plot_data,
    aes(y = true_location),
    linetype = 2,
    linewidth = 1
  ) +
  geom_line(
    data = plot_data,
    aes(y = location_hat),
    color = "red",
    linewidth = 1
  ) +
  facet_wrap(~ Scenario) +
  theme_classic() +
  labs(
    x = "Covariate x",
    y = "Response, true mean, and fitted location",
    title = "Fitted response-scale location by simulation case"
  )

