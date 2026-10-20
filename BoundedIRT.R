##-------------------------------------------------------------------------------
## This script implement the Chapter Item Response model of
## the Minicourse
## Bazán, J. L. (2026). Regression and Latent Variable Models for Bounded Outcomes.
## XLVIII Jornadas Nacionales de Estadística 2026. Universidad de Valparaíso, 
## Octuber 20-23, 2026.
##--------------------------------------------------------------------------------




## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  message = FALSE,
  warning = FALSE,
  fig.align = "center"
)

local_lib <- normalizePath("Rlib", mustWork = FALSE)
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

boundirt_available <- requireNamespace("BoundIRT", quietly = TRUE)
if (!boundirt_available) {
  stop(
    "The BoundIRT package is required. Install it with ",
    "install.packages('BoundIRT', repos = 'https://cloud.r-project.org')."
  )
}
set.seed(2022)



## ----custom-priors, eval=FALSE------------------------------------------------
# fit <- BoundIRT::fitBIRT(
#   z = response_matrix,
#   model = "beta",
#   inflated = "auto",
#   prior = list(sA = 2, sB = 5, sP = 3, sG = 5)
# )


## ----detailed-fit-call, eval=FALSE--------------------------------------------
# fit <- BoundIRT::fitBIRT(
#   z = response_matrix,
#   model = "beta",
#   inflated = "auto",
#   iter = 4000,
#   warmup = 2000,
#   nchains = 4,
#   prior = list(sA = 2, sB = 5, sP = 3, sG = 5),
#   seed = 2022,
#   control = list(adapt_delta = 0.95, max_treedepth = 12)
# )


## ----posterior-predictive-check, eval=FALSE-----------------------------------
# y_rep <- BoundIRT::repBIRT(fit, nrep = 500)
# 
# # Example checks for the first item
# observed_zero_rate <- mean(response_matrix[, 1] == 0, na.rm = TRUE)
# replicated_zero_rate <- apply(y_rep[, , 1] == 0, 1, mean, na.rm = TRUE)
# 
# observed_one_rate <- mean(response_matrix[, 1] == 1, na.rm = TRUE)
# replicated_one_rate <- apply(y_rep[, , 1] == 1, 1, mean, na.rm = TRUE)


## ----predictive-comparison-code, eval=FALSE-----------------------------------
# log_lik <- BoundIRT::loglikBIRT(
#   fit,
#   marginal = TRUE,
#   nquad = 15
# )
# 
# loo_result <- loo::loo(log_lik)
# waic_result <- loo::waic(log_lik)


## ----bridge-comparison-code, eval=FALSE---------------------------------------
# bridge_beta <- BoundIRT::bridgeBIRT(
#   beta_fit,
#   niter = 10000,
#   tol = 1e-7
# )
# 
# bridge_simplex <- BoundIRT::bridgeBIRT(simplex_fit)
# 
# log_bf_beta_simplex <-
#   bridge_beta$log_ml - bridge_simplex$log_ml


## ----abasement-install-boundirt, eval=FALSE-----------------------------------
# install.packages("BoundIRT", repos = "https://cloud.r-project.org")
 library(BoundIRT)
#packageVersion("BoundIRT")




## ----abasement-item-analysis, eval=boundirt_available, include=FALSE----------
data("Abasement", package = "BoundIRT")
abasement_complete <- Abasement[stats::complete.cases(Abasement), , drop = FALSE]
abasement_n_items <- ncol(abasement_complete)
abasement_alpha <- abasement_n_items / (abasement_n_items - 1) *
  (1 - sum(apply(abasement_complete, 2, stats::var)) /
     stats::var(rowSums(abasement_complete)))

#abasement_boundaries <- boundary_summary(Abasement)
abasement_total <- rowSums(abasement_complete)
abasement_descriptives <- data.frame(
  item = colnames(abasement_complete),
  Min = vapply(abasement_complete, min, numeric(1)),
  Max = vapply(abasement_complete, max, numeric(1)),
  Mean = vapply(abasement_complete, mean, numeric(1)),
  SD = vapply(abasement_complete, stats::sd, numeric(1)),
  corrected_item_total = vapply(
    abasement_complete,
    function(item) stats::cor(item, abasement_total - item),
    numeric(1)
  ),
  row.names = NULL
)

abasement_item_scores <- data.frame(
  item = rep(colnames(Abasement), each = nrow(Abasement)),
  response = as.vector(as.matrix(Abasement))
)


## ----abasement-boundary-table, eval=boundirt_available, echo=FALSE------------
knitr::kable(
  abasement_boundaries,
  digits = 3,
  row.names = FALSE,
  caption = "Boundary summary of the original Abasement item scores."
)


## ----abasement-descriptives, eval=boundirt_available, echo=FALSE--------------
knitr::kable(
  abasement_descriptives,
  digits = 3,
  row.names = FALSE,
  col.names = c("Item", "Min", "Max", "Mean", "SD", "Corrected item-total r"),
  caption = "Item descriptive statistics and corrected item-total correlations."
)


## ----abasement-item-histograms, eval=boundirt_available, echo=FALSE, fig.width=7.3, fig.height=7.2, fig.cap="Observed Abasement item-score distributions before boundary transformation. Each panel contains 244 responses; zero responses are included."----
ggplot2::ggplot(
  abasement_item_scores,
  ggplot2::aes(x = response)
) +
  ggplot2::geom_histogram(
    binwidth = 0.05,
    boundary = 0,
    closed = "left",
    colour = "white",
    fill = "#007C83"
  ) +
  ggplot2::facet_wrap(~ item, ncol = 3) +
  ggplot2::coord_cartesian(xlim = c(0, 1)) +
  ggplot2::scale_x_continuous(breaks = c(0, 0.5, 1)) +
  ggplot2::labs(x = "Original item score", y = "Number of respondents") +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "grey92", colour = "grey75")
  )


## ----abasement-data, eval=boundirt_available, echo=FALSE----------------------
library(BoundIRT)
data("Abasement", package = "BoundIRT")

abasement_fit_control <- list(
  iter = 2000,
  warmup = 1000,
  nchains = 1,
  silent = TRUE,
  seed = 1310,
  control = list(adapt_delta = 0.95)
)

fit_abasement_model <- function(data, model, inflated) {
  set.seed(abasement_fit_control$seed)
  do.call(
    fitBIRT,
    c(
      list(z = data, model = model, inflated = inflated),
      abasement_fit_control
    )
  )
}

out_beta <- fit_abasement_model(Abasement, model = "beta", inflated = "zero")
out_jsb <- fit_abasement_model(Abasement, model = "samejima", inflated = "zero")
out_simplex <- fit_abasement_model(Abasement, model = "simplex", inflated = "zero")

n <- nrow(Abasement)
Abasement_transformed <- (Abasement * (n - 1) + 0.5) / n
transformed_fits <- setNames(
  lapply(
    c("beta", "samejima", "simplex"),
    function(model) fit_abasement_model(
      Abasement_transformed,
      model = model,
      inflated = "none"
    )
  ),
  c("beta", "samejima", "simplex")
)

loglik_noninflated <- function(object, marginal = TRUE, nquad = 15,
                               silent = TRUE) {
  nit <- ncol(object$z)
  zero_mask <- one_mask <- rep(0, nit)
  zero_mask_matrix <- one_mask_matrix <- diag(0, nit)
  eval(body(BoundIRT::loglikBIRT))
}

format_item_parameters <- function(object) {
  estimates <- coef(object, sd = TRUE)
  data.frame(
    item = rownames(estimates),
    a = estimates[, 1],
    SD_a = estimates[, 2],
    b = estimates[, 3],
    SD_b = estimates[, 4],
    phi = estimates[, 5],
    SD_phi = estimates[, 6],
    gamma0 = estimates[, 7],
    SD_gamma0 = estimates[, 8],
    row.names = NULL
  )
}


## ----refit-abasement, eval=FALSE, echo=FALSE----------------------------------
# set.seed(1310)
# abasement_beta <- fitBIRT(
#   Abasement,
#   model = "beta"
# )
# print(abasement_beta)
# coef(abasement_beta, sd = FALSE)
# summary(abasement_beta, c("alpha", "beta"))


## ----abasement-predictive-fit, eval=boundirt_available, echo=FALSE------------
if (requireNamespace("loo", quietly = TRUE)) {
  abasement_log_lik_beta <- BoundIRT::loglikBIRT(out_beta)
  abasement_log_lik_jsb <- BoundIRT::loglikBIRT(out_jsb)
  abasement_log_lik_simplex <- BoundIRT::loglikBIRT(out_simplex)
  transformed_log_lik <- lapply(transformed_fits, loglik_noninflated)

  abasement_loo_beta <- loo::loo(abasement_log_lik_beta)
  abasement_waic_beta <- loo::waic(abasement_log_lik_beta)
  abasement_loo_jsb <- loo::loo(abasement_log_lik_jsb)
  abasement_waic_jsb <- loo::waic(abasement_log_lik_jsb)
  abasement_loo_simplex <- loo::loo(abasement_log_lik_simplex)
  abasement_waic_simplex <- loo::waic(abasement_log_lik_simplex)
  transformed_loo <- lapply(transformed_log_lik, loo::loo)
  transformed_waic <- lapply(transformed_log_lik, loo::waic)
  abasement_loo_comparison <- loo::loo_compare(
    list(Johnson_SB = abasement_loo_jsb, Beta = abasement_loo_beta,
         Simplex = abasement_loo_simplex)
  )

  abasement_fit_indices <- data.frame(
    data = rep(c("Original; zero-inflated", "Transformed; non-inflated"),
               each = 3),
    model = rep(c("Beta", "Johnson S_B", "Simplex"), 2),
    LOOIC = c(
      abasement_loo_beta$estimates["looic", "Estimate"],
      abasement_loo_jsb$estimates["looic", "Estimate"],
      abasement_loo_simplex$estimates["looic", "Estimate"],
      vapply(transformed_loo, function(x) x$estimates["looic", "Estimate"],
             numeric(1))
    ),
    LOOIC_SE = c(
      abasement_loo_beta$estimates["looic", "SE"],
      abasement_loo_jsb$estimates["looic", "SE"],
      abasement_loo_simplex$estimates["looic", "SE"],
      vapply(transformed_loo, function(x) x$estimates["looic", "SE"],
             numeric(1))
    ),
    WAIC = c(
      abasement_waic_beta$estimates["waic", "Estimate"],
      abasement_waic_jsb$estimates["waic", "Estimate"],
      abasement_waic_simplex$estimates["waic", "Estimate"],
      vapply(transformed_waic, function(x) x$estimates["waic", "Estimate"],
             numeric(1))
    ),
    WAIC_SE = c(
      abasement_waic_beta$estimates["waic", "SE"],
      abasement_waic_jsb$estimates["waic", "SE"],
      abasement_waic_simplex$estimates["waic", "SE"],
      vapply(transformed_waic, function(x) x$estimates["waic", "SE"],
             numeric(1))
    )
  )
  best_model <- c("beta", "jsb", "simplex")[
    which.min(abasement_fit_indices$LOOIC[1:3])
  ]
  best_fit <- list(beta = out_beta, jsb = out_jsb,
                   simplex = out_simplex)[[best_model]]
  knitr::kable(
    abasement_fit_indices,
    digits = 2,
    row.names = FALSE,
    caption = "Person-marginal predictive criteria. Compare values only within each data group."
  )
} else {
  best_model <- "jsb"
  best_fit <- out_jsb
}


## ----abasement-best-parameters, eval=boundirt_available, echo=FALSE-----------
knitr::kable(
  format_item_parameters(best_fit), digits = 3, row.names = FALSE,
  caption = "Item-parameter posterior means and standard deviations for the selected Johnson S_B model."
)


## ----abasement-best-icc, eval=boundirt_available, echo=FALSE, fig.width=7.2, fig.height=4.8, fig.cap="Item characteristic curves for selected Abasement items under the zero-inflated Johnson $S_B$ IRT model."----
invisible(plot(
  best_fit, what = "icc", items = c(1, 3, 6, 9),
  col = c("#007C83", "#C44E52", "#222222", "#A36A00"),
  lwd = 2
))
legend(
  "topleft", legend = colnames(Abasement)[c(1, 3, 6, 9)],
  col = c("#007C83", "#C44E52", "#222222", "#A36A00"),
  lty = 1, lwd = 2, bty = "n", cex = 0.85
)


## ----abasement-jsb-item-map, eval=boundirt_available, echo=FALSE, fig.width=7.2, fig.height=5.4, fig.cap="Posterior mean discrimination and easiness parameters under the selected Johnson $S_B$ IRT model. Dashed lines mark the respective item-parameter means."----
jsb_item_parameters <- format_item_parameters(out_jsb)
jsb_item_parameters$label <- colnames(Abasement)

ggplot2::ggplot(
  jsb_item_parameters,
  ggplot2::aes(x = a, y = b, label = label)
) +
  ggplot2::geom_vline(
    xintercept = mean(jsb_item_parameters$a),
    colour = "grey55",
    linetype = "dashed",
    linewidth = 0.45
  ) +
  ggplot2::geom_hline(
    yintercept = mean(jsb_item_parameters$b),
    colour = "grey55",
    linetype = "dashed",
    linewidth = 0.45
  ) +
  ggplot2::geom_point(colour = "#007C83", size = 3) +
  ggrepel::geom_text_repel(
    colour = "black",
    size = 3.2,
    box.padding = 0.35,
    point.padding = 0.25,
    min.segment.length = 0,
    seed = 1310,
    max.overlaps = Inf
  ) +
  ggplot2::labs(
    x = expression("Discrimination " * a[j]),
    y = expression("Easiness " * b[j])
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank())


## ----abasement-jsb-theta-estimates, eval=boundirt_available, echo=FALSE-------
jsb_theta_estimates <- data.frame(
  person = seq_len(nrow(best_fit$theta)),
  observed_score = rowSums(Abasement),
  posterior_mean = best_fit$theta[, "theta_p"],
  posterior_sd = best_fit$theta[, "sd"]
)
write.csv(
  jsb_theta_estimates,
  "section8_theta_estimates.csv",
  row.names = FALSE
)


## ----abasement-jsb-theta-summary, eval=boundirt_available, echo=FALSE---------
theta_posterior_means <- jsb_theta_estimates$posterior_mean
theta_descriptive_statistics <- data.frame(
  statistic = c(
    "N", "Mean", "Standard deviation", "Minimum", "2.5%", "First quartile",
    "Median", "Third quartile", "97.5%", "Maximum"
  ),
  value = c(
    length(theta_posterior_means),
    mean(theta_posterior_means),
    sd(theta_posterior_means),
    min(theta_posterior_means),
    unname(quantile(theta_posterior_means, 0.025)),
    unname(quantile(theta_posterior_means, 0.25)),
    median(theta_posterior_means),
    unname(quantile(theta_posterior_means, 0.75)),
    unname(quantile(theta_posterior_means, 0.975)),
    max(theta_posterior_means)
  )
)
knitr::kable(
  theta_descriptive_statistics,
  digits = 3,
  row.names = FALSE,
  caption = "Descriptive statistics of the person-parameter posterior means."
)


## ----abasement-jsb-theta-histogram, eval=boundirt_available, echo=FALSE, fig.width=6.7, fig.height=4.6, fig.cap="Distribution of posterior mean person parameters under the selected Johnson $S_B$ IRT model."----
ggplot2::ggplot(
  jsb_theta_estimates,
  ggplot2::aes(x = posterior_mean)
) +
  ggplot2::geom_histogram(
    bins = 20,
    boundary = 0,
    colour = "black",
    fill = "#007C83",
    linewidth = 0.35
  ) +
  ggplot2::geom_vline(
    xintercept = mean(theta_posterior_means),
    colour = "black",
    linetype = "dashed",
    linewidth = 0.65
  ) +
  ggplot2::labs(
    x = expression("Posterior mean of " * theta[i]),
    y = "Number of persons"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank())


## ----abasement-score-theta, eval=boundirt_available, include=FALSE------------
score_theta_correlation <- stats::cor(
  jsb_theta_estimates$observed_score,
  jsb_theta_estimates$posterior_mean
)

# Group scores after rounding only to remove floating-point summation noise.
score_groups <- split(
  jsb_theta_estimates$person,
  round(jsb_theta_estimates$observed_score, 8)
)
score_groups <- score_groups[lengths(score_groups) >= 2]
theta_spread <- vapply(
  score_groups,
  function(ids) diff(range(jsb_theta_estimates$posterior_mean[ids])),
  numeric(1)
)
selected_score_groups <- names(sort(theta_spread, decreasing = TRUE))[1:2]
profile_people <- do.call(rbind, lapply(selected_score_groups, function(key) {
  ids <- score_groups[[key]]
  ids <- ids[order(jsb_theta_estimates$posterior_mean[ids])]
  selected <- c(ids[1], ids[length(ids)])
  data.frame(
    score_group = sprintf("Total score = %.3f", as.numeric(key)),
    person = selected,
    observed_score = jsb_theta_estimates$observed_score[selected],
    posterior_mean = jsb_theta_estimates$posterior_mean[selected],
    posterior_sd = jsb_theta_estimates$posterior_sd[selected],
    theta_group = c("Lower theta", "Higher theta")
  )
}))

profile_responses <- do.call(rbind, lapply(seq_len(nrow(profile_people)), function(k) {
  id <- profile_people$person[k]
  data.frame(
    score_group = profile_people$score_group[k],
    theta_group = profile_people$theta_group[k],
    item = factor(colnames(Abasement), levels = colnames(Abasement)),
    response = as.numeric(Abasement[id, ])
  )
}))


## ----abasement-score-theta-plot, eval=boundirt_available, echo=FALSE, fig.width=6.7, fig.height=4.6, fig.cap="Observed ten-item total score versus posterior mean theta under the selected Johnson $S_B$ IRT model."----
ggplot2::ggplot(
  jsb_theta_estimates,
  ggplot2::aes(x = observed_score, y = posterior_mean)
) +
  ggplot2::geom_point(colour = "#007C83", alpha = 0.65, size = 1.8) +
  ggplot2::geom_smooth(
    method = "lm", se = FALSE, colour = "#C44E52", linewidth = 0.8
  ) +
  ggplot2::labs(
    x = "Observed total score",
    y = expression("Posterior mean of " * theta[i])
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank())


## ----abasement-tied-score-people, eval=boundirt_available, echo=FALSE---------
knitr::kable(
  profile_people[, c("score_group", "person", "observed_score",
                     "posterior_mean", "posterior_sd")],
  digits = 3,
  row.names = FALSE,
  col.names = c("Score group", "Person", "Total score", "Theta mean", "Theta SD"),
  caption = "Respondents with matching total scores but different theta estimates."
)


## ----abasement-tied-score-profiles, eval=boundirt_available, echo=FALSE, fig.width=7.3, fig.height=5.8, fig.cap="Item-response profiles for two pairs of respondents with equal total scores but different posterior mean theta values."----
ggplot2::ggplot(
  profile_responses,
  ggplot2::aes(
    x = item, y = response,
    colour = theta_group, linetype = theta_group, group = theta_group
  )
) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 2) +
  ggplot2::facet_wrap(~ score_group, ncol = 1) +
  ggplot2::scale_colour_manual(
    values = c("Lower theta" = "#007C83", "Higher theta" = "#C44E52")
  ) +
  ggplot2::scale_linetype_manual(
    values = c("Lower theta" = "solid", "Higher theta" = "dashed")
  ) +
  ggplot2::scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  ggplot2::labs(x = "Item", y = "Original response", colour = NULL, linetype = NULL) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    legend.position = "top",
    panel.grid.minor = ggplot2::element_blank()
  )


## ----session-information------------------------------------------------------
sessionInfo()
