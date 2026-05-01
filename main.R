# install.packages(c("tidyverse", "sandwich", "lmtest", "crossfit", "ggplot2", "ranger"))

library(tidyverse)
library(crossfit)
library(ggplot2)

if (!requireNamespace("ranger", quietly = TRUE)) {
  stop("Please install the 'ranger' package: install.packages('ranger')", call. = FALSE)
}

# 1) Load the What If / NHEFS data
data_path <- if (file.exists("Data/nhefs.csv")) {
  "Data/nhefs.csv"
} else {
  "data/nhefs.csv"
}

nhefs <- read_csv(data_path, show_col_types = FALSE)

# 2) Quick look
#glimpse(nhefs)
#summary(nhefs)

# 3) Keep complete cases for a first pass
dat <- nhefs %>%
  filter(
    !is.na(qsmk),
    !is.na(wt82_71),
    !is.na(sex),
    !is.na(race),
    !is.na(age),
    !is.na(education),
    !is.na(smokeintensity),
    !is.na(smokeyrs),
    !is.na(exercise),
    !is.na(active),
    !is.na(wt71)
  )

# Treatment: qsmk = quit smoking
# Outcome: wt82_71 = weight change from 1971 to 1982

# 4) Naive difference in means
dat %>%
  group_by(qsmk) %>%
  summarise(mean_outcome = mean(wt82_71), n = n())

naive_ate <- with(dat, mean(wt82_71[qsmk == 1]) - mean(wt82_71[qsmk == 0]))
naive_ate

naive_res <- with(
  dat,
  tibble(
    estimate = naive_ate,
    variance = var(wt82_71[qsmk == 1]) / sum(qsmk == 1) +
      var(wt82_71[qsmk == 0]) / sum(qsmk == 0)
  )
) %>%
  mutate(
    std_error = sqrt(variance),
    ci_lower = estimate - 1.96 * std_error,
    ci_upper = estimate + 1.96 * std_error
  )

# 5) Outcome regression (g-computation style)
outcome_model <- lm(
  wt82_71 ~ qsmk + sex + race + age + I(age^2) + education +
    smokeintensity + I(smokeintensity^2) + smokeyrs + I(smokeyrs^2) +
    exercise + active + wt71 + I(wt71^2),
  data = dat
)

dat1 <- dat %>% mutate(qsmk = 1)
dat0 <- dat %>% mutate(qsmk = 0)

mu1_hat <- predict(outcome_model, newdata = dat1)
mu0_hat <- predict(outcome_model, newdata = dat0)

gcomp_ate <- mean(mu1_hat - mu0_hat)
gcomp_ate
gcomp_se <- summary(outcome_model)$coefficients["qsmk", "Std. Error"]
gcomp_res <- tibble(
  estimate = gcomp_ate,
  variance = gcomp_se^2,
  std_error = gcomp_se,
  ci_lower = gcomp_ate - 1.96 * gcomp_se,
  ci_upper = gcomp_ate + 1.96 * gcomp_se
)

# 6) Propensity score model
ps_model <- glm(
  qsmk ~ sex + race + age + I(age^2) + education +
    smokeintensity + I(smokeintensity^2) + smokeyrs + I(smokeyrs^2) +
    exercise + active + wt71 + I(wt71^2),
  family = binomial(),
  data = dat
)

dat$ps <- predict(ps_model, type = "response")

# Optional trimming for stability
dat$ps <- pmin(pmax(dat$ps, 0.01), 0.99)

# 7) AIPW / one-step estimator
dat$mu1 <- predict(outcome_model, newdata = dat1)
dat$mu0 <- predict(outcome_model, newdata = dat0)

dat$phi <- with(dat,
  qsmk / ps * (wt82_71 - mu1) -
    (1 - qsmk) / (1 - ps) * (wt82_71 - mu0) +
    mu1 - mu0
)

aipw_ate <- mean(dat$phi)
aipw_se  <- sd(dat$phi) / sqrt(nrow(dat))
aipw_var <- aipw_se^2
aipw_ci  <- c(aipw_ate - 1.96 * aipw_se, aipw_ate + 1.96 * aipw_se)

# 8) Cross-fitted AIPW using the crossfit package
nuis_ps <- create_nuisance(
  fit = function(data, ...) {
    glm(
      qsmk ~ sex + race + age + I(age^2) + education +
        smokeintensity + I(smokeintensity^2) + smokeyrs + I(smokeyrs^2) +
        exercise + active + wt71 + I(wt71^2),
      family = binomial(),
      data = data
    )
  },
  predict = function(model, data, ...) {
    p <- as.numeric(predict(model, newdata = data, type = "response"))
    pmin(pmax(p, 0.01), 0.99)
  }
)

nuis_mu1 <- create_nuisance(
  fit = function(data, ...) {
    lm(
      wt82_71 ~ qsmk + sex + race + age + I(age^2) + education +
        smokeintensity + I(smokeintensity^2) + smokeyrs + I(smokeyrs^2) +
        exercise + active + wt71 + I(wt71^2),
      data = data
    )
  },
  predict = function(model, data, ...) {
    data1 <- data
    data1$qsmk <- 1
    as.numeric(predict(model, newdata = data1))
  }
)

nuis_mu0 <- create_nuisance(
  fit = function(data, ...) {
    lm(
      wt82_71 ~ qsmk + sex + race + age + I(age^2) + education +
        smokeintensity + I(smokeintensity^2) + smokeyrs + I(smokeyrs^2) +
        exercise + active + wt71 + I(wt71^2),
      data = data
    )
  },
  predict = function(model, data, ...) {
    data0 <- data
    data0$qsmk <- 0
    as.numeric(predict(model, newdata = data0))
  }
)

target_aipw <- function(data, ps_hat, mu1_hat, mu0_hat, ...) {
  mean(
    data$qsmk / ps_hat * (data$wt82_71 - mu1_hat) -
      (1 - data$qsmk) / (1 - ps_hat) * (data$wt82_71 - mu0_hat) +
      mu1_hat - mu0_hat
  )
}

method_cf_aipw <- create_method(
  target = target_aipw,
  list_nuisance = list(
    ps_hat = nuis_ps,
    mu1_hat = nuis_mu1,
    mu0_hat = nuis_mu0
  ),
  folds = 5,
  repeats = 1,
  eval_fold = 1,
  mode = "estimate",
  fold_allocation = "independence",
  aggregate_panels = mean_estimate,
  aggregate_repeats = mean_estimate
)

set.seed(20260429)
cf_res <- crossfit(dat, method_cf_aipw)
cf_aipw_ate <- cf_res$estimates[[1]]

estimate_glm_cf_aipw <- function(data, n_folds = 5, seed = 20260429) {
  set.seed(seed)
  fold_id <- sample(rep(seq_len(n_folds), length.out = nrow(data)))
  nuisances <- tibble(
    ps_hat = rep(NA_real_, nrow(data)),
    mu1_hat = rep(NA_real_, nrow(data)),
    mu0_hat = rep(NA_real_, nrow(data))
  )

  for (fold in seq_len(n_folds)) {
    train <- data[fold_id != fold, ]
    holdout <- data[fold_id == fold, ]

    ps_fit <- glm(
      qsmk ~ sex + race + age + I(age^2) + education +
        smokeintensity + I(smokeintensity^2) + smokeyrs + I(smokeyrs^2) +
        exercise + active + wt71 + I(wt71^2),
      family = binomial(),
      data = train
    )

    outcome_fit <- lm(
      wt82_71 ~ qsmk + sex + race + age + I(age^2) + education +
        smokeintensity + I(smokeintensity^2) + smokeyrs + I(smokeyrs^2) +
        exercise + active + wt71 + I(wt71^2),
      data = train
    )

    holdout1 <- holdout
    holdout0 <- holdout
    holdout1$qsmk <- 1
    holdout0$qsmk <- 0

    nuisances[fold_id == fold, ] <- tibble(
      ps_hat = pmin(
        pmax(as.numeric(predict(ps_fit, newdata = holdout, type = "response")), 0.01),
        0.99
      ),
      mu1_hat = as.numeric(predict(outcome_fit, newdata = holdout1)),
      mu0_hat = as.numeric(predict(outcome_fit, newdata = holdout0))
    )
  }

  phi <- with(
    bind_cols(data, nuisances),
    qsmk / ps_hat * (wt82_71 - mu1_hat) -
      (1 - qsmk) / (1 - ps_hat) * (wt82_71 - mu0_hat) +
      mu1_hat - mu0_hat
  )

  ate <- mean(phi)
  se <- sd(phi) / sqrt(nrow(data))

  tibble(
    estimate = ate,
    variance = se^2,
    std_error = se,
    ci_lower = ate - 1.96 * se,
    ci_upper = ate + 1.96 * se
  )
}

cf_aipw_res <- estimate_glm_cf_aipw(dat)
cf_aipw_ate <- cf_aipw_res$estimate

# 9) Cross-fitted AIPW with random forest nuisance models
rf_dat <- dat %>%
  select(
    wt82_71, qsmk, sex, race, age, education,
    smokeintensity, smokeyrs, exercise, active, wt71
  ) %>%
  mutate(
    qsmk_num = qsmk,
    qsmk = factor(qsmk, levels = c(0, 1), labels = c("No quit", "Quit")),
    sex = factor(sex),
    race = factor(race),
    education = factor(education),
    exercise = factor(exercise),
    active = factor(active)
  )

rf_covariates <- c(
  "sex", "race", "age", "education", "smokeintensity",
  "smokeyrs", "exercise", "active", "wt71"
)

rf_ps_formula <- as.formula(paste("qsmk ~", paste(rf_covariates, collapse = " + ")))
rf_outcome_formula <- as.formula(
  paste("wt82_71 ~ qsmk +", paste(rf_covariates, collapse = " + "))
)

clip_ps <- function(ps, lower = 0.01, upper = 0.99) {
  pmin(pmax(ps, lower), upper)
}

fit_rf_nuisance <- function(train_data, seed) {
  ps_fit <- ranger::ranger(
    formula = rf_ps_formula,
    data = train_data,
    probability = TRUE,
    num.trees = 500,
    mtry = 4,
    min.node.size = 5,
    importance = "permutation",
    seed = seed
  )

  outcome_fit <- ranger::ranger(
    formula = rf_outcome_formula,
    data = train_data,
    num.trees = 500,
    mtry = 4,
    min.node.size = 5,
    importance = "permutation",
    seed = seed + 1
  )

  list(ps = ps_fit, outcome = outcome_fit)
}

predict_rf_nuisance <- function(fits, new_data) {
  ps_hat <- predict(fits$ps, data = new_data)$predictions[, "Quit"]

  data1 <- new_data
  data0 <- new_data
  data1$qsmk <- factor("Quit", levels = levels(rf_dat$qsmk))
  data0$qsmk <- factor("No quit", levels = levels(rf_dat$qsmk))

  tibble(
    ps_hat = clip_ps(as.numeric(ps_hat)),
    mu1_hat = as.numeric(predict(fits$outcome, data = data1)$predictions),
    mu0_hat = as.numeric(predict(fits$outcome, data = data0)$predictions)
  )
}

estimate_rf_aipw <- function(data, nuisances) {
  phi <- with(
    bind_cols(data, nuisances),
    qsmk_num / ps_hat * (wt82_71 - mu1_hat) -
      (1 - qsmk_num) / (1 - ps_hat) * (wt82_71 - mu0_hat) +
      mu1_hat - mu0_hat
  )

  ate <- mean(phi)
  se <- sd(phi) / sqrt(nrow(data))

  tibble(
    estimate = ate,
    variance = se^2,
    std_error = se,
    ci_lower = ate - 1.96 * se,
    ci_upper = ate + 1.96 * se
  )
}

set.seed(20260429)
rf_n_folds <- 5
rf_fold_id <- sample(rep(seq_len(rf_n_folds), length.out = nrow(rf_dat)))
rf_cf_nuisances <- tibble(
  ps_hat = rep(NA_real_, nrow(rf_dat)),
  mu1_hat = rep(NA_real_, nrow(rf_dat)),
  mu0_hat = rep(NA_real_, nrow(rf_dat))
)

for (fold in seq_len(rf_n_folds)) {
  rf_train <- rf_dat[rf_fold_id != fold, ]
  rf_holdout <- rf_dat[rf_fold_id == fold, ]
  rf_fits <- fit_rf_nuisance(rf_train, seed = 20260429 + fold * 10)
  rf_cf_nuisances[rf_fold_id == fold, ] <- predict_rf_nuisance(rf_fits, rf_holdout)
}

rf_aipw_res <- estimate_rf_aipw(rf_dat, rf_cf_nuisances)
rf_cf_aipw_ate <- rf_aipw_res$estimate

results_tbl <- tibble(
  method = c(
    "Naive",
    "G-computation",
    "AIPW",
    "Cross-fitted AIPW",
    "Cross-fitted RF AIPW"
  ),
  estimate = c(naive_ate, gcomp_ate, aipw_ate, cf_aipw_ate, rf_cf_aipw_ate),
  variance = c(
    naive_res$variance,
    gcomp_res$variance,
    aipw_var,
    cf_aipw_res$variance,
    rf_aipw_res$variance
  ),
  std_error = c(
    naive_res$std_error,
    gcomp_res$std_error,
    aipw_se,
    cf_aipw_res$std_error,
    rf_aipw_res$std_error
  ),
  ci_lower = c(
    naive_res$ci_lower,
    gcomp_res$ci_lower,
    aipw_ci[1],
    cf_aipw_res$ci_lower,
    rf_aipw_res$ci_lower
  ),
  ci_upper = c(
    naive_res$ci_upper,
    gcomp_res$ci_upper,
    aipw_ci[2],
    cf_aipw_res$ci_upper,
    rf_aipw_res$ci_upper
  )
)

print(results_tbl)
print(nrow(dat))
