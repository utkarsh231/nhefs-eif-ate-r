# install.packages(c("tidyverse", "sandwich", "lmtest", "crossfit", "ggplot2"))

library(tidyverse)
library(crossfit)
library(ggplot2)

# 1) Load the What If / NHEFS data
nhefs <- read_csv("data/nhefs.csv", show_col_types = FALSE)

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

cf_res <- crossfit(dat, method_cf_aipw)
cf_aipw_ate <- cf_res$estimates[[1]]

results_tbl <- tibble(
  method = c("Naive", "G-computation", "AIPW", "Cross-fitted AIPW"),
  estimate = c(naive_ate, gcomp_ate, aipw_ate, cf_aipw_ate),
  std_error = c(NA_real_, NA_real_, aipw_se, NA_real_),
  ci_lower = c(NA_real_, NA_real_, aipw_ci[1], NA_real_),
  ci_upper = c(NA_real_, NA_real_, aipw_ci[2], NA_real_)
)

print(results_tbl)

ps_overlap_plot <- ggplot(dat, aes(x = ps, fill = factor(qsmk))) +
  geom_density(alpha = 0.4) +
  labs(
    title = "Propensity Score Overlap by Treatment Group",
    x = "Estimated propensity score",
    y = "Density",
    fill = "qsmk"
  )

print(nrow(dat))

print(ps_overlap_plot)

print(naive_ate)
print(gcomp_ate)
print(aipw_ate)
print(aipw_se)
print(aipw_ci)
print(cf_aipw_ate)