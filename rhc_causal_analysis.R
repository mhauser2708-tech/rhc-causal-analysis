# Libraries
library(tidyverse)
library(fixest)
library(causalweight)
library(ggplot2)
library(DoubleML)
library(mlr3)
library(mlr3learners)
library(data.table)
library(grf)


# Load and prepare data
df <- read.csv("data/rhc.csv", row.names = 1)

# Treatment and outcome
df$treat <- as.integer(df$swang1 == "RHC")
df$died  <- as.integer(df$death == "Yes")

# Convert Yes/No variables to 0/1
yn_cols <- c(
  "resp","card","neuro","gastr","renal","meta","hema",
  "seps","trauma","ortho","dnr1","dth30"
)

for (col in yn_cols) {
  df[[col]] <- as.integer(df[[col]] == "Yes")
}

# Covariates
covs_cont <- c(
  "age","edu","surv2md1","das2d3pc","aps1","scoma1",
  "meanbp1","wblc1","hrt1","resp1","temp1","pafi1",
  "alb1","hema1","bili1","crea1","sod1","pot1",
  "paco21","ph1","wtkilo1"
)

covs_bin <- c(
  "cardiohx","chfhx","dementhx","psychhx","chrpulhx","renalhx",
  "liverhx","gibledhx","malighx","immunhx","transhx","amihx",
  "resp","card","neuro","gastr","renal","meta","hema","seps",
  "trauma","ortho","dnr1"
)

all_covs <- c(covs_cont, covs_bin)

df <- df[complete.cases(df[, covs_cont]), ]

cat(
  "N =", nrow(df),
  "| Treated =", sum(df$treat),
  "| Control =", sum(1L - df$treat),
  "\n"
)

# Balance table
fmla_ps <- as.formula(
  paste("treat ~", paste(all_covs, collapse = " + "))
)

balance_table <- map_dfr(all_covs, function(v) {
  
  ctrl_vals <- df[[v]][df$treat == 0]
  trt_vals  <- df[[v]][df$treat == 1]
  
  tibble(
    Variable  = v,
    Ctrl_Mean = mean(ctrl_vals, na.rm = TRUE),
    Trt_Mean  = mean(trt_vals, na.rm = TRUE),
    Std_Diff  = (
      mean(trt_vals, na.rm = TRUE) -
        mean(ctrl_vals, na.rm = TRUE)
    ) /
      sqrt(
        (
          var(ctrl_vals, na.rm = TRUE) +
            var(trt_vals, na.rm = TRUE)
        ) / 2
      ),
    p_value = t.test(ctrl_vals, trt_vals)$p.value
  )
  
})

print(balance_table, n = 44)

# OLS without controls
ols0 <- feols(
  died ~ treat,
  data = df,
  vcov = "HC3"
)

etable(ols0)

# OLS with controls
fmla_ols <- as.formula(
  paste("died ~ treat +", paste(all_covs, collapse = " + "))
)

ols1 <- feols(
  fmla_ols,
  data = df,
  vcov = "HC3"
)

etable(ols1)

etable(
  ols0,
  ols1,
  vcov = "HC3",
  dict = c(treat = "RHC (treated)"),
  title = "OLS estimates of RHC effect on 180-day mortality"
)

# IPW
set.seed(42)

ipw_fit <- treatweight(
  y    = df$died,
  d    = df$treat,
  x    = as.matrix(df[, all_covs]),
  trim = 0.05,
  boot = 999
)

cat("\n=== IPW (treatweight) ===\n")
cat("ATE:  ", round(ipw_fit$effect, 4), "\n")
cat("SE:   ", round(ipw_fit$se,     4), "\n")
cat("95% CI:", round(ipw_fit$effect - 1.96*ipw_fit$se, 4),
    "to", round(ipw_fit$effect + 1.96*ipw_fit$se, 4), "\n")

# Propensity score distribution
df$ps <- glm(
  fmla_ps,
  data = df,
  family = binomial
)$fitted.values

ggplot(df, aes(x = ps, fill = factor(treat))) +
  geom_histogram(
    bins = 40,
    alpha = 0.65,
    position = "identity"
  ) +
  scale_fill_manual(
    values = c("#d62728", "#1f77b4"),
    labels = c("Control (No RHC)", "Treated (RHC)")
  ) +
  labs(
    title = "Figure 1: Propensity Score Distribution",
    x = "Estimated Propensity Score",
    y = "Count",
    fill = ""
  ) +
  theme_bw()

# Balance plot
df$ipw_weight <- ifelse(
  df$treat == 1,
  1 / df$ps,
  1 / (1 - df$ps)
)

bal_after <- map_dfr(all_covs, function(v) {
  
  ctrl_w <- df[[v]][df$treat == 0]
  trt_w  <- df[[v]][df$treat == 1]
  
  wc <- df$ipw_weight[df$treat == 0]
  wt <- df$ipw_weight[df$treat == 1]
  
  mc_w <- weighted.mean(ctrl_w, wc, na.rm = TRUE)
  mt_w <- weighted.mean(trt_w, wt, na.rm = TRUE)
  
  sc <- sd(ctrl_w, na.rm = TRUE)
  st <- sd(trt_w, na.rm = TRUE)
  
  tibble(
    Variable = v,
    Std_Diff_After = (mt_w - mc_w) /
      sqrt((sc^2 + st^2) / 2)
  )
  
})

bal_plot <- balance_table %>%
  select(Variable, Std_Diff) %>%
  rename(Before = Std_Diff) %>%
  left_join(
    bal_after %>% rename(After = Std_Diff_After),
    by = "Variable"
  ) %>%
  pivot_longer(
    c(Before, After),
    names_to = "Timing",
    values_to = "Std_Diff"
  )

ggplot(
  bal_plot,
  aes(
    x = Std_Diff,
    y = Variable,
    color = Timing,
    shape = Timing
  )
) +
  geom_point(size = 2) +
  geom_vline(
    xintercept = c(-0.1, 0.1),
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_vline(
    xintercept = 0,
    color = "black",
    linewidth = 0.5
  ) +
  scale_color_manual(
    values = c(
      Before = "#d62728",
      After = "#1f77b4"
    )
  ) +
  scale_shape_manual(
    values = c(
      Before = 16,
      After = 18
    )
  ) +
  labs(
    title = "Figure 2: Covariate Balance Before and After IPW",
    x = "Standardised Difference",
    y = NULL,
    color = "",
    shape = ""
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 7)
  )

# AIPW Double ML
set.seed(42)

dml_data <- data.table(
  died = df$died,
  treat = df$treat,
  df[, all_covs]
)

dml_obj <- DoubleMLData$new(
  data = dml_data,
  y_col = "died",
  d_cols = "treat",
  x_cols = all_covs
)

ml_g <- lrn(
  "regr.ranger",
  num.trees = 500
)

ml_m <- lrn(
  "classif.ranger",
  num.trees = 500,
  predict_type = "prob"
)

dml_fit <- DoubleMLIRM$new(
  data = dml_obj,
  ml_g = ml_g,
  ml_m = ml_m,
  n_folds = 5,
  score = "ATE"
)

dml_fit$fit()

cat("\n=== AIPW Double ML (DoubleML) ===\n")
cat("ATE: ", round(as.numeric(dml_fit$coef), 4), "\n")
cat("SE: ", round(as.numeric(dml_fit$se), 4), "\n")

cat(
  "95% CI:",
  round(as.numeric(dml_fit$coef) - 1.96 * as.numeric(dml_fit$se), 4),
  "to",
  round(as.numeric(dml_fit$coef) + 1.96 * as.numeric(dml_fit$se), 4),
  "\n"
)

# CATE
set.seed(42)

X <- as.matrix(df[, all_covs])

cf <- causal_forest(
  X = X,
  Y = df$died,
  W = df$treat,
  num.trees = 2000
)

df$cate_pred <- predict(cf)$predictions

df$age_group <- cut(
  df$age,
  breaks = c(0, 50, 65, 80, Inf),
  labels = c("< 50", "50–65", "65–80", "> 80")
)

cate_age <- df %>%
  group_by(age_group) %>%
  summarise(
    CATE = mean(cate_pred),
    SE = sd(cate_pred) / sqrt(n()),
    lower = CATE - 1.96 * SE,
    upper = CATE + 1.96 * SE,
    n = n(),
    .groups = "drop"
  )

print(cate_age)

ggplot(
  cate_age,
  aes(
    x = age_group,
    y = CATE,
    ymin = lower,
    ymax = upper,
    fill = CATE > 0
  )
) +
  geom_col(width = 0.5, alpha = 0.75) +
  geom_errorbar(width = 0.2, linewidth = 1) +
  geom_hline(yintercept = 0) +
  scale_fill_manual(
    values = c("steelblue", "#d62728"),
    guide = "none"
  ) +
  labs(
    title = "Figure 3: Estimated CATE by Age Group",
    x = "Age Group",
    y = "Estimated CATE"
  ) +
  theme_bw()

# Results summary
results <- tibble(
  Method = c(
    "OLS (no controls)",
    "OLS (with controls)",
    "IPW",
    "AIPW Double ML"
  ),
  ATE = c(
    coef(ols0)["treat"],
    coef(ols1)["treat"],
    ipw_fit$effect,
    as.numeric(dml_fit$coef)
  ),
  SE = c(
    se(ols0)["treat"],
    se(ols1)["treat"],
    ipw_fit$se,
    as.numeric(dml_fit$se)
  )
) %>%
  mutate(
    CI_low = ATE - 1.96 * SE,
    CI_high = ATE + 1.96 * SE
  )
print(results, digits = 4)

# Forest plot
ggplot(
  results,
  aes(
    x = ATE,
    xmin = CI_low,
    xmax = CI_high,
    y = fct_rev(Method),
    color = Method == "AIPW Double ML"
  )
) +
  geom_pointrange(
    size = 0.8,
    linewidth = 1.2
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  scale_color_manual(
    values = c("grey40", "#1f77b4"),
    guide = "none"
  ) +
  labs(
    title = "Figure 4: Summary of Causal Effect EST",
    subtitle = "Point estimates with 95% confidence intervals",
    x = "ATE (probability scale)",
    y = ""
  ) +
  theme_bw()