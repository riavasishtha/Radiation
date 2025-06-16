

# Load required libraries
library(lmtp)
library(survival)
library(dplyr)
library(ggplot2)
library(nnls)

# Assuming your data is loaded as 'df'
clean_data <- read.csv("real_HNSCC_data.csv")


# Verify data structure
print("Data structure:")
str(clean_data)

# Define confounders (all numeric, using original variable names)
confounders <- c("Sex", "Age", "Smoking_PY", "Stage_numeric", 
                 "HPV_Positive", "HPV_Unknown", "Chemotherapy", 
                 "RT_year")

# Define treatment variable
treatment <- "BED_DD"

# Define outcome and censoring
outcome <- "Survival_time"
censor_var <- "Death_event"

# Create modified treatment policies
# For continuous treatments, shift functions should take (data, trt) and return modified trt

# Policy 1: Increase by 20%
policy_increase <- function(data, trt) {
  return(data[[trt]] * 1.20)
}

# Policy 2: Decrease by 20%
policy_decrease <- function(data, trt) {
  return(data[[trt]] * 0.80)
}

# Policy 3: Add 10 units
policy_add <- function(data, trt) {
  return(data[[trt]] + 10)
}

# Policy 4: Set to median
policy_median <- function(data, trt) {
  median_val <- median(data[[trt]], na.rm = TRUE)
  return(rep(median_val, nrow(data)))
}




# For survival analysis with lmtp, we need to create pseudo-outcomes
# One approach is to use restricted mean survival time or survival probability at fixed time

# Calculate restricted mean survival time (RMST) at a specific time point
# Choose a reasonable follow-up time 
tau <- 36  # months

# Create survival object
surv_obj <- Surv(clean_data$Survival_time, clean_data$Death_event)

# Fit Kaplan-Meier to get baseline survival function
km_fit <- survfit(surv_obj ~ 1, data = clean_data)

# Use survival probability at fixed time as binary outcome
clean_data$surv_indicator <- ifelse(clean_data$Survival_time >= tau & clean_data$Death_event == 0, 1,
                                    ifelse(clean_data$Survival_time >= tau & clean_data$Death_event == 1, 0,
                                           ifelse(clean_data$Survival_time < tau & clean_data$Death_event == 1, 0, NA)))

# Remove rows with missing survival indicator (administrative censoring before tau)
analysis_data <- clean_data[!is.na(clean_data$surv_indicator), ]


print(paste("Analysis dataset has", nrow(analysis_data), "observations"))
print(paste("Events observed:", sum(analysis_data$Death_event)))
print(paste("Survival rate at", tau, "months:", mean(analysis_data$surv_indicator, na.rm = TRUE)))

summary(analysis_data)

# LMTP Analysis for survival probability at tau months
# Using the survival indicator as binary outcome
# Set mtp = TRUE for continuous treatments

lmtp_fit_increase <- lmtp_tmle(
  data = analysis_data,
  trt = treatment,
  outcome = "surv_indicator",
  baseline = confounders,
  shift = policy_increase,
  outcome_type = "binomial",
  mtp = TRUE,                   # Important: Set to TRUE for continuous treatments
  learners_trt = "SL.glm",      # Treatment model
  learners_outcome = "SL.glm"   # Outcome model
)

lmtp_fit_decrease <- lmtp_tmle(
  data = analysis_data,
  trt = treatment,
  outcome = "surv_indicator", 
  baseline = confounders,
  shift = policy_decrease,
  outcome_type = "binomial",
  mtp = TRUE,                   # Important: Set to TRUE for continuous treatments
  learners_trt = "SL.glm",
  learners_outcome = "SL.glm"
)


lmtp_fit_median <- lmtp_tmle(
  data = analysis_data,
  trt = treatment,
  outcome = "surv_indicator", 
  baseline = confounders,
  shift = policy_median,
  outcome_type = "binomial",
  mtp = TRUE,                   # Important: Set to TRUE for continuous treatments
  learners_trt = "SL.glm",
  learners_outcome = "SL.glm"
)

# Additive shift policy
lmtp_fit_add <- lmtp_tmle(
  data = analysis_data,
  trt = treatment,
  outcome = "surv_indicator",
  baseline = confounders, 
  shift = policy_add,
  outcome_type = "binomial",
  mtp = TRUE,
  learners_trt = "SL.glm",
  learners_outcome = "SL.glm"
)

# Natural course (observed treatment)
lmtp_fit_natural <- lmtp_tmle(
  data = analysis_data,
  trt = treatment,
  outcome = "surv_indicator",
  baseline = confounders,
  shift = NULL,  # Natural course
  outcome_type = "binomial",
  mtp = TRUE,    # Set to TRUE even for natural course with continuous treatment
  learners_trt = "SL.glm",
  learners_outcome = "SL.glm"
)

# Extract results using the nested estimate object
results_summary <- data.frame(
  Policy = c("Observed BED_DD", "Median BED_DD", "BED_DD +20%", "BED_DD -20%", "BED_DD +10 units"),
  Survival_Prob = c(
    lmtp_fit_natural$estimate@x,
    lmtp_fit_median$estimate@x,
    lmtp_fit_increase$estimate@x,
    lmtp_fit_decrease$estimate@x,
    lmtp_fit_add$estimate@x
  ),
  SE = c(
    lmtp_fit_natural$estimate@std_error,
    lmtp_fit_median$estimate@std_error,
    lmtp_fit_increase$estimate@std_error,
    lmtp_fit_decrease$estimate@std_error,
    lmtp_fit_add$estimate@std_error
  )
)

# Compute 95% confidence intervals
results_summary$CI_Lower <- results_summary$Survival_Prob - 1.96 * results_summary$SE
results_summary$CI_Upper <- results_summary$Survival_Prob + 1.96 * results_summary$SE

# Print results
print("LMTP Results - Survival Probability at tau months:")
print(results_summary)

ggplot(results_summary, aes(x = Policy, y = Survival_Prob)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2) +
  labs(title = "Estimated Survival Probability at tau Months",
       y = "Survival Probability",
       x = "Treatment Policy") +
  theme_bw()


# Calculate causal contrasts
contrast_median <- lmtp_contrast(lmtp_fit_median, ref = lmtp_fit_natural, type = "additive")
contrast_increase <- lmtp_contrast(lmtp_fit_increase, ref = lmtp_fit_natural, type = "additive")
contrast_decrease <- lmtp_contrast(lmtp_fit_decrease, ref = lmtp_fit_natural, type = "additive")
contrast_add      <- lmtp_contrast(lmtp_fit_add,      ref = lmtp_fit_natural, type = "additive")

# Print contrast results
cat("\nCausal Contrasts:\n")

cat("Effect of median BED_DD vs natural course:\n")
cat(sprintf("Risk Difference: %.4f\n", contrast_median$vals$theta))
cat(sprintf("95%% CI: (%.4f, %.4f)\n\n", contrast_median$vals$conf.low, contrast_median$vals$conf.high))

cat("Effect of 20% BED_DD increase vs natural course:\n")
cat(sprintf("Risk Difference: %.4f\n", contrast_increase$vals$theta))
cat(sprintf("95%% CI: (%.4f, %.4f)\n\n", contrast_increase$vals$conf.low, contrast_increase$vals$conf.high))

cat("Effect of 20% BED_DD decrease vs natural course:\n")
cat(sprintf("Risk Difference: %.4f\n", contrast_decrease$vals$theta))
cat(sprintf("95%% CI: (%.4f, %.4f)\n\n", contrast_decrease$vals$conf.low, contrast_decrease$vals$conf.high))

cat("Effect of BED_DD +10 units vs natural course:\n")
cat(sprintf("Risk Difference: %.4f\n", contrast_add$vals$theta))
cat(sprintf("95%% CI: (%.4f, %.4f)\n\n", contrast_add$vals$conf.low, contrast_add$vals$conf.high))


# Plot 1: Survival probabilities under each policy
surv_plot <- ggplot(results_summary, aes(x = Policy, y = Survival_Prob)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2) +
  ylim(0, 1) + 
  labs(
    title = paste("Estimated Survival Probability at", tau, "Months"),
    subtitle = "Causal Effects of BED_DD Treatment Policies",
    y = "Survival Probability",
    x = "Treatment Policy"
  ) +
  theme_bw()  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save Plot 1
ggsave("survival_probabilities.png", plot = surv_plot, dpi = 600, width = 8, height = 6)

# Summary statistics of treatment variable
cat("\nBED_DD Treatment Variable Summary:\n")
print(summary(clean_data$BED_DD))

cat("\nTreatment distribution by key confounders:\n")
cat("Sex distribution:\n")
print(table(clean_data$Sex, useNA = "ifany"))
cat("BED_DD by Sex:\n")
print(by(clean_data$BED_DD, clean_data$Sex, summary))


# Additional analysis: Dose-response curve
bed_levels <- seq(quantile(clean_data$BED_DD, 0.05, na.rm = TRUE),
                  quantile(clean_data$BED_DD, 0.95, na.rm = TRUE),
                  length.out = 20)

dose_response_results <- data.frame(
  BED_DD_Level = bed_levels,
  Survival_Prob = NA,
  SE = NA
)

for(i in seq_along(bed_levels)) {
  policy_fixed <- function(data, trt) rep(bed_levels[i], length(trt))
  
  fit_temp <- lmtp_tmle(
    data = analysis_data,
    trt = treatment,
    outcome = "surv_indicator",
    baseline = confounders,
    shift = policy_fixed,
    outcome_type = "binomial",
    mtp = TRUE,
    learners_trt = "SL.glm",
    learners_outcome = "SL.glm"
  )
  
  dose_response_results$Survival_Prob[i] <- fit_temp$estimate@x
  dose_response_results$SE[i] <- fit_temp$estimate@std_error
}

dose_response_results$CI_Lower <- dose_response_results$Survival_Prob - 1.96 * dose_response_results$SE
dose_response_results$CI_Upper <- dose_response_results$Survival_Prob + 1.96 * dose_response_results$SE

# Plot 2: Dose-response curve
dose_plot <- ggplot(dose_response_results, aes(x = BED_DD_Level, y = Survival_Prob)) +
  geom_line() +
  geom_point() +
  geom_ribbon(aes(ymin = CI_Lower, ymax = CI_Upper), alpha = 0.3) +
  labs(
    title = "Dose-Response Curve: BED_DD vs Survival Probability",
    x = "BED_DD Level",
    y = paste("Survival Probability at", tau, "Months")
  ) + 
  theme_bw()

# Save Plot 2
ggsave("dose_response_curve.png", plot = dose_plot, dpi = 600, width = 8, height = 6)

# Print final dose-response data
cat("\nDose-Response Results:\n")
print(dose_response_results)



