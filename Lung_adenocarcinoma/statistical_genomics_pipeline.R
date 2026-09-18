# ==============================================================================
# PIPELINE: Translational Oncology & Biomarker Diagnostics Pipeline
# Context: KRAS Amplification & Sotorasib Efficacy in Lung Adenocarcinoma
# Outputs: 1 CSV Data Matrix, 4 Console & PNG Charts, 1 Formatted Word Report
# ==============================================================================

# 1. ENVIRONMENT CONFIGURATION & AUTOMATED PACKAGES CHECK
required_packages <- c("ggplot2", "officer", "flextable")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(ggplot2)
library(officer)
library(flextable)

# 2. CLINICAL ONCOLOGY COHORT DATASET
n_patients <- 200

# Categorical Clinical Metrics
clinical_stage <- sample(c("Stage_I", "Stage_II", "Stage_III"), n_patients, replace = TRUE)
treatment <- sample(c("Standard_Chemo", "Sotorasib_Targeted"), n_patients, replace = TRUE)
sex <- sample(c("Male", "Female"), n_patients, replace = TRUE)

# Continuous Demographic Covariate
age <- round(runif(n_patients, min = 40, max = 78))

# Transcriptomic Expressions (Biologically intertwined pathways)
# KRAS expression driven higher by advanced disease stage severity
kras_expression <- 4.2 + (ifelse(clinical_stage == "Stage_II", 1.2, 0) + 
                            ifelse(clinical_stage == "Stage_III", 2.9, 0)) + rnorm(n_patients, 0, 0.8)

# MAPK1 is downstream of KRAS and shares redundant biological pathway information
mapk1_expression <- (kras_expression * 0.65) + rnorm(n_patients, 1.2, 0.5)

# ACTB is a standard constant housekeeping control gene (independent baseline noise)
actb_expression <- rnorm(n_patients, mean = 7.0, sd = 0.4) 

# Continuous Phenotypic Outcome: Macroscopic Tumor Volume (mm^3)
# Sotorasib targeted drug therapy severely halts tumor volume specifically in late Stage III disease
tumor_volume <- 200 + (age * 1.5) + (kras_expression * 45) + 
  (ifelse(clinical_stage == "Stage_III" & treatment == "Sotorasib_Targeted", -120, 0)) + 
  rnorm(n_patients, 0, 25)

# Compile into master Clinical Genomics Matrix
lung_cancer_cohort <- data.frame(
  Patient_ID = sprintf("LUNG_%04d", 1:n_patients),
  Age = age,
  Sex = as.factor(sex),
  Clinical_Stage = factor(clinical_stage, levels = c("Stage_I", "Stage_II", "Stage_III")),
  Treatment_Arm = as.factor(treatment),
  KRAS_Expression = kras_expression,
  MAPK1_Expression = mapk1_expression,
  ACTB_Control = actb_expression,
  Tumor_Volume_mm3 = tumor_volume
)

# Export copy of data for reproducibility
write.csv(lung_cancer_cohort, "lung_adenocarcinoma_cohort.csv", row.names = FALSE)
cat("SUCCESS: 'lung_adenocarcinoma_cohort.csv' generated.\n\n")


# ==============================================================================
# 3. PERFORM STATISTICAL CALCULATIONS & FORMAT NUMBERS
# ==============================================================================

# Analysis 1: Two-Way ANOVA 
anova_model <- aov(KRAS_Expression ~ Clinical_Stage * Treatment_Arm, data = lung_cancer_cohort)
anova_raw <- (summary(anova_model))[[1]] # Extract underlying statistical table safely

anova_df <- data.frame(
  Source   = trimws(rownames(anova_raw)),
  Df       = anova_raw[["Df"]],
  Sum_Sq   = round(anova_raw[["Sum Sq"]], 3),
  Mean_Sq  = round(anova_raw[["Mean Sq"]], 3),
  F_Value  = round(anova_raw[["F value"]], 3),
  P_Value  = format.pval(anova_raw[["Pr(>F)"]], digits = 3, eps = 0.001)
)
# Clean row names mapping where R leaves trailing characters or NAs
anova_df$F_Value[is.na(anova_df$F_Value)] <- ""
anova_df$P_Value[anova_df$P_Value == "NA"] <- ""

# Analysis 2: Multiple Linear Regression
regression_model <- lm(Tumor_Volume_mm3 ~ KRAS_Expression + MAPK1_Expression + ACTB_Control + Age + Sex, data = lung_cancer_cohort)
reg_summary <- summary(regression_model)
reg_raw <- reg_summary$coefficients

reg_df <- data.frame(
  Predictor = rownames(reg_raw),
  Estimate  = round(reg_raw[, "Estimate"], 3),
  Std_Error = round(reg_raw[, "Std. Error"], 3),
  T_Value   = round(reg_raw[, "t value"], 3),
  P_Value   = format.pval(reg_raw[, "Pr(>|t|)"], digits = 3, eps = 0.001)
)

# Extract Model Predictions and Residuals for diagnostics
lung_cancer_cohort$Predicted_Volume <- predict(regression_model)
lung_cancer_cohort$Residuals <- residuals(regression_model)

# Analysis 3: Diagnostics (Variance Inflation Factor Calculation)
X <- model.matrix(regression_model)[,-1]
vif_values <- sapply(1:ncol(X), function(i) 1 / (1 - summary(lm(X[,i] ~ X[,-i]))$r.squared))
vif_df <- data.frame(
  Biomarker = colnames(X), 
  VIF_Value = round(unname(vif_values), 3)
)


# ==============================================================================
# 4. GENERATE FULL VISUAL SUITE & DISPLAY TO CONSOLE (GGPLOT2)
# ==============================================================================

# Universal canvas theme to prevent cut-off text or clipped boundaries
safe_margins <- theme_minimal(base_size = 11) + 
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, vjust = 1.5, size = 11),
    axis.title.x = element_text(vjust = -0.5, size = 10),
    axis.title.y = element_text(vjust = 1.5, size = 10),
    plot.margin = margin(t = 20, r = 25, b = 20, l = 25), 
    legend.position = "bottom",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8)
  )

# Plot 1: Boxplot for ANOVA
p1 <- ggplot(lung_cancer_cohort, aes(x = Clinical_Stage, y = KRAS_Expression, fill = Treatment_Arm)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(0.8)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), size = 0.8, alpha = 0.3, color = "darkgray") +
  scale_fill_manual(values = c("Standard_Chemo" = "#E64B35FF", "Sotorasib_Targeted" = "#4DBBD5FF")) +
  labs(title = "KRAS Expression Across Lung Pathology", x = "Lung Cancer Stage", y = "Log2 KRAS Expression", fill = "Clinical Arm") +
  safe_margins
ggsave("plot1_anova_boxplot.png", plot = p1, width = 5.5, height = 3.8, dpi = 300)
print(p1)

# Plot 2: Actual vs Predicted Regression Fit
p2 <- ggplot(lung_cancer_cohort, aes(x = Predicted_Volume, y = Tumor_Volume_mm3)) +
  geom_point(color = "#3C5488FF", alpha = 0.5, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed", size = 0.8) +
  labs(title = "Prognostic Accuracy: Predicted vs. Actual", x = "Predicted Volume (mm³)", y = "Actual Observed Volume (mm³)") +
  safe_margins
ggsave("plot2_regression_fit.png", plot = p2, width = 5.5, height = 3.8, dpi = 300)
print(p2)

# Plot 3: Residuals Fitted Homoscedasticity Analysis
p3 <- ggplot(lung_cancer_cohort, aes(x = Predicted_Volume, y = Residuals)) +
  geom_point(color = "#00A087FF", alpha = 0.5, size = 1.5) +
  geom_hline(yintercept = 0, color = "red", linetype = "solid", size = 0.8) +
  labs(title = "Residual Homoscedasticity Verification", x = "Fitted (Predicted) Values", y = "Residual Model Error") +
  safe_margins
ggsave("plot3_residuals_diagnostic.png", plot = p3, width = 5.5, height = 3.8, dpi = 300)
print(p3)

# Plot 4: VIF Multicollinearity Bars
p4 <- ggplot(vif_df, aes(x = reorder(Biomarker, VIF_Value), y = VIF_Value, fill = VIF_Value > 4)) +
  geom_bar(stat = "identity", width = 0.5, show.legend = FALSE) +
  geom_hline(yintercept = 4, color = "red", linetype = "dashed", size = 0.8) +
  scale_fill_manual(values = c("FALSE" = "#4DBBD5FF", "TRUE" = "#E64B35FF")) +
  coord_flip() + 
  labs(title = "Collinearity Profiling (VIF Pathway Check)", x = "Model Predictor", y = "Variance Inflation Factor") +
  safe_margins
ggsave("plot4_vif_bars.png", plot = p4, width = 5.5, height = 3.8, dpi = 300)
print(p4)


# ==============================================================================
# 5. AUTOMATED MS WORD REPORT COMPILATION ENGINE
# ==============================================================================
make_clean_table <- function(df) {
  ft <- flextable(df) |> 
    theme_vanilla() |> 
    align(align = "center", part = "all") |> 
    autofit() |> 
    fit_to_width(max_width = 6.5) # Keeps tables safely locked inside margins
  return(ft)
}

doc <- read_docx() |>
  body_add_par("KRAS Targeted Therapeutics Analytics Report", style = "heading 1") |>
  body_add_par("Lung Adenocarcinoma Stratification Pipeline", style = "Normal") |>
  body_add_par("", style = "Normal") |>
  
  # Section 1: ANOVA
  body_add_par("1. Biomarker Expression Profiling (Two-Way ANOVA)", style = "heading 2") |>
  body_add_flextable(make_clean_table(anova_df)) |>
  body_add_par("", style = "Normal") |>
  body_add_img(src = "plot1_anova_boxplot.png", width = 5.0, height = 3.4) |>
  body_add_par("", style = "Normal") |>
  
  # Section 2: Multiple Regression
  body_add_par("2. Predictive Prognostics & Multiple Regression Outputs", style = "heading 2") |>
  body_add_flextable(make_clean_table(reg_df)) |>
  body_add_par(paste("Model Performance: Multiple R-squared =", round(reg_summary$r.squared, 4)), style = "Normal") |>
  body_add_par("", style = "Normal") |>
  body_add_img(src = "plot2_regression_fit.png", width = 5.0, height = 3.4) |>
  body_add_img(src = "plot3_residuals_diagnostic.png", width = 5.0, height = 3.4) |>
  body_add_par("", style = "Normal") |>
  
  # Section 3: VIF Diagnostics
  body_add_par("3. Pathway Redundancy Assessment (VIF Diagnostics)", style = "heading 2") |>
  body_add_flextable(make_clean_table(vif_df)) |>
  body_add_par("", style = "Normal") |>
  body_add_img(src = "plot4_vif_bars.png", width = 5.0, height = 3.4)

# Output final compiled document asset
print(doc, target = "Statistical_Genomics_Report.docx")

cat("\n========================================================================\n")
cat(" SUCCESS: ALL GRAPHICS PRINTED AND MARGIN-SAFE REPORT DOCUMENT SAVED!\n")
cat("========================================================================\n")
