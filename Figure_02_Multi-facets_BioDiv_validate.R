# READNE: validate BioDiv generated from eDNA datasets, 09.04.2025
# Created by Yan Zhang (Nanjing University, China)
# ALL RIGHTS RESERVED
wd <- '../data/'
setwd(wd)

fun_path <- '../code/00_functions.R' ## CHANGE ME
source(fun_path)

useCors <- c('#D42985', '#0D9248','#0DB0AA','#E18B27')
all_BioDiv_para <- c("SR_Re","Redundancy_Re","MNTD_Re","GD_Re",
                     "SR_Ori","Redundancy_Ori","MNTD_Tre","GDvalue.COI_Tre")
col_num <- rep(1:4, 2)

fac_names <- c("Reanalyzed Species Richness", "Reanalyzed Func. Redundancy", "ASV-based MNTD", "ASV-based Genetic Diversity",
               "Reported Species Richness", "Reported Func. Redundancy", "Tree-based MNTD", "Barcode-derived Genetic Diversity"  )
#############################################################
####A01: Validate Taxon diversity: gamma
#############################################################
BioDiv <- read.csv('Use_All_Rivers_Multifacet_Div_all.csv')
head(BioDiv)

spepool <- read.csv('Use_SpePool_Basin.csv')
head(spepool)

data <- merge.data.frame(BioDiv[c(2,8,9,13,15,19)], spepool, by = 'RivID')
data$SR_Re <- log10(data$SR_Re+1)
data$SR_Ori <- log10(data$SR_Ori+1)
data$nSpes <- log10(data$nSpes+1)

# General function to fit GLMMs and generate plots
fit_models_plot <- function(response_var, rand_var, data, color) {
  form_full <- as.formula(paste0(response_var, " ~ nSpes + TotVolume + (1|", rand_var, ")"))
  form_null <- as.formula(paste0(response_var, " ~ TotVolume + (1|", rand_var, ")"))
  
  data_sub <- data[!is.na(data[[response_var]]), ]
  
  fit_full <- glmmTMB(form_full, data = data_sub)
  fit_null <- glmmTMB(form_null, data = data_sub)
  
  r2_full <- r.squaredGLMM(fit_full)
  r2_null <- r.squaredGLMM(fit_null)
  sig_test <- anova(fit_full, fit_null)
  
  # plotting
  pval <- sig_test[2, 8]
  r2_diff <- r2_full[1] - r2_null[1]
  lab <- if (pval < 0.001) {
    paste0('Marginal R² = ', sprintf('%.02f', r2_diff), ' p < 0.001')
  } else {
    paste0('Marginal R² = ', sprintf('%.02f', r2_diff), ' p = ', sprintf('%.03f', pval))
  }
  
  data_sub$fitted_y <- predict(fit_full)
  
  x = min(data_sub$nSpes) + (max(data_sub$nSpes) - min(data_sub$nSpes)) / 4
  y = max(data_sub[[response_var]]) - (max(data_sub[[response_var]]) - min(data_sub[[response_var]])) / 10
  
  p <- ggplot(data_sub, aes(x = nSpes, y = .data[[response_var]])) +
    geom_point(alpha = 0.8, color = color) +
    labs(y = 'Species richness', x = "Species pool") +
    scale_x_continuous(breaks = c(1,2,3), labels = c(10, 100, 1000)) +
    scale_y_continuous(breaks = c(0,1,2), labels = c(1, 10, 100)) +
    geom_smooth(mapping = aes(x = nSpes, y = fitted_y), method = 'glm', formula = y ~ x, se = TRUE,
                lwd = 1, alpha = 0.1, color = color, fill = color) +
    annotate('text', x = x, y = y, label = lab, size = 3) +
    theme_bw() +
    theme(panel.grid = element_blank(), plot.background = element_blank())
  
  return(list(
    model_full = fit_full,
    model_null = fit_null,
    plot = p,
    r2_diff = r2_diff,
    p_value = pval
  ))
}

# Loop over all combinations
results <- list()
param_table <- data.frame()
responses <- c("SR_Re", "SR_Ori")
randoms <- c("Primer", "BarcodeRegion")

for (resp in responses) {
  for (rand in randoms) {
    tag <- paste0(resp, "|", rand)
    res <- fit_models_plot(response_var = resp, rand_var = rand, data = data, color = useCors[1])
    
    # Store results
    results[[tag]] <- list(full = res$model_full, null = res$model_null, plot = res$plot)
    
    # Store summary statistics
    summary_fit <- summary(res$model_full)
    coefs <- summary_fit$coefficients$cond
    
    param_table <- rbind(param_table, data.frame(
      Response = resp,
      Random = rand,
      Intercept_Est = round(coefs["(Intercept)", "Estimate"], 3),
      Intercept_SE = round(coefs["(Intercept)", "Std. Error"], 3),
      Intercept_p = round(coefs["(Intercept)", "Pr(>|z|)"], 4),
      nSpes_Est = round(coefs["nSpes", "Estimate"], 3),
      nSpes_SE = round(coefs["nSpes", "Std. Error"], 3),
      nSpes_p = round(coefs["nSpes", "Pr(>|z|)"], 4),
      TotVol_Est = round(coefs["TotVolume", "Estimate"], 3),
      TotVol_SE = round(coefs["TotVolume", "Std. Error"], 3),
      TotVol_p = round(coefs["TotVolume", "Pr(>|z|)"], 4),
      Marginal_R2 = round(r.squaredGLMM(res$model_full)[1], 3),
      Conditional_R2 = round(r.squaredGLMM(res$model_full)[2], 3),
      AIC = AIC(res$model_full),
      LogLik = logLik(res$model_full)[1],
      P_value = res$p_value,
      N_obs = nobs(res$model_full)
    ))
  }
}

saveRDS(results, "../result/Res_01_Validate_TD_all.RDS")
write.csv(param_table, "../result/Res_01_Validate_TD_all.csv", row.names = FALSE)

# Example: view plots individually
results[["SR_Re|Primer"]]$plot
results[["SR_Re|BarcodeRegion"]]$plot
results[["SR_Ori|Primer"]]$plot
results[["SR_Ori|BarcodeRegion"]]$plot

#############################################################
####A02: Validate Phylo diversity: gamma
#############################################################
BioDiv <- read.csv('Use_All_Rivers_Multifacet_Div_all.csv')
head(BioDiv)

data <- BioDiv[,c(1:14,17,21)]
data <- na.omit(data)
data[,15:16] <- log10(data[,15:16])

# define model
fit_models_plot_xy <- function(y_var, x_var, rand_var, data, color) {
  # Build model formulas
  form_full <- as.formula(paste0(y_var, " ~ ", x_var, " + TotVolume + (1|", rand_var, ")"))
  form_null <- as.formula(paste0(y_var, " ~ TotVolume + (1|", rand_var, ")"))
  
  # Remove NA values
  data_sub <- data[!is.na(data[[y_var]]) & !is.na(data[[x_var]]), ]
  
  # Fit models
  fit_full <- glmmTMB(form_full, data = data_sub)
  fit_null <- glmmTMB(form_null, data = data_sub)
  
  # Model stats
  r2_full <- r.squaredGLMM(fit_full)
  r2_null <- r.squaredGLMM(fit_null)
  sig_test <- anova(fit_full, fit_null)
  
  r2_diff <- r2_full[1] - r2_null[1]
  pval <- sig_test[2, 8]
  lab <- if (pval < 0.001) {
    paste0('Marginal R² = ', sprintf('%.02f', r2_diff), ' p < 0.001')
  } else {
    paste0('Marginal R² = ', sprintf('%.02f', r2_diff), ' p = ', sprintf('%.03f', pval))
  }
  
  # Predictions
  data_sub$fitted_y <- predict(fit_full)
  
  # Annotation position
  x = min(data_sub[[x_var]]) + (max(data_sub[[x_var]]) - min(data_sub[[x_var]])) / 4
  y = max(data_sub[[y_var]]) - (max(data_sub[[y_var]]) - min(data_sub[[y_var]])) / 10
  
  # Plot
  p <- ggplot(data_sub, aes_string(x = x_var, y = y_var)) +
    geom_point(alpha = 0.8, color = color) +
    labs(y = y_var, x = x_var) +
    geom_smooth(aes(y = fitted_y), method = 'glm', formula = y ~ x, se = TRUE,
                lwd = 1, alpha = 0.1, color = color, fill = color) +
    annotate('text', x = x, y = y, label = lab, size = 3) +
    theme_bw() +
    theme(panel.grid = element_blank(), plot.background = element_blank())
  
  # Return list
  return(list(
    model_full = fit_full,
    model_null = fit_null,
    plot = p,
    r2_diff = r2_diff,
    p_value = pval,
    summary = summary(fit_full)
  ))
}
# Define variables
y_vars <- c("MNTD_Re")
x_vars <- c("MNTD_Tre")
randoms <- c("Primer", "BarcodeRegion")

results <- list()
param_table <- data.frame()

for (i in seq_along(y_vars)) {
  y <- y_vars[i]
  x <- x_vars[i]
  
  for (rand in randoms) {
    tag <- paste0(y, "|", rand)
    res <- fit_models_plot_xy(y_var = y, x_var = x, rand_var = rand, data = data, color = useCors[3])
    
    # Store models and plot
    results[[tag]] <- list(
      full = res$model_full,
      null = res$model_null,
      plot = res$plot
    )
    
    # Get fixed effects
    coefs <- res$summary$coefficients$cond
    
    # Add to parameter table
    param_table <- rbind(param_table, data.frame(
      Response = y,
      Predictor = x,
      Random = rand,
      Intercept_Est = round(coefs["(Intercept)", "Estimate"], 3),
      Intercept_SE = round(coefs["(Intercept)", "Std. Error"], 3),
      Intercept_p = round(coefs["(Intercept)", "Pr(>|z|)"], 4),
      X_Est = round(coefs[x, "Estimate"], 3),
      X_SE = round(coefs[x, "Std. Error"], 3),
      X_p = round(coefs[x, "Pr(>|z|)"], 4),
      TotVol_Est = round(coefs["TotVolume", "Estimate"], 3),
      TotVol_SE = round(coefs["TotVolume", "Std. Error"], 3),
      TotVol_p = round(coefs["TotVolume", "Pr(>|z|)"], 4),
      Marginal_R2 = round(r.squaredGLMM(res$model_full)[1], 3),
      Conditional_R2 = round(r.squaredGLMM(res$model_full)[2], 3),
      AIC = AIC(res$model_full),
      LogLik = logLik(res$model_full)[1],
      P_value = res$p_value,
      N_obs = nobs(res$model_full)
    ))
  }
}

saveRDS(results, "../result/Res_01_Validate_PD_gamma.RDS")
write.csv(param_table, "../result/Res_01_Validate_PD_gamma.csv", row.names = FALSE)

#############################################################
####A02: Validate Phylo diversity: alpha
#############################################################
BioDiv <- read.csv('Use_All_Sites_Multifacet_Div_all.csv')
head(BioDiv)

data <- BioDiv[,c(1:17,20,24)]
data <- na.omit(data)
data[,18:19] <- log10(data[,18:19])

# define model
fit_models_plot_xy <- function(y_var, x_var, rand_var, data, color) {
  # Build model formulas
  form_full <- as.formula(paste0(y_var, " ~ ", x_var, " + SamVolume + (1|", rand_var, ")"))
  form_null <- as.formula(paste0(y_var, " ~ SamVolume + (1|", rand_var, ")"))
  
  # Remove NA values
  data_sub <- data[!is.na(data[[y_var]]) & !is.na(data[[x_var]]), ]
  
  # Fit models
  fit_full <- glmmTMB(form_full, data = data_sub)
  fit_null <- glmmTMB(form_null, data = data_sub)
  
  # Model stats
  r2_full <- r.squaredGLMM(fit_full)
  r2_null <- r.squaredGLMM(fit_null)
  sig_test <- anova(fit_full, fit_null)
  
  r2_diff <- r2_full[1] - r2_null[1]
  pval <- sig_test[2, 8]
  lab <- if (pval < 0.001) {
    paste0('Marginal R² = ', sprintf('%.02f', r2_diff), ' p < 0.001')
  } else {
    paste0('Marginal R² = ', sprintf('%.02f', r2_diff), ' p = ', sprintf('%.03f', pval))
  }
  
  # Predictions
  data_sub$fitted_y <- predict(fit_full)
  
  # Annotation position
  x = min(data_sub[[x_var]]) + (max(data_sub[[x_var]]) - min(data_sub[[x_var]])) / 4
  y = max(data_sub[[y_var]]) - (max(data_sub[[y_var]]) - min(data_sub[[y_var]])) / 10
  
  # Plot
  p <- ggplot(data_sub, aes_string(x = x_var, y = y_var)) +
    geom_point(alpha = 0.8, color = color) +
    labs(y = y_var, x = x_var) +
    geom_smooth(aes(y = fitted_y), method = 'glm', formula = y ~ x, se = TRUE,
                lwd = 1, alpha = 0.1, color = color, fill = color) +
    annotate('text', x = x, y = y, label = lab, size = 3) +
    theme_bw() +
    theme(panel.grid = element_blank(), plot.background = element_blank())
  
  # Return list
  return(list(
    model_full = fit_full,
    model_null = fit_null,
    plot = p,
    r2_diff = r2_diff,
    p_value = pval,
    summary = summary(fit_full)
  ))
}
# Define variables
y_vars <- c("MNTD_Re")
x_vars <- c("MNTD_Tre")
randoms <- c("Primer", "BarcodeRegion")

results <- list()
param_table <- data.frame()

for (i in seq_along(y_vars)) {
  y <- y_vars[i]
  x <- x_vars[i]
  
  for (rand in randoms) {
    tag <- paste0(y, "|", rand)
    res <- fit_models_plot_xy(y_var = y, x_var = x, rand_var = rand, data = data, color = useCors[3])
    
    # Store models and plot
    results[[tag]] <- list(
      full = res$model_full,
      null = res$model_null,
      plot = res$plot
    )
    
    # Get fixed effects
    coefs <- res$summary$coefficients$cond
    
    # Add to parameter table
    param_table <- rbind(param_table, data.frame(
      Response = y,
      Predictor = x,
      Random = rand,
      Intercept_Est = round(coefs["(Intercept)", "Estimate"], 3),
      Intercept_SE = round(coefs["(Intercept)", "Std. Error"], 3),
      Intercept_p = round(coefs["(Intercept)", "Pr(>|z|)"], 4),
      X_Est = round(coefs[x, "Estimate"], 3),
      X_SE = round(coefs[x, "Std. Error"], 3),
      X_p = round(coefs[x, "Pr(>|z|)"], 4),
      SamVol_Est = round(coefs["SamVolume", "Estimate"], 3),
      SamVol_SE = round(coefs["SamVolume", "Std. Error"], 3),
      SamVol_p = round(coefs["SamVolume", "Pr(>|z|)"], 4),
      Marginal_R2 = round(r.squaredGLMM(res$model_full)[1], 3),
      Conditional_R2 = round(r.squaredGLMM(res$model_full)[2], 3),
      AIC = AIC(res$model_full),
      LogLik = logLik(res$model_full)[1],
      P_value = res$p_value,
      N_obs = nobs(res$model_full)
    ))
  }
}

saveRDS(results, "../result/Res_01_Validate_PD_alpha.RDS")
write.csv(param_table, "../result/Res_01_Validate_PD_alpha.csv", row.names = FALSE)
#############################################################
####A03: Validate gen diversity: gamma
#############################################################
BioDiv <- read.csv('Use_All_Rivers_Multifacet_Div_all.csv')
head(BioDiv)

data <- BioDiv[,c(1:14,18,22:24)]
data <- na.omit(data)
data[,c(15,18)] <- log10(data[,c(15,18)] + 1e-3)
data[,c(16,17)] <- log10(data[,c(16,17)] + 1e-4)

# define model
fit_models_plot_xy <- function(y_var, x_var, rand_var, data, color) {
  # Build model formulas
  form_full <- as.formula(paste0(y_var, " ~ ", x_var, " + TotVolume + (1|", rand_var, ")"))
  form_null <- as.formula(paste0(y_var, " ~ TotVolume + (1|", rand_var, ")"))
  
  # Remove NA values
  data_sub <- data[!is.na(data[[y_var]]) & !is.na(data[[x_var]]), ]
  
  # Fit models
  fit_full <- glmmTMB(form_full, data = data_sub)
  fit_null <- glmmTMB(form_null, data = data_sub)
  
  # Model stats
  r2_full <- r.squaredGLMM(fit_full)
  r2_null <- r.squaredGLMM(fit_null)
  sig_test <- anova(fit_full, fit_null)
  
  r2_diff <- r2_full[1] - r2_null[1]
  pval <- sig_test[2, 8]
  lab <- if (pval < 0.001) {
    paste0('Marginal R² = ', sprintf('%.02f', r2_diff), ' p < 0.001')
  } else {
    paste0('Marginal R² = ', sprintf('%.02f', r2_diff), ' p = ', sprintf('%.03f', pval))
  }
  
  # Predictions
  data_sub$fitted_y <- predict(fit_full)
  
  # Annotation position
  x = min(data_sub[[x_var]]) + (max(data_sub[[x_var]]) - min(data_sub[[x_var]])) / 4
  y = max(data_sub[[y_var]]) - (max(data_sub[[y_var]]) - min(data_sub[[y_var]])) / 10
  
  # Plot
  p <- ggplot(data_sub, aes_string(x = x_var, y = y_var)) +
    geom_point(alpha = 0.8, color = color) +
    labs(y = y_var, x = x_var) +
    geom_smooth(aes(y = fitted_y), method = 'glm', formula = y ~ x, se = TRUE,
                lwd = 1, alpha = 0.1, color = color, fill = color) +
    annotate('text', x = x, y = y, label = lab, size = 3) +
    theme_bw() +
    theme(panel.grid = element_blank(), plot.background = element_blank())
  
  # Return list
  return(list(
    model_full = fit_full,
    model_null = fit_null,
    plot = p,
    r2_diff = r2_diff,
    p_value = pval,
    summary = summary(fit_full)
  ))
}
# Define variables
y <- "GD_Re"
x_vars <- c("GD.12S_Tre", "GD.16S_Tre", "GDvalue.COI_Tre")
randoms <- c("Primer", "BarcodeRegion")

results <- list()
param_table <- data.frame()

for (i in seq_along(x_vars)) {
  x <- x_vars[i]
  
  for (rand in randoms) {
    tag <- paste0(y, "_", x, "_", rand)
    
    res <- fit_models_plot_xy(y_var = y, x_var = x, rand_var = rand, data = data, color = useCors[4])
    
    results[[tag]] <- list(full = res$model_full,null = res$model_null,plot = res$plot)
    
    coefs <- res$summary$coefficients$cond
    
    param_table <- rbind(param_table, data.frame(
      Response = y,
      Predictor = x,
      Random = rand,
      Intercept_Est = round(coefs["(Intercept)", "Estimate"], 3),
      Intercept_SE = round(coefs["(Intercept)", "Std. Error"], 3),
      Intercept_p = round(coefs["(Intercept)", "Pr(>|z|)"], 4),
      X_Est = round(coefs[x, "Estimate"], 3),
      X_SE = round(coefs[x, "Std. Error"], 3),
      X_p = round(coefs[x, "Pr(>|z|)"], 4),
      TotVol_Est = round(coefs["TotVolume", "Estimate"], 3),
      TotVol_SE = round(coefs["TotVolume", "Std. Error"], 3),
      TotVol_p = round(coefs["TotVolume", "Pr(>|z|)"], 4),
      Marginal_R2 = round(r.squaredGLMM(res$model_full)[1], 3),
      Conditional_R2 = round(r.squaredGLMM(res$model_full)[2], 3),
      AIC = AIC(res$model_full),
      LogLik = logLik(res$model_full)[1],
      P_value = ifelse(res$p_value < 0.001, "<0.001", round(res$p_value, 3)),
      N_obs = nobs(res$model_full)
    ))
  }
}
saveRDS(results, "../result/Res_01_Validate_GD_gamma.RDS")
write.csv(param_table, "../result/Res_01_Validate_GD_gamma.csv", row.names = FALSE)
#############################################################
####A03: Validate gen diversity: alpha
#############################################################
BioDiv <- read.csv('Use_All_Sites_Multifacet_Div_all.csv')
head(BioDiv)

data <- BioDiv[,c(1:17,21,25:27)]
data <- na.omit(data)
data[,c(18,21)] <- log10(data[,c(18,21)] + 1e-3)
data[,c(19,20)] <- log10(data[,c(19,20)] + 1e-4)

# define model
fit_models_plot_xy <- function(y_var, x_var, rand_var, data, color) {
  # Build model formulas
  form_full <- as.formula(paste0(y_var, " ~ ", x_var, " + SamVolume + (1|", rand_var, ")"))
  form_null <- as.formula(paste0(y_var, " ~ SamVolume + (1|", rand_var, ")"))
  
  # Remove NA values
  data_sub <- data[!is.na(data[[y_var]]) & !is.na(data[[x_var]]), ]
  
  # Fit models
  fit_full <- glmmTMB(form_full, data = data_sub)
  fit_null <- glmmTMB(form_null, data = data_sub)
  
  # Model stats
  r2_full <- r.squaredGLMM(fit_full)
  r2_null <- r.squaredGLMM(fit_null)
  sig_test <- anova(fit_full, fit_null)
  
  r2_diff <- r2_full[1] - r2_null[1]
  pval <- sig_test[2, 8]
  lab <- if (pval < 0.001) {
    paste0('Marginal R² = ', sprintf('%.02f', r2_diff), ' p < 0.001')
  } else {
    paste0('Marginal R² = ', sprintf('%.02f', r2_diff), ' p = ', sprintf('%.03f', pval))
  }
  
  # Predictions
  data_sub$fitted_y <- predict(fit_full)
  
  # Annotation position
  x = min(data_sub[[x_var]]) + (max(data_sub[[x_var]]) - min(data_sub[[x_var]])) / 4
  y = max(data_sub[[y_var]]) - (max(data_sub[[y_var]]) - min(data_sub[[y_var]])) / 10
  
  # Plot
  p <- ggplot(data_sub, aes_string(x = x_var, y = y_var)) +
    geom_point(alpha = 0.8, color = color) +
    labs(y = y_var, x = x_var) +
    geom_smooth(aes(y = fitted_y), method = 'glm', formula = y ~ x, se = TRUE,
                lwd = 1, alpha = 0.1, color = color, fill = color) +
    annotate('text', x = x, y = y, label = lab, size = 3) +
    theme_bw() +
    theme(panel.grid = element_blank(), plot.background = element_blank())
  
  # Return list
  return(list(
    model_full = fit_full,
    model_null = fit_null,
    plot = p,
    r2_diff = r2_diff,
    p_value = pval,
    summary = summary(fit_full)
  ))
}
# Define variables
y <- "GD_Re"
x_vars <- c("GD.12S_Tre", "GD.16S_Tre", "GDvalue.COI_Tre")
randoms <- c("Primer", "BarcodeRegion")

results <- list()
param_table <- data.frame()

for (i in seq_along(x_vars)) {
  x <- x_vars[i]
  
  for (rand in randoms) {
    tag <- paste0(y, "_", x, "_", rand)
    
    res <- fit_models_plot_xy(y_var = y, x_var = x, rand_var = rand, data = data, color = useCors[4])
    
    results[[tag]] <- list(full = res$model_full,null = res$model_null,plot = res$plot)
    
    coefs <- res$summary$coefficients$cond
    
    param_table <- rbind(param_table, data.frame(
      Response = y,
      Predictor = x,
      Random = rand,
      Intercept_Est = round(coefs["(Intercept)", "Estimate"], 3),
      Intercept_SE = round(coefs["(Intercept)", "Std. Error"], 3),
      Intercept_p = round(coefs["(Intercept)", "Pr(>|z|)"], 4),
      X_Est = round(coefs[x, "Estimate"], 3),
      X_SE = round(coefs[x, "Std. Error"], 3),
      X_p = round(coefs[x, "Pr(>|z|)"], 4),
      SamVol_Est = round(coefs["SamVolume", "Estimate"], 3),
      SamVol_SE = round(coefs["SamVolume", "Std. Error"], 3),
      SamVol_p = round(coefs["SamVolume", "Pr(>|z|)"], 4),
      Marginal_R2 = round(r.squaredGLMM(res$model_full)[1], 3),
      Conditional_R2 = round(r.squaredGLMM(res$model_full)[2], 3),
      AIC = AIC(res$model_full),
      LogLik = logLik(res$model_full)[1],
      P_value = ifelse(res$p_value < 0.001, "<0.001", round(res$p_value, 3)),
      N_obs = nobs(res$model_full)
    ))
  }
}
saveRDS(results, "../result/Res_01_Validate_GD_alpha.RDS")
write.csv(param_table, "../result/Res_01_Validate_GD_alpha.csv", row.names = FALSE)
#############################################################
####***combine plots and tables***
#############################################################
##Figure 02
TD <- readRDS('../result/Res_01_Validate_TD_all.RDS')
PD <- readRDS('../result/Res_01_Validate_PD_alpha.RDS')
GD <- readRDS('../result/Res_01_Validate_GD_alpha.RDS')

p1 <- TD$`SR_Re|Primer`[[3]]+
  xlab('Basin Species Records')+ylab('Reanalyzed Species Richness') +
  ggtitle('a')

p2 <- PD$`MNTD_Re|BarcodeRegion`[[3]]+
  xlab('Tree-based MNTD (Log10)') + ylab('ASV-based MNTD (Log10)')+
  ggtitle('b')
p2

p3 <- GD$GD_Re_GDvalue.COI_Tre_BarcodeRegion[[3]]+
  xlab('Barcode-based Genetic Diversity (Log10)') + ylab('ASV-based Genetic Diversity (Log10)')+
  ggtitle('c')
p3
g <- ggarrange(p1, p2, p3, ncol = 3, align = 'hv')
g

pdf('../image/Figure_02_validate_BioDiv.pdf', width = 15, height = 4)
g
dev.off()
##Figure SI XX for main figure
TD <- readRDS('../result/Res_01_Validate_TD_all.RDS')
PD <- readRDS('../result/Res_01_Validate_PD_gamma.RDS')
GD <- readRDS('../result/Res_01_Validate_GD_gamma.RDS')

p1 <- TD$`SR_Ori|Primer`[[3]]+
  xlab('Basin Species Records')+ylab('Reanalyzed Species Richness') +
  ggtitle('a')

p2 <- PD$`MNTD_Re|BarcodeRegion`[[3]]+
  xlab('Tree-based MNTD (Log10)') + ylab('ASV-based MNTD (Log10)')+
  ggtitle('b')
p2

p3 <- GD$GD_Re_GDvalue.COI_Tre_BarcodeRegion[[3]]+
  xlab('Barcode-based Genetic Diversity (Log10)') + ylab('ASV-based Genetic Diversity (Log10)')+
  ggtitle('c')
p3
g <- ggarrange(p1, p2, p3, ncol = 3, align = 'hv')
g

pdf('../image/Figure_SI_XX_validate_BioDiv.pdf', width = 15, height = 4)
g
dev.off()
##Figure SI XX for main figure
PD1 <- readRDS('../result/Res_01_Validate_PD_alpha.RDS')
PD2 <- readRDS('../result/Res_01_Validate_PD_gamma.RDS')

p1 <- PD1$`PD_Re|BarcodeRegion`[[3]]+
  xlab('Tree-based PD (Log10)') + ylab('ASV-based PD (Log10)')+
  ggtitle('a')

p2 <- PD1$`MPD_Re|BarcodeRegion`[[3]]+
  xlab('Tree-based MPD (Log10)') + ylab('ASV-based MPD (Log10)')+
  ggtitle('b')

p3 <- PD2$`PD_Re|BarcodeRegion`[[3]]+
  xlab('Tree-based PD (Log10)') + ylab('ASV-based PD (Log10)')+
  ggtitle('c')

p4 <- PD2$`MPD_Re|BarcodeRegion`[[3]]+
  xlab('Tree-based MPD (Log10)') + ylab('ASV-based MPD (Log10)')+
  ggtitle('d')

g <- ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, align = 'hv')
g

pdf('../image/Figure_SI_XX_validate_BioDiv_4additional.pdf', width = 10, height = 8)
g
dev.off()

##Figure SI XX for main figure
GD1 <- readRDS('../result/Res_01_Validate_GD_gamma.RDS')
GD2 <- readRDS('../result/Res_01_Validate_GD_alpha.RDS')

p1 <- GD1$GD_Re_GD.12S_Tre_BarcodeRegion[[3]]+
  xlab('12S-based Genetic Diversity (Log10)') + ylab('ASV-based Genetic Diversity (Log10)')+
  ggtitle('a')

p2 <- GD1$GD_Re_GD.16S_Tre_BarcodeRegion[[3]]+
  xlab('16S-based Genetic Diversity (Log10)') + ylab('ASV-based Genetic Diversity (Log10)')+
  ggtitle('b')

p3 <- GD2$GD_Re_GD.12S_Tre_BarcodeRegion[[3]]+
  xlab('12S-based Genetic Diversity (Log10)') + ylab('ASV-based Genetic Diversity (Log10)')+
  ggtitle('c')

p4 <- GD2$GD_Re_GD.16S_Tre_BarcodeRegion[[3]]+
  xlab('16S-based Genetic Diversity (Log10)') + ylab('ASV-based Genetic Diversity (Log10)')+
  ggtitle('d')
g <- ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, align = 'hv')
g

pdf('../image/Figure_SI_XX_validate_BioDiv_GD.pdf', width = 10, height = 8)
g
dev.off()
##Table S01
TD <- read.csv('../result/Res_01_Validate_TD_all.csv')
TD$Predictor <- 'SpePool'; colnames(TD)[6:8] <- c("X_Est","X_SE","X_p"); TD$type = 'gamma'
PD1 <- read.csv('../result/Res_01_Validate_PD_gamma.csv'); PD1$type = 'gamma'
PD2 <- read.csv('../result/Res_01_Validate_PD_alpha.csv'); PD2$type = 'alpha'
GD1 <- read.csv('../result/Res_01_Validate_GD_gamma.csv'); GD1$type = 'gamma'
GD2 <- read.csv('../result/Res_01_Validate_GD_alpha.csv'); GD2$type = 'alpha'

keep_cols <- c("type","Response","Predictor","Random","X_Est", "X_SE", "X_p", "Marginal_R2", "Conditional_R2", "N_obs")

data <- rbind.data.frame(TD[keep_cols], PD1[keep_cols], PD2[keep_cols], GD1[keep_cols], GD2[keep_cols])
write.table(data, '../result/Res_01_Validate_BioDiv_model_summary.csv', row.names = F, sep = ',', quote = F)
