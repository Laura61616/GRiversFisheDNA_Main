#############################################################
# loading libraries
#############################################################
## reading and formatting data
library(dplyr)
library(reshape2)
library(purrr)
library(broom)
library(taxize)
## plotting
library(ggplot2)
library(ggsci)
library(ggpubr)
library(patchwork)
library(ggpmisc)
library(ggspatial)
library(ggtext)
library(hrbrthemes)
## geo mapping
library(sf)
## biodiversity
library(vegan)
library(ape)
## modelling
library(glmmTMB)
library(MuMIn)
library(performance)
library(lmodel2)
library(emmeans)
library(mgcv)
library(report)
library(partR2)

library(corrplot)
#############################################################
# loading functions
#############################################################
#########basic functions
standardize_pair <- function(df) {
  df_re <- df %>%
    mutate(across(c(Var1, Var2), as.character)) %>%
    rowwise() %>%
    mutate(
      Var1_tmp = min(Var1, Var2),
      Var2_tmp = max(Var1, Var2)
    ) %>%
    ungroup()
  
  df_re <- df_re[, !(names(df_re) %in% c("Var1", "Var2"))]
  names(df_re)[names(df_re) == "Var1_tmp"] <- "Var1"
  names(df_re)[names(df_re) == "Var2_tmp"] <- "Var2"
  df_re <- df_re[,colnames(df)]
  return(df_re)
}

summary_para <- function(x, lab = 'NULL'){
  require(moments)
  if (length(x) > 40000) {
    t <- agostino.test(sample(x,40000))
  }else{
    t <- agostino.test(x)
  }

  add <- data.frame(factor = lab,
                    integer = is.integer(x),
                    Min = min(x),
                    Qu1 = quantile(x, 0.25),
                    Median = median(x),
                    Mean = mean(x),
                    Qu3 = quantile(x,0.75),
                    Max = max(x),
                    Var = var(x),
                    skewness = t$statistic[1],
                    skew_p = t$p.value,
                    dispersion_ratio = var(x)/mean(x),
                    zero_ratio = sum(x == 0)/length(x))
  rownames(add) <- NULL; return(add)

}

scaleme <- function(x, nm) {
  if (!exists('allscalings')) stop('Have to load allscalings')
  if (!(nm %in% allscalings[, 'factor'])) stop('nm not found in allscalings')
  
  current <- allscalings[allscalings$factor == nm, ]
  trans <- current$trans
  plus <- current$plus
  center <- current$center
  scale <- current$scale
  
  if (trans == 'log10') {
    x.sc <- (log10(x + plus) - center) / scale
  } else if (trans == 'sqrt') {
    x.sc <- (sqrt(x + plus) - center) / scale
  } else {
    x.sc <- (x + plus - center) / scale
  }
  return(x.sc)
}
unscaleme <- function(x.sc, nm) {
  if (!exists('allscalings')) stop('Have to load allscalings')
  if (!(nm %in% allscalings[, 'factor'])) stop('nm not found in allscalings')
  
  current <- allscalings[allscalings$factor == nm, ]
  trans <- current$trans
  plus <- current$plus
  center <- current$center
  scale <- current$scale
  
  if (trans == 'log10') {
    x <- 10^(x.sc * scale + center) - plus
  } else if (trans == 'sqrt') {
    x <- (x.sc * scale + center)^2 - plus
  } else {
    x <- x.sc * scale + center - plus
  }
  return(x)
}

sqrt_transform <- function(x) {
  ifelse(x < 0, -sqrt(abs(x)), sqrt(x))
}

shrink_01 <- function(x, eps = 1e-6) {
  x * (1 - 2 * eps) + eps
}

corvif <- function(dataz) {
  dataz <- as.data.frame(dataz)
  
  #vif part
  form    <- formula(paste("fooy ~ ",paste(strsplit(names(dataz)," "),collapse=" + ")))
  dataz   <- data.frame(fooy=1 + rnorm(nrow(dataz)) ,dataz)
  lm_mod  <- lm(form,dataz)
  
  cat("\n\nVariance inflation factors\n\n")
  print(myvif(lm_mod))
}
myvif <- function(mod) {
  v <- vcov(mod)
  assign <- attributes(model.matrix(mod))$assign
  tmp_cor <- matrix()
  if (names(coefficients(mod)[1]) == "(Intercept)") {
    v <- v[-1, -1]
    assign <- assign[-1]
  } else warning("No intercept: vifs may not be sensible.")
  terms <- labels(terms(mod))
  n.terms <- length(terms)
  if (n.terms < 2) stop("The model contains fewer than 2 terms")
  if (length(assign) > dim(v)[1] ) {
    diag(tmp_cor)<-0
    if (any(tmp_cor==1.0)){
      return("Sample size is too small, 100% collinearity is present")
    } else {
      return("Sample size is too small")
    }
  }
  R <- cov2cor(v)
  detR <- det(R)
  result <- matrix(0, n.terms, 3)
  rownames(result) <- terms
  colnames(result) <- c("GVIF", "Df", "GVIF^(1/2Df)")
  for (term in 1:n.terms) {
    subs <- which(assign == term)
    result[term, 1] <- det(as.matrix(R[subs, subs])) * det(as.matrix(R[-subs, -subs])) / detR
    result[term, 2] <- length(subs)
  }
  if (all(result[, 2] == 1)) {
    result <- data.frame(GVIF=result[, 1])
  } else {
    result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
  }
  invisible(result)
}
#########spatial autocorrelation
spatial_autocor <- function(lon, lat){
  require(geosphere)
  if (length(lon) == 1) {
    return(0)
  }
  df <- data.frame(lon = lon, lat = lat)
  dist_matrix <- outer(
    1:nrow(df), 1:nrow(df),
    Vectorize(function(i, j) distVincentySphere(df[i, c("lon", "lat")], df[j, c("lon", "lat")])))
  return(mean(apply(dist_matrix, 1, function(row) min(row[row > 0])), na.rm = T))
}
#########generating label from model
gam_fit_lab2 <- function(fit, fac){
  p <- fit$p.table[fac,4]
  r <- fit$r.sq

  ind <- case_when(p < 0.001 ~ '***',
                   (p >= 0.001)&(p < 0.01) ~ '**', 
                   (p >= 0.01)&(p < 0.05) ~ '*', 
                   (p >= 0.05)&(p < 0.1) ~ '.', 
                   p >= 0.1 ~ '')
  if (is.nan(p)) {
    lab <- paste0('R2 = ', round(r,2), 
                  ', p = NA')
    return(lab)
    
  }
  if (p < 0.001) {
    lab <- paste0('R2 = ', round(r,2), 
                  ', p < 0.001', ' ',ind)
  }else{
    lab <- paste0('R2 = ', round(r,2), 
                  ', p = ', round(p,3), ' ',ind)
  }
  return(lab)
}
lme_fit_lab <- function(fit, fac, explained = NULL){
  p <- fit$coefficients[fac,'Pr(>|t|)']
  r <- fit$AICtab
  if (is.null(explained)) {
    expl <- NULL
  }else{
    expl <- explained
  }
  ind <- case_when(p < 0.001 ~ '***',
                   (p >= 0.001)&(p < 0.01) ~ '**', 
                   (p >= 0.01)&(p < 0.05) ~ '*', 
                   (p >= 0.05)&(p < 0.1) ~ '.', 
                   p >= 0.1 ~ '')
  if (p < 0.001) {
    lab <- paste0('AIC = ', round(r,2), 
                  '\nexplained = ', round(expl*100, 2),
                  '%\np < 0.001', ' ',ind)
  }else{
    lab <- paste0('AIC = ', round(r,2), 
                  '\nexplained = ', round(expl*100, 2),
                  '%\np = ', round(p,3), ' ',ind)
  }
  return(lab)
}
glme_fit_lab <- function(fit, fac, explained = NULL){
  p <- fit$coefficients$cond[fac,'Pr(>|z|)']
  r <- fit$AICtab[1]
  if (is.null(explained)) {
    expl <- NULL
  }else{
    expl <- explained
  }
  ind <- case_when(p < 0.001 ~ '***',
                   (p >= 0.001)&(p < 0.01) ~ '**', 
                   (p >= 0.01)&(p < 0.05) ~ '*', 
                   (p >= 0.05)&(p < 0.1) ~ '.', 
                   p >= 0.1 ~ '')
  if (p < 0.001) {
    lab <- paste0('AIC = ', round(r,2), 
                  '\nexplained = ', round(expl*100, 2),
                  '%\np < 0.001', ' ',ind)
  }else{
    lab <- paste0('AIC = ', round(r,2), 
                  '\nexplained = ', round(expl*100, 2),
                  '%\np = ', round(p,3), ' ',ind)
  }
  return(lab)
}
summary_lme <- function(rep, r, useFactors, other_para = 'Continent'){
  est <- as.data.frame(r$Ests)
  rownames(est) <- est$term
  est <- est[c(useFactors, grep(other_para, rownames(est), value = T)),]
  
  explained <- as.data.frame(r$R2)
  rownames(explained) <- explained$term
  explained <- explained[c(useFactors, other_para),]
  
  rep <- rep[!is.na(rep$Parameter),]
  rownames(rep) <- gsub('_×_',':',gsub(' ','_',rep$Parameter))
  cond_r2 <- na.omit(rep[rep$Parameter=='R2 (conditional)','Fit'])
  if (length(cond_r2) == 0) {
    cond_r2 <- NA
  }
  
  paras <- c(useFactors, grep(other_para,rownames(rep),value = T))
  paras[grep('Continent',paras)[1]] <- 'Continent'
  full_res <- data.frame(
    'AIC' = rep[rep$Parameter=='AIC','Fit'],
    'cond_r2' = cond_r2,
    'marg_r2' = rep[rep$Parameter=='R2 (marginal)','Fit'],
    rep[c(useFactors, grep(other_para,rownames(rep),value = T)),],
    'estimate' = est[, 'estimate'],
    'expained' = explained[paras, 'estimate']
  )
  return(full_res)
}
