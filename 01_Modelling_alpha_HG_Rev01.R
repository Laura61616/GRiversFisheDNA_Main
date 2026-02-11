# READNE: 01_Modelling_alpha_HG.R, 27.03.2025
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

############################################################
####A01: running models
############################################################
data <- readRDS('Use_All_Sites_modelData_BioDiv_Rev01.RDS')
env <- data$env; BD <- data$BD
env[env$Primer=='Tele02','Primer'] <- 'MiFish'
env_factors <- c("catch_area","bio01","bio12","HFP_Lights2009","V_lc1",
                 "autocor","SR_init_re","SR_init_ori",
                 "Primer","BarcodeRegion","BarcodeLength","SamVolume")
bio_factors <- intersect(all_BioDiv_para, colnames(BD))

models <- aics <- anova_test <- NULL
for (i in 1:length(bio_factors)) {
    cat(paste0('[processing] ', i, ' out of ', length(bio_factors), ' bio indices: ', bio_factors[i], '...\n'))
    dt <- data.frame(y = BD[,bio_factors[i]], env[,c('ItemID','Continent',env_factors)])
    dt <- dt[!is.na(dt$y),]
    
    mn <- c('Primer:null', 
            'Primer:basic', 'Primer:basic_minus', 
            'Primer:basic*temp', 'Primer:basic+temp',
            'Primer:basic*(temp+human)','Primer:basic*temp+human',
            'Primer:basic*(temp+human_add)','Primer:basic*temp+human_add',
            'Primer:basic*(temp+human_inter)',
            'Barcode:null', 
            'Barcode:basic', 'Barcode:basic_minus', 
            'Barcode:basic*temp', 'Barcode:basic+temp',
            'Barcode:basic*(temp+human)','Barcode:basic*temp+human',
            'Barcode:basic*(temp+human_add)','Barcode:basic*temp+human_add',
            'Barcode:basic*(temp+human_inter)' )
    ##models
    b_null1 <- glmmTMB(y ~ SamVolume + autocor + (1|Primer), data = dt, family = beta_family())
    b_basic1 <- glmmTMB(y ~ catch_area * autocor + SamVolume + (1|Primer), data = dt,family = beta_family())
    b_basic1_minus <- glmmTMB(y ~ catch_area + autocor + SamVolume + (1|Primer), data = dt,family = beta_family())
    b_temp1 <- glmmTMB(y ~ catch_area * autocor + catch_area * bio01 + SamVolume + (1|Primer), data = dt,family = beta_family())
    b_temp1_minus <- glmmTMB(y ~ catch_area * autocor + bio01 + SamVolume + (1|Primer), data = dt,family = beta_family())
    b_full1 <- glmmTMB(y ~ catch_area * autocor + catch_area * bio01 + catch_area * HFP_Lights2009 + SamVolume + (1|Primer), data = dt,family = beta_family())
    b_full1_minus <- glmmTMB(y ~ catch_area * autocor + catch_area * bio01 + HFP_Lights2009 + SamVolume + (1|Primer), data = dt,family = beta_family())
    b_full1_add <- glmmTMB(y ~ catch_area * autocor + catch_area * bio01 + catch_area * HFP_Lights2009 + catch_area * V_lc1 + SamVolume + (1|Primer), data = dt,family = beta_family())
    b_full1_add_minus <- glmmTMB(y ~ catch_area * autocor + catch_area * bio01 + catch_area * HFP_Lights2009 + V_lc1 + SamVolume + (1|Primer), data = dt,family = beta_family())
    b_full1_inter <- glmmTMB(y ~ catch_area * autocor + catch_area * bio01 + catch_area * HFP_Lights2009 + catch_area * V_lc1 + HFP_Lights2009 * V_lc1 + SamVolume + (1|Primer), data = dt,family = beta_family())
    
    #--
    b_null2 <- glmmTMB(y ~ SamVolume + autocor + (1|BarcodeRegion), data = dt, family = beta_family())
    b_basic2 <- glmmTMB(y ~ catch_area * autocor + SamVolume + (1|BarcodeRegion), data = dt,family = beta_family())
    b_basic2_minus <- glmmTMB(y ~ catch_area + autocor + SamVolume + (1|BarcodeRegion), data = dt,family = beta_family())
    b_temp2 <- glmmTMB(y ~ catch_area * autocor + catch_area * bio01 + SamVolume + (1|BarcodeRegion), data = dt,family = beta_family())
    b_temp2_minus <- glmmTMB(y ~ catch_area * autocor + bio01 + SamVolume + (1|BarcodeRegion), data = dt,family = beta_family())
    b_full2 <- glmmTMB(y ~ catch_area * autocor + catch_area * bio01 + catch_area * HFP_Lights2009 + SamVolume + (1|BarcodeRegion), data = dt,family = beta_family())
    b_full2_minus <- glmmTMB(y ~ catch_area * autocor + catch_area * bio01 + HFP_Lights2009 + SamVolume + (1|BarcodeRegion), data = dt,family = beta_family())
    b_full2_add <- glmmTMB(y ~ catch_area * autocor + catch_area * bio01 + catch_area * HFP_Lights2009 + catch_area * V_lc1 + SamVolume + (1|BarcodeRegion), data = dt,family = beta_family())
    b_full2_add_minus <- glmmTMB(y ~ catch_area * autocor + catch_area * bio01 + catch_area * HFP_Lights2009 + V_lc1 + SamVolume + (1|BarcodeRegion), data = dt,family = beta_family())
    b_full2_inter <- glmmTMB(y ~ catch_area * autocor + catch_area * bio01 + catch_area * HFP_Lights2009 + catch_area * V_lc1 + HFP_Lights2009 * V_lc1 + SamVolume + (1|BarcodeRegion), data = dt,family = beta_family())
    
    ##
    models_add <- list(b_null1,b_basic1,b_basic1_minus,b_temp1,b_temp1_minus,b_full1,b_full1_minus,b_full1_add,b_full1_add_minus,b_full1_inter,
                       b_null2,b_basic2,b_basic2_minus,b_temp2,b_temp2_minus,b_full2,b_full2_minus,b_full2_add,b_full2_add_minus,b_full2_inter )
    names(models_add) <- mn
    models$add <- models_add; names(models)[length(models)] <- bio_factors[i]
    
    ##aics
    aics_add <- cbind.data.frame(BD = bio_factors[i],
                                 model = mn,
                                 convergence = sapply(models_add, performance::check_convergence),
                                 singularity = sapply(models_add, performance::check_singularity),
                                 n_obs = sapply(models_add, nobs),
                                 n_fixef = sapply(models_add, function(x) length(fixef(x))),
                                 ratio_fixef = sapply(models_add, function(x) length(fixef(x)))/sapply(models_add, nobs),
                                 do.call(rbind, lapply(models_add, function(m) summary(m)$AICtab)),
                                 do.call(rbind, lapply(models_add, MuMIn::r.squaredGLMM)) )
    
    colnames(aics_add) <- c('BD','model','convergence','sigularity','n_obs','n_fixef','ratio_fixef','AIC','BIC','logLik','deviance','df.resid','R2m','R2c')
    rownames(aics_add) <- NULL
    aics <- rbind.data.frame(aics, aics_add)
    
    ##anova
    tmp <- list(as.data.frame(anova(b_null1, b_basic1_minus)), #test plus 'catch_area'
                as.data.frame(anova(b_basic1, b_basic1_minus)), #test plus 'catch_area:autocor'
                as.data.frame(anova(b_temp1_minus, b_basic1)), # test 'bio01'
                as.data.frame(anova(b_temp1, b_temp1_minus)), # test 'catch_area:bio01'
                as.data.frame(anova(b_temp1, b_full1_minus)), # test human
                as.data.frame(anova(b_full1, b_full1_minus)), # test catch_area:human
                as.data.frame(anova(b_full1_add_minus, b_full1)), # test %crop
                as.data.frame(anova(b_full1_add, b_full1_add_minus)), # test %crop*area
                as.data.frame(anova(b_full1_inter, b_full1_add)), # test %crop*Light
                
                as.data.frame(anova(b_null2, b_basic2_minus)), #test plus 'catch_area'
                as.data.frame(anova(b_basic2, b_basic2_minus)), #test plus 'catch_area:autocor'
                as.data.frame(anova(b_temp2_minus, b_basic2)), # test 'bio01'
                as.data.frame(anova(b_temp2, b_temp2_minus)), # test 'catch_area:bio01'
                as.data.frame(anova(b_temp2, b_full2_minus)), # test human
                as.data.frame(anova(b_full2, b_full2_minus)), # test catch_area:human
                as.data.frame(anova(b_full2_add_minus, b_full2)), # test %crop
                as.data.frame(anova(b_full2_add, b_full2_add_minus)), # test %crop*area
                as.data.frame(anova(b_full2_inter, b_full2_add)) # test %crop*Light
    )
    
    anova_test_add <- data.frame(
      BD = bio_factors[i],
      random = rep(c('Primer','Barcode'), each = 9),
      test_term = rep(c('catch_area','catch_area:autocor','bio01','catch_area:bio01','NTL','catch_area:NTL', '%crop', 'catch_area:%crop', 'NTL:%crop'),2),
      mod1 = sapply(tmp, function(x) rownames(x)[1]),
      mod2 = sapply(tmp, function(x) rownames(x)[2]),
      do.call(rbind, lapply(tmp, function(x) x[2,])) )
    rownames(anova_test_add) <- NULL
    anova_test <- rbind.data.frame(anova_test, anova_test_add)
  }
write.table(aics, '../result/Res_02_AllSites_BD_MultiFactor_modelFit_HG_aics_all.csv', row.names = F, sep = ',', quote = F)
write.table(anova_test, '../result/Res_02_AllSites_BD_MultiFactor_modelFit_HG_anova_all.csv', row.names = F, sep = ',', quote = F)
saveRDS(models, '../result/Res_02_AllSites_BD_MultiFactor_modelFit_HG_models_all.RDS')
############################################################
####A01: counting n
############################################################
data <- readRDS('Use_All_Sites_modelData_BioDiv_Rev01.RDS')
BD <- data$BD

bio_factors <- intersect(all_BioDiv_para, colnames(BD))

res <- NULL
for (i in 1:length(bio_factors)) {
  cat(paste0('[processing] ', i, ' out of ', length(bio_factors), ' bio indices: ', bio_factors[i], '...\n'))
  dt <- data.frame(y = BD[,bio_factors[i]], BD[,c('ItemID','RivID','SiteID')])
  dt <- dt[!is.na(dt$y),]
  
  add <- data.frame(BD = bio_factors[i], n = nrow(dt), nRiv = length(unique(dt$RivID)), nSite = nrow(unique(dt[3:4])))
  res <- rbind.data.frame(res, add)
}
write.table(res, '../result/Res_02_AllSites_BD_MultiFactor_modelFit_HG_n_all.csv', row.names = F, sep = ',', quote = F)
############################################################
####A02: effectsize of full models
############################################################
models <- readRDS('../result/Res_02_AllSites_BD_MultiFactor_modelFit_HG_models_all.RDS')

bio_factors <- names(models)
res <- NULL
for (i in 1:length(bio_factors)) {
  cat(paste0('[processing] ', i, ' out of ', length(bio_factors), ' bio indices: ', bio_factors[i], '...\n'))
  mn <- names(models[[i]])
  for (m in mn) {
    mod <- summary(models[[i]][[m]])
    add <- data.frame(BD = bio_factors[i],
                      model = m,
                      factor = rownames(mod$coefficients$cond),
                      mod$coefficients$cond)
    add <- add[add$factor!="(Intercept)",]; rownames(add) <- NULL
    colnames(add)[5:7] <- c('Std','z','p')
    
    res <- rbind.data.frame(res, add)
  }
}

res$BD2 <- gsub('\\..*','',gsub('_.*','',res$BD)); res[res$BD2=='GDvalue','BD2'] <- 'GD'
res$type <- case_when(grepl('_Re',res$BD)~'reclassified',TRUE~'traditional')
res$random <- gsub(':.*','',res$model)
res$model2 <- gsub('.*:','',res$model)
write.table(res, '../result/Res_02_AllSites_BD_MultiFactor_modelFit_HG_effectsize_all.csv', row.names = F, sep = ',', quote = F)