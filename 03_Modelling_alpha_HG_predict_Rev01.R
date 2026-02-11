# READNE: 03_Modelling_alpha_HG_predict.R, 27.03.2025
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
####A00: range of new data
############################################################
range_w <- read.csv('Use_All_Sites_modelData_BioDiv_newdata_range_Rev01.csv', row.names = 1)
allscalings <- read.csv('Use_All_Sites_modelData_BioDiv_transform_scaling_Rev01.csv')

temps <- seq(ceiling(range_w['bio01','Min']), floor(range_w['bio01','Max']), by = 1)
temps_sc <- scaleme(temps, 'bio01')

pres <- seq(ceiling(range_w['bio12','Min']), floor(range_w['bio12','Max']), by = 100)
pres_sc <- scaleme(pres, 'bio12')

catch_areas <- exp(seq(log(range_w['catch_area','Min']), log(range_w['catch_area','Max']), length.out = 10))
catch_areas_sc <- scaleme(catch_areas, 'catch_area')

humans <- seq(ceiling(range_w['HFP_Lights2009','Min']), floor(range_w['HFP_Lights2009','Max']), by = 1)
humans_sc <- scaleme(humans, 'HFP_Lights2009')

humans2 <- 10^(seq(log10(range_w['V_lc1','Min']+1e-3), log10(range_w['V_lc1','Max']+1e-3), by = .5))-1e-3
humans2_sc <- scaleme(humans2, 'V_lc1')

newdata <- list(bio01 = temps_sc, bio12 = pres_sc, catch_areas = catch_areas_sc, human = humans_sc, human2 = humans2_sc)
saveRDS(newdata, 'Use_All_Sites_modelData_BioDiv_newdata_Rev01.RDS')
############################################################
####A01: prediction of full models
############################################################
allscalings <- read.csv('Use_All_Sites_modelData_BioDiv_transform_scaling_Rev01.csv')
newdata <- readRDS('Use_All_Sites_modelData_BioDiv_newdata_Rev01.RDS')
models <- readRDS('../result/Res_02_AllSites_BD_MultiFactor_modelFit_HG_models_all.RDS')

bio_factors <- names(models)

res <- NULL
for (i in 1:length(bio_factors)) {
  cat(paste0('[processing] ', i, ' out of ', length(bio_factors), ' bio indices: ', bio_factors[i], '...\n'))
  
  ##Primer as random
  mod <- models[[i]]$`Primer:basic*temp+human_add`
  dt <- mod$frame
  
  t1 <- emtrends(mod, ~bio01, var = 'catch_area', type = "response", at = list(bio01=newdata$bio01))
  t2 <- emtrends(mod, ~HFP_Lights2009, var = 'catch_area', type = "response", at = list(HFP_Lights2009=newdata$human))
  t3 <- emtrends(mod, ~catch_area, var = 'HFP_Lights2009', type = "response", at = list(catch_area=newdata$catch_areas))
  
  t1_95 <- summary(t1); t1_99 <- summary(t1, level = 0.99)
  tt1 <- data.frame(t1_95, t1_99[,5:6]); colnames(tt1)[5:8] <- c('LCL_95','UCL_95','LCL_99','UCL_99')
  t2_95 <- summary(t2); t2_99 <- summary(t2, level = 0.99)
  tt2 <- data.frame(t2_95, t2_99[,5:6]); colnames(tt2)[5:8] <- c('LCL_95','UCL_95','LCL_99','UCL_99')
  t3_95 <- summary(t3); t3_99 <- summary(t3, level = 0.99)
  tt3 <- data.frame(t3_95, t3_99[,5:6]); colnames(tt3)[5:8] <- c('LCL_95','UCL_95','LCL_99','UCL_99')
  
  ##Barcode as random
  mod <- models[[i]]$`Barcode:basic*temp+human_add`
  dt <- mod$frame
  
  t4 <- emtrends(mod, ~bio01, var = 'catch_area', type = "response", at = list(bio01=newdata$bio01))
  t5 <- emtrends(mod, ~HFP_Lights2009, var = 'catch_area', type = "response", at = list(HFP_Lights2009=newdata$human))
  t6 <- emtrends(mod, ~catch_area, var = 'HFP_Lights2009', type = "response", at = list(catch_area=newdata$catch_areas))
  
  t4_95 <- summary(t4); t4_99 <- summary(t4, level = 0.99)
  tt4 <- data.frame(t4_95, t4_99[,5:6]); colnames(tt4)[5:8] <- c('LCL_95','UCL_95','LCL_99','UCL_99')
  t5_95 <- summary(t5); t5_99 <- summary(t5, level = 0.99)
  tt5 <- data.frame(t5_95, t5_99[,5:6]); colnames(tt5)[5:8] <- c('LCL_95','UCL_95','LCL_99','UCL_99')
  t6_95 <- summary(t6); t6_99 <- summary(t6, level = 0.99)
  tt6 <- data.frame(t6_95, t6_99[,5:6]); colnames(tt6)[5:8] <- c('LCL_95','UCL_95','LCL_99','UCL_99')
  
  ##archive
  add <- data.frame(BD = bio_factors[i],
                    random = c(rep('Primer', nrow(tt1)+nrow(tt2)+nrow(tt3)), rep('Barcode', nrow(tt4)+nrow(tt5)+nrow(tt6))),
                    model = 'area*temp+human_add',
                    x_name = c(rep('bio01', nrow(tt1)), rep('HFP_Lights2009', nrow(tt2)), rep('catch_area',nrow(tt3)),rep('bio01',nrow(tt4)),rep('HFP_Lights2009', nrow(tt5)), rep('catch_area',nrow(tt6)) ),
                    x_sc = c(tt1$bio01, tt2$HFP_Lights2009, tt3$catch_area, tt4$bio01, tt5$HFP_Lights2009, tt6$catch_area),
                    y_name = c(rep('catch_area.trend', nrow(tt1)+nrow(tt2)), rep('HFP_Lights2009.trend',nrow(tt3)),rep('catch_area.trend', nrow(tt4)+nrow(tt5)), rep('HFP_Lights2009.trend',nrow(tt6)) ),
                    y = c(tt1$catch_area.trend, tt2$catch_area.trend, tt3$HFP_Lights2009.trend, tt4$catch_area.trend, tt5$catch_area.trend, tt6$HFP_Lights2009.trend ),
                    rbind.data.frame(tt1[,3:8],tt2[,3:8],tt3[,3:8],tt4[,3:8],tt5[,3:8],tt6[,3:8] ) )
  
  res <- rbind.data.frame(res, add)
}

res$x <- case_when(res$x_name == 'bio01' ~ unscaleme(res$x_sc, 'bio01'),
                   res$x_name == 'HFP_Lights2009' ~ unscaleme(res$x_sc, 'HFP_Lights2009'),
                   res$x_name == 'catch_area' ~ unscaleme(res$x_sc, 'catch_area'))

res$BD2 <- gsub('\\..*','',gsub('_.*','',res$BD)); res[res$BD2=='GDvalue','BD2'] <- 'GD'
res$type <- case_when(grepl('_Re',res$BD)~'reclassified',TRUE~'traditional')
write.table(res, '../result/Res_02_SAR_based_Full.csv', row.names = F, sep = ',', quote = F)
############################################################
####A01: prediction for single ItemID/RiverID
############################################################
# prepare bio01, HFP_Light, autocor for single RiverID
data <- readRDS('Use_All_Sites_modelData_BioDiv_Rev01.RDS')
env <- data$env
env_use <- env %>%
  group_by(N_RivID) %>%
  summarise(bio01 = median(bio01),
            HFP_Lights2009 = median(HFP_Lights2009),
            autocor = median(autocor)) %>%
  ungroup()
# predicting
allscalings <- read.csv('Use_All_Sites_modelData_BioDiv_transform_scaling_Rev01.csv')
models <- readRDS('../result/Res_02_AllSites_BD_MultiFactor_modelFit_HG_models_all.RDS')

bio_factors <- names(models)

res <- NULL
for (i in 1:length(bio_factors)) {
  cat(paste0('[processing] ', i, ' out of ', length(bio_factors), ' bio indices: ', bio_factors[i], '...\n'))
  
  ##Primer as random
  mod <- models[[i]]$`Primer:basic*(temp+human)`
  dt <- mod$frame
  t1 <- do.call(rbind, lapply(seq_len(nrow(env_use)), 
                              function(i) { row <- env_use[i, ]
                              tr_95 <- emtrends(mod, ~bio01+HFP_Lights2009+autocor, var = "catch_area", type = "response", at = list(bio01 = row$bio01, HFP_Lights2009 = row$HFP_Lights2009, autocor = row$autocor), level = 0.95); 
                              tr_99 <- emtrends(mod, ~bio01+HFP_Lights2009+autocor, var = "catch_area", type = "response", at = list(bio01 = row$bio01, HFP_Lights2009 = row$HFP_Lights2009, autocor = row$autocor), level = 0.99); 
                              df_95 <- as.data.frame(tr_95); df_99 <- as.data.frame(tr_99); 
                              data.frame(RivID = row$N_RivID, 
                                         bio01 = row$bio01, HFP_Lights2009 = row$HFP_Lights2009, autocor = row$autocor, 
                                         trend = df_95$catch_area.trend, SE = df_95$SE, df = df_95$df, 
                                         LCL_95 = df_95$asymp.LCL, UCL_95 = df_95$asymp.UCL, 
                                         LCL_99 = df_99$asymp.LCL, UCL_99 = df_99$asymp.UCL, row.names = NULL) }))
  
  ##Barcode as random
  mod <- models[[i]]$`Barcode:basic*(temp+human)`
  dt <- mod$frame
  t2 <- do.call(rbind, lapply(seq_len(nrow(env_use)), 
                              function(i) { row <- env_use[i, ]
                              tr_95 <- emtrends(mod, ~bio01+HFP_Lights2009+autocor, var = "catch_area", type = "response", at = list(bio01 = row$bio01, HFP_Lights2009 = row$HFP_Lights2009, autocor = row$autocor), level = 0.95); 
                              tr_99 <- emtrends(mod, ~bio01+HFP_Lights2009+autocor, var = "catch_area", type = "response", at = list(bio01 = row$bio01, HFP_Lights2009 = row$HFP_Lights2009, autocor = row$autocor), level = 0.99); 
                              df_95 <- as.data.frame(tr_95); df_99 <- as.data.frame(tr_99); 
                              data.frame(RivID = row$N_RivID, 
                                         bio01 = row$bio01, HFP_Lights2009 = row$HFP_Lights2009, autocor = row$autocor, 
                                         trend = df_95$catch_area.trend, SE = df_95$SE, df = df_95$df, 
                                         LCL_95 = df_95$asymp.LCL, UCL_95 = df_95$asymp.UCL, 
                                         LCL_99 = df_99$asymp.LCL, UCL_99 = df_99$asymp.UCL, row.names = NULL) }))
  
  ##archive
  add <- data.frame(BD = bio_factors[i],
                    random = c(rep('Primer', nrow(t1)), rep('Barcode', nrow(t2))),
                    model = 'area*temp+human_add',
                    rbind.data.frame(t1,t2 ) )
  
  res <- rbind.data.frame(res, add)
}
res$BD2 <- gsub('\\..*','',gsub('_.*','',res$BD)); res[res$BD2=='GDvalue','BD2'] <- 'GD'
res$type <- case_when(grepl('_Re',res$BD)~'reclassified',TRUE~'traditional')
res$bio01_unscale <- unscaleme(res$bio01, 'bio01')
res$HFP_Lights2009_unscale <- unscaleme(res$HFP_Lights2009, 'HFP_Lights2009')

head(res)
write.table(res, '../result/Res_02_SAR_based_Full_perRiver.csv', row.names = F, sep = ',', quote = F)
############################################################
####A02: prediction of downsampled models
############################################################
allscalings <- read.csv('Use_All_Sites_modelData_BioDiv_transform_scaling_Rev01.csv')
newdata <- readRDS('Use_All_Sites_modelData_BioDiv_newdata_Rev01.RDS')
pf <- list.files('../result/Res_02_downsample/', full.names = T)

out <- paste0('../result/Res_02_SAR_based_Full_downsampled_',Sys.Date(),'.csv')
#--!!be careful, will remove the output while the same outdir existed!!
if (file.exists(out)) {
  file.remove(out)
}

for (i in 1:length(pf)) {
  cat(paste0('[processing] ', i, ' out of ', length(pf), ' reshuffle IDs...\n'))
  models <- readRDS(pf[i])
  
  bio_factors <- names(models)
  for (j in 1:length(models)) {
    cat(paste0('\r\t[processing] ', bio_factors[j], '...\r'))
    ##Primer as random
    mod <- models[[i]]$`Primer:basic*temp+human_add`
    dt <- mod$frame
    
    t1 <- emtrends(mod, ~bio01, var = 'catch_area', type = "response", at = list(bio01=newdata$bio01))
    t2 <- emtrends(mod, ~HFP_Lights2009, var = 'catch_area', type = "response", at = list(HFP_Lights2009=newdata$human))
    t3 <- emtrends(mod, ~catch_area, var = 'HFP_Lights2009', type = "response", at = list(catch_area=newdata$catch_areas))
    
    t1_95 <- summary(t1); t1_99 <- summary(t1, level = 0.99)
    tt1 <- data.frame(t1_95, t1_99[,5:6]); colnames(tt1)[5:8] <- c('LCL_95','UCL_95','LCL_99','UCL_99')
    t2_95 <- summary(t2); t2_99 <- summary(t2, level = 0.99)
    tt2 <- data.frame(t2_95, t2_99[,5:6]); colnames(tt2)[5:8] <- c('LCL_95','UCL_95','LCL_99','UCL_99')
    t3_95 <- summary(t3); t3_99 <- summary(t3, level = 0.99)
    tt3 <- data.frame(t3_95, t3_99[,5:6]); colnames(tt3)[5:8] <- c('LCL_95','UCL_95','LCL_99','UCL_99')
    
    ##Barcode as random
    mod <- models[[i]]$`Barcode:basic*temp+human_add`
    dt <- mod$frame
    
    t4 <- emtrends(mod, ~bio01, var = 'catch_area', type = "response", at = list(bio01=newdata$bio01))
    t5 <- emtrends(mod, ~HFP_Lights2009, var = 'catch_area', type = "response", at = list(HFP_Lights2009=newdata$human))
    t6 <- emtrends(mod, ~catch_area, var = 'HFP_Lights2009', type = "response", at = list(catch_area=newdata$catch_areas))
    
    t4_95 <- summary(t4); t4_99 <- summary(t4, level = 0.99)
    tt4 <- data.frame(t4_95, t4_99[,5:6]); colnames(tt4)[5:8] <- c('LCL_95','UCL_95','LCL_99','UCL_99')
    t5_95 <- summary(t5); t5_99 <- summary(t5, level = 0.99)
    tt5 <- data.frame(t5_95, t5_99[,5:6]); colnames(tt5)[5:8] <- c('LCL_95','UCL_95','LCL_99','UCL_99')
    t6_95 <- summary(t6); t6_99 <- summary(t6, level = 0.99)
    tt6 <- data.frame(t6_95, t6_99[,5:6]); colnames(tt6)[5:8] <- c('LCL_95','UCL_95','LCL_99','UCL_99')
    
    ##archive
    add <- data.frame(DSid = i,
                      BD = bio_factors[i],
                      random = c(rep('Primer', nrow(tt1)+nrow(tt2)+nrow(tt3)), rep('Barcode', nrow(tt4)+nrow(tt5)+nrow(tt6))),
                      model = 'area*(temp+human)',
                      x_name = c(rep('bio01', nrow(tt1)), rep('HFP_Lights2009', nrow(tt2)), rep('catch_area',nrow(tt3)),rep('bio01',nrow(tt4)),rep('HFP_Lights2009', nrow(tt5)), rep('catch_area',nrow(tt6)) ),
                      x_sc = c(tt1$bio01, tt2$HFP_Lights2009, tt3$catch_area, tt4$bio01, tt5$HFP_Lights2009, tt6$catch_area),
                      y_name = c(rep('catch_area.trend', nrow(tt1)+nrow(tt2)), rep('HFP_Lights2009.trend',nrow(tt3)),rep('catch_area.trend', nrow(tt4)+nrow(tt5)), rep('HFP_Lights2009.trend',nrow(tt6)) ),
                      y = c(tt1$catch_area.trend, tt2$catch_area.trend, tt3$HFP_Lights2009.trend, tt4$catch_area.trend, tt5$catch_area.trend, tt6$HFP_Lights2009.trend ),
                      rbind.data.frame(tt1[,3:8],tt2[,3:8],tt3[,3:8],tt4[,3:8],tt5[,3:8],tt6[,3:8] ) )
    
    add$x <- case_when(add$x_name == 'bio01' ~ unscaleme(add$x_sc, 'bio01'),
                       add$x_name == 'HFP_Lights2009' ~ unscaleme(add$x_sc, 'HFP_Lights2009'),
                       add$x_name == 'catch_area' ~ unscaleme(add$x_sc, 'catch_area'))
    
    add$BD2 <- gsub('\\..*','',gsub('_.*','',add$BD)); add[add$BD2=='GDvalue','BD2'] <- 'GD'
    add$type <- case_when(grepl('_Re',add$BD)~'reclassified',TRUE~'traditional')
    
    if (!file.exists(out)) {
      write.table(add, out, sep = ',', row.names = F, col.names = T, quote = F, append = F)
    }else{
      write.table(add, out, sep = ',', row.names = F, col.names = F, quote = F, append = T)
      
    }
  }
}
############################################################
####A02: prediction of reshuffled models
############################################################
allscalings <- read.csv('Use_All_Sites_modelData_BioDiv_transform_scaling_Rev01.csv')
newdata <- readRDS('Use_All_Sites_modelData_BioDiv_newdata_Rev01.RDS')

pf <- list.files('../result/Res_02_reshuffle/', full.names = T)
out <- paste0('../result/Res_02_SAR_based_Human_reshuffled_',Sys.Date(),'.csv')
##!!be careful, will remove the output while the same outdir existed!!
if (file.exists(out)) {
  file.remove(out)
}
for (i in 1:length(pf)) {
  cat(paste0('[processing] ', i, ' out of ', length(pf), ' reshuffle IDs...\n'))
  models <- readRDS(pf[i])
  
  bio_factors <- names(models)
  for (j in 1:length(models)) {
    cat(paste0('\r\t[processing] ', bio_factors[j], '...\r'))
    ##Primer as random
    mod <- models[[i]]$`Primer:basic*temp+human_add`
    dt <- mod$frame
    
    t1 <- emtrends(mod, ~bio01, var = 'catch_area', type = "response", at = list(bio01=newdata$bio01))
    t2 <- emtrends(mod, ~HFP_Lights2009, var = 'catch_area', type = "response", at = list(HFP_Lights2009=newdata$human))
    t3 <- emtrends(mod, ~catch_area, var = 'HFP_Lights2009', type = "response", at = list(catch_area=newdata$catch_areas))
    
    t1_95 <- summary(t1); t1_99 <- summary(t1, level = 0.99)
    tt1 <- data.frame(t1_95, t1_99[,5:6]); colnames(tt1)[5:8] <- c('LCL_95','UCL_95','LCL_99','UCL_99')
    t2_95 <- summary(t2); t2_99 <- summary(t2, level = 0.99)
    tt2 <- data.frame(t2_95, t2_99[,5:6]); colnames(tt2)[5:8] <- c('LCL_95','UCL_95','LCL_99','UCL_99')
    t3_95 <- summary(t3); t3_99 <- summary(t3, level = 0.99)
    tt3 <- data.frame(t3_95, t3_99[,5:6]); colnames(tt3)[5:8] <- c('LCL_95','UCL_95','LCL_99','UCL_99')
    
    ##Barcode as random
    mod <- models[[i]]$`Barcode:basic*temp+human_add`
    dt <- mod$frame
    
    t4 <- emtrends(mod, ~bio01, var = 'catch_area', type = "response", at = list(bio01=newdata$bio01))
    t5 <- emtrends(mod, ~HFP_Lights2009, var = 'catch_area', type = "response", at = list(HFP_Lights2009=newdata$human))
    t6 <- emtrends(mod, ~catch_area, var = 'HFP_Lights2009', type = "response", at = list(catch_area=newdata$catch_areas))
    
    t4_95 <- summary(t4); t4_99 <- summary(t4, level = 0.99)
    tt4 <- data.frame(t4_95, t4_99[,5:6]); colnames(tt4)[5:8] <- c('LCL_95','UCL_95','LCL_99','UCL_99')
    t5_95 <- summary(t5); t5_99 <- summary(t5, level = 0.99)
    tt5 <- data.frame(t5_95, t5_99[,5:6]); colnames(tt5)[5:8] <- c('LCL_95','UCL_95','LCL_99','UCL_99')
    t6_95 <- summary(t6); t6_99 <- summary(t6, level = 0.99)
    tt6 <- data.frame(t6_95, t6_99[,5:6]); colnames(tt6)[5:8] <- c('LCL_95','UCL_95','LCL_99','UCL_99')
    
    ##archive
    add <- data.frame(RSid = i,
                      BD = bio_factors[i],
                      random = c(rep('Primer', nrow(tt1)+nrow(tt2)+nrow(tt3)), rep('Barcode', nrow(tt4)+nrow(tt5)+nrow(tt6))),
                      model = 'area*temp+human_add',
                      x_name = c(rep('bio01', nrow(tt1)), rep('HFP_Lights2009', nrow(tt2)), rep('catch_area',nrow(tt3)),rep('bio01',nrow(tt4)),rep('HFP_Lights2009', nrow(tt5)), rep('catch_area',nrow(tt6)) ),
                      x_sc = c(tt1$bio01, tt2$HFP_Lights2009, tt3$catch_area, tt4$bio01, tt5$HFP_Lights2009, tt6$catch_area),
                      y_name = c(rep('catch_area.trend', nrow(tt1)+nrow(tt2)), rep('HFP_Lights2009.trend',nrow(tt3)),rep('catch_area.trend', nrow(tt4)+nrow(tt5)), rep('HFP_Lights2009.trend',nrow(tt6)) ),
                      y = c(tt1$catch_area.trend, tt2$catch_area.trend, tt3$HFP_Lights2009.trend, tt4$catch_area.trend, tt5$catch_area.trend, tt6$HFP_Lights2009.trend ),
                      rbind.data.frame(tt1[,3:8],tt2[,3:8],tt3[,3:8],tt4[,3:8],tt5[,3:8],tt6[,3:8] ) )
    
    add$x <- case_when(add$x_name == 'bio01' ~ unscaleme(add$x_sc, 'bio01'),
                       add$x_name == 'HFP_Lights2009' ~ unscaleme(add$x_sc, 'HFP_Lights2009'),
                       add$x_name == 'catch_area' ~ unscaleme(add$x_sc, 'catch_area'))
    
    add$BD2 <- gsub('\\..*','',gsub('_.*','',add$BD)); add[add$BD2=='GDvalue','BD2'] <- 'GD'
    add$type <- case_when(grepl('_Re',add$BD)~'reclassified',TRUE~'traditional')
    
    if (!file.exists(out)) {
      write.table(add, out, sep = ',', row.names = F, col.names = T, quote = F, append = F)
    }else{
      write.table(add, out, sep = ',', row.names = F, col.names = F, quote = F, append = T)
    }
  }
}

############################################################
####A03: combine data
############################################################
allscalings <- read.csv('Use_All_Sites_modelData_BioDiv_transform_scaling_Rev01.csv')

data_RS <- read.csv('../result/Res_02_SAR_based_Full_reshuffled.csv')
data_DS <- read.csv('../result/Res_02_SAR_based_Full_downsampled.csv')

dt <- rbind.data.frame(data.frame(DT = 'reshuffled', id = data_RS$DSid, data_RS[-1]),
                       data.frame(DT = 'downsampled', id = data_DS$DSid, data_DS[-1]) )
dt$x <- case_when(dt$x_name == 'bio01' ~ unscaleme(dt$x_sc, 'bio01'),
                   dt$x_name == 'HFP_Lights2009' ~ unscaleme(dt$x_sc, 'HFP_Lights2009'),
                   dt$x_name == 'catch_area' ~ unscaleme(dt$x_sc, 'catch_area'))

res <- list('SAR_bio01' = subset(dt, x_name == 'bio01'),
            'SAR_human' = subset(dt, x_name == 'HFP_Lights2009'),
            'SHR_area' = subset(dt, x_name == 'catch_area'))
saveRDS(res, '../result/Res_02_SAR_based_Full_all.RDS')
