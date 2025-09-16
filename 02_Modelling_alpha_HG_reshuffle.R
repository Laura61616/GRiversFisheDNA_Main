# READNE: 02_Modelling_alpha_HG_reshuffle.R, 27.03.2025
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
####A00: running models--reshuffling catch_area
############################################################
data <- readRDS('Use_All_Sites_modelData_BioDiv.RDS')
env <- data$env; BD <- data$BD; rownames(BD) <- BD$samID
env_factors <- c("catch_area","bio01","HFP_Lights2009","autocor",
                 "Primer","BarcodeRegion","BarcodeLength","SamVolume")
bio_factors <- intersect(all_BioDiv_para, colnames(BD))

mn <- c('Primer:basic*temp', 'Primer:basic*(temp+human)', 'Barcode:basic*temp', 'Barcode:basic*(temp+human)' )
RStimes <- 10

pf <- '../result/Res_02_reshuffle/'
dir.create(pf)
for (j in 1:RStimes) {
  cat(paste0('[processing] ', j, ' out of ', RStimes, ' reshufflings...\n'))
  set.seed(j)
  
  models <- NULL
  for (i in 1:length(bio_factors)) {
    cat('\r\t[processing] ', bio_factors[i], '...\r')
    dt_tmp <- data.frame(y = BD[,bio_factors[i]], env[,c('ItemID','Continent',env_factors)])
    dt_tmp <- dt_tmp[!is.na(dt_tmp$y),]
    dt <- dt_tmp; dt$bio01 <- sample(dt$bio01, nrow(dt)); dt$HFP_Lights2009 <- sample(dt$HFP_Lights2009, nrow(dt))
    ##models
    #--
    b_temp1 <- glmmTMB(y ~ catch_area * autocor + catch_area * bio01 + SamVolume + (1|Primer), data = dt,family = beta_family())
    b_full1 <- glmmTMB(y ~ catch_area * autocor + catch_area * bio01 + catch_area * HFP_Lights2009 + SamVolume + (1|Primer), data = dt,family = beta_family())
    
    b_temp2 <- glmmTMB(y ~ catch_area * autocor + catch_area * bio01 + SamVolume + (1|BarcodeRegion), data = dt,family = beta_family())
    b_full2 <- glmmTMB(y ~ catch_area * autocor + catch_area * bio01 + catch_area * HFP_Lights2009 + SamVolume + (1|BarcodeRegion), data = dt,family = beta_family())
    
    ##
    models_add <- list(b_temp1, b_full1, b_temp2, b_full2 )
    names(models_add) <- mn
    models$add <- models_add; names(models)[length(models)] <- bio_factors[i]
  }
  
  out <- paste0(pf, 'Res_02_AllSites_BD_models_reshuffle_', sprintf('%04d',j), '.RDS')
  saveRDS(models, out)
}