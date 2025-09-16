# READNE: 02_Modelling_alpha_HG_downsampling.R, 27.03.2025
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
####A01: running models--downsampling 70% of datasets
############################################################
data <- readRDS('Use_All_Sites_modelData_BioDiv.RDS')
env <- data$env; BD <- data$BD; rownames(BD) <- BD$samID
env_factors <- c("catch_area","bio01","HFP_Lights2009","autocor",
                 "Primer","BarcodeRegion","BarcodeLength","SamVolume")
bio_factors <- intersect(all_BioDiv_para, colnames(BD))

mn <- c('Primer:basic*temp', 'Primer:basic*(temp+human)', 'Barcode:basic*temp', 'Barcode:basic*(temp+human)' )
DStimes <- 1000
sample_rate <- 0.7

pf <- '../result/Res_02_downsample/'
dir.create(pf)

for (j in 1:DStimes) {
  cat(paste0('[processing] ', j, ' out of ', DStimes, ' downsampling...\n'))
  set.seed(j)
  
  keep.ids <- sample(rownames(BD), round(nrow(BD)*sample_rate, 0))
  bio_factors <- intersect(all_BioDiv_para, colnames(BD))
  
  models <- NULL
  for (i in 1:length(bio_factors)) {
    cat('\r\t[processing] ', bio_factors[i], '...\r')
    dt <- data.frame(y = BD[keep.ids,bio_factors[i]], env[keep.ids,c('ItemID','Continent',env_factors)])
    dt <- dt[!is.na(dt$y),]
    
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
  
  out <- paste0(pf, 'Res_02_AllSites_BD_models_downsampled_', sprintf('%04d',j), '.RDS')
  saveRDS(models, out)
}