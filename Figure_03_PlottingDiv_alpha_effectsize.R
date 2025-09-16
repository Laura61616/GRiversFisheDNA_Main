# READNE: summary models, 27.03.2025
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
####Figure S: R2 of full models
############################################################
data <- read.csv('../result/Res_02_AllSites_BD_MultiFactor_modelFit_HG_aics_all.csv')
data$BD2 <- gsub('\\..*','',gsub('_.*','',data$BD)); data[data$BD2=='GDvalue','BD2'] <- 'GD'
data$type <- case_when(grepl('_Re',data$BD)~'reclassified',TRUE~'traditional')
data$random <- gsub(':.*','',data$model)
data$model2 <- gsub('.*:','',data$model)

data <- subset(data, BD %in% all_BioDiv_para & model2 == 'basic*(temp+human)')
data$BD <- factor(data$BD, levels = all_BioDiv_para)
data$BD2 <- factor(data$BD2, levels = (unique(data$BD2)))

dt1 <- melt(data[c(1,13:18)], id.vars = c('BD','BD2','type','random','model2'))
p1 <- ggplot(dt1, aes(x = BD, y = value, fill = variable))+
  facet_wrap(~random, nrow = 2)+
  geom_col(position = 'dodge2') +
  geom_text(aes(label = sprintf("%.2f", value)), 
            position = position_dodge2(width = 0.9, preserve = "single"),
            vjust = -0.5, size = 3) +
  scale_x_discrete(labels = fac_names)+
  xlab('') + ylab('R2 of GLMM models') +
  theme_bw() +
  theme(legend.position = 'right',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = 'grey'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10, face = 'bold'),
        legend.text = element_text(size = 10),
        legend.title = element_blank())
p1

p2 <- ggplot(data, aes(x = BD, y = AIC, fill = random))+
  geom_col(position = 'dodge2') +
  scale_x_discrete(labels = fac_names)+
  xlab('') + ylab('AIC of GLMM models') +
  theme_bw() +
  theme(legend.position = 'right',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = 'grey'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10, face = 'bold'),
        legend.text = element_text(size = 10),
        legend.title = element_blank())
p2

pdf('../image/Figure_SI_XX_MultiFactor_modelFit_HG_aics_all.pdf', width = 8, height = 10)
p1;p2
dev.off()
############################################################
####A02: ovserved and fitted
############################################################
res <- read.csv('../result/Res_02_AllSites_BD_MultiFactor_modelFit_HG_fitted_all.csv')

data <- subset(res, BD %in% all_BioDiv_para)
data <- subset(data, (BD2 %in% c('SR','Redundancy') & random == 'Primer')|(BD2 %in% c('MNTD','GD') & random == 'Barcode'))
data$BD <- factor(data$BD, levels = all_BioDiv_para, labels = fac_names)
data$BD2 <- factor(data$BD2, levels = (unique(data$BD2)))

p <- ggplot(data, aes(x = observed, y = fitted))+
  facet_wrap(~BD, scale = 'free', nrow = 4, dir = 'v')+
  geom_point(color = 'grey') +
  geom_smooth(color = 'black', fill = 'grey60', method = 'glm', formula = y ~ x, 
              level = 0.95, lwd = 1, alpha = 0.1)+
  stat_poly_eq(aes(label =  paste(after_stat(p.value.label), after_stat(adj.rr.label), sep = "*\" \"*")),
               formula = y ~ x, parse = TRUE,label.y = "top", label.x = "left", p.digits = 3, rsquared.conf.level = NA)+
  ylab('Fitted') + xlab('Obseved') +
  theme_bw() +
  theme(legend.position = 'none',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = 'grey'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10, face = 'bold'),
        legend.text = element_text(size = 10),
        legend.title = element_blank())
p

pdf('../image/Figure_SI_XX_MultiFactor_modelFit_HG_fitted.pdf', width = 8, height = 12)
p
dev.off()
############################################################
####A02: effectsize of full models
############################################################
res <- read.csv('../result/Res_02_AllSites_BD_MultiFactor_modelFit_HG_effectsize_all.csv')

factors <- c("catch_area","bio01","HFP_Lights2009","catch_area:HFP_Lights2009","catch_area:bio01","catch_area:autocor","autocor")
fac_names <- c("Catchment Area", "Temperature", "Nighttime Lights", "Catchment Area\n×Nighttime Lights", "Catchment Area\n×Temperature", "Catchment Area\n×Spatial autocor", "Spatial autocor")
bio_factors <- all_BioDiv_para
data <- subset(res, model2 == "basic*(temp+human)" & BD %in% bio_factors & factor %in% factors)
data <- subset(data, (BD2 %in% c('SR','Redundancy') & random == 'Primer')|(BD2 %in% c('MNTD','GD') & random == 'Barcode'))
data$note <- case_when(data$p < 0.05 ~ 1,data$p >= 0.05 ~ 2)
data$factor <- factor(data$factor, levels = rev(factors))
data$BD <- factor(data$BD, levels = bio_factors)
##plotting
model_nudges <- c(SR = 0.1, Redundancy = 0.05, MNTD = -0.05, GD = -0.1)
data$nudge_y = model_nudges[as.character(data$BD2)]
data$BD2 <- factor(data$BD2,
                   levels = c('SR','Redundancy','MNTD','GD'),
                   labels = c('Species Richness', 'Functional Redundancy', 'MNTD', 'Genetic Sequence Diversity'))
data$type <- factor(data$type,
                   levels = c('reclassified', 'traditional'),
                   labels = c('Reanalyzed','Original/Conventional'))
data$fill_group <- ifelse(data$note == 1, as.character(data$BD2),'none')
data$fill_group <- factor(data$fill_group, levels = c(levels(data$BD2), 'none'))
p <- ggplot(data, aes(x = Estimate, y = factor, color = BD2))+
  facet_wrap(~type)+
  geom_errorbar(aes(xmin = Estimate - 1.96*Std, xmax = Estimate + 1.96*Std), width = 0, size = 0.5, alpha = 0.3, position = position_nudge(y = data$nudge_y)) +
  geom_errorbar(aes(xmin = Estimate - 1.64*Std, xmax = Estimate + 1.64*Std), width = 0, size = 1, alpha = 0.3, position = position_nudge(y = data$nudge_y)) +
  geom_errorbar(aes(xmin = Estimate - 1.28*Std, xmax = Estimate + 1.28*Std), width = 0, size = 2, alpha = 0.5, position = position_nudge(y = data$nudge_y)) +
  # geom_rect(aes(xmin = Estimate - 1.64 * Std, xmax = Estimate + 1.64 * Std, ymin = as.numeric(factor) - 0.02 + nudge_y, ymax = as.numeric(factor) + 0.02 + nudge_y, fill = fill_group), size = 0.1, alpha = 0.7) +
  geom_point(size = 2, position = position_nudge(y = data$nudge_y)) +
  ylab('') + xlab('Estimate') +
  scale_y_discrete(labels = rev(fac_names))+
  scale_color_manual(values = useCors) +
  # scale_fill_manual(values = c(useCors, 'transparent')) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  theme_bw() +
  theme(legend.position = 'right',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = 'grey'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10, face = 'bold'),
        legend.text = element_text(size = 10),
        legend.title = element_blank())
p
pdf('../image/Figure_03_MultiFactor_modelFit_HG_effectsize.pdf', width = 10, height = 8)
p
dev.off()