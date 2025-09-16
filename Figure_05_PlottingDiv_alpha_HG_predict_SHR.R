# READNE: plotting species-human relationship, 17.06.2025
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
####Figure_05: plotting ZSHR based on downsampling & reshuffles
############################################################
##Human verus ZSAR
res <- readRDS('../result/Res_02_SAR_based_Full_all.RDS')
res <- subset(res$SHR_area, BD %in% all_BioDiv_para & DT == "downsampled")
res <- subset(res, (BD2 %in% c('SR','Redundancy') & random == 'Primer')|(BD2 %in% c('MNTD','GD') & random == 'Barcode') )
res$BD2 <- factor(res$BD2,
                  levels = c('SR','Redundancy','MNTD','GD'),
                  labels = c('Species Richness', 'Functional Redundancy', 'MNTD', 'Genetic Diversity'))
res$type <- factor(res$type,
                   levels = c('reclassified', 'traditional'),
                   labels = c('Reanalyzed','Reported/Conventional'))

dt_avg  <- res %>%
  group_by(BD, BD2, x_name, x, x_sc, y_name, type) %>%
  summarise(
    y = mean(y, na.rm = TRUE),
    SE = mean(SE, na.rm = TRUE),
    LCL_95 = mean(LCL_95, na.rm = TRUE),
    UCL_95 = mean(UCL_95, na.rm = TRUE),
    LCL_99 = mean(LCL_99, na.rm = TRUE),
    UCL_99 = mean(UCL_99, na.rm = TRUE),
    .groups = "drop"  )

p_values <- read.csv('../result/Res_02_SHR_based_Full_AreaEffects_significance.csv')
p_values <- subset(p_values, (BD2 %in% c('SR','Redundancy') & random == 'Primer')|(BD2 %in% c('MNTD','GD') & random == 'Barcode') )
p_values$BD2 <- factor(p_values$BD2,
                       levels = c('SR','Redundancy','MNTD','GD'),
                       labels = c('Species Richness', 'Functional Redundancy', 'MNTD', 'Genetic Diversity'))
p_values$type <- factor(p_values$type,
                        levels = c('reclassified', 'traditional'),
                        labels = c('Reanalyzed','Reported/Conventional'))

y_max_df <- res %>%
  group_by(BD2) %>%
  summarise(y_max = max(y, na.rm = TRUE)) %>%
  mutate(y_reclassified = y_max * 0.75,
         y_traditional  = y_max * 0.65)
p_values2 <- left_join(p_values, y_max_df, by = "BD2")
p_values2 <- p_values2 %>%
  mutate(y_pos = ifelse(type == "Reanalyzed", y_reclassified, y_traditional))

p <- ggplot(res, aes(x = log10(x), y = y, color = BD2, fill = BD2, linetype = type, shape = type)) +
  facet_wrap(~BD2, scales = 'free_y', ncol = 4)+
  # geom_point(color = 'grey90', fill = 'transparent')+
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_ribbon(aes(ymin = LCL_95, ymax = UCL_95), alpha = 0.2, data = dt_avg) +
  geom_smooth(method = 'glm', formula = y ~ x, level = 0.95, lwd = 1, alpha = 0.1)+
  geom_text(data = p_values2, aes(x = -Inf, y = y_pos, label = p_label), hjust = -0.1, vjust = 0,inherit.aes = FALSE, size = 4)+
  scale_color_manual(values = useCors)+
  scale_linetype_manual(values = c('solid', 'dashed')) +
  scale_shape_manual(values = 1:2)+
  scale_x_continuous(breaks = c(0,2,4,6), labels = c(1,100,10^4,10^6))+
  xlab('Catchment Area') +
  ylab('Slope (Biodiversity/Human activity)') +
  theme_bw() +
  theme(
    legend.position = 'none',
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = 'grey')
  )
p
pdf('../image/Figure_05_SHR.pdf', height = 4, width = 16)
p
dev.off()
############################################################
####Figure_S_XX: plotting ZSAR based on downsampling & reshuffles:human
############################################################
res <- readRDS('../result/Res_02_SAR_based_Full_all.RDS')
res <- subset(res$SHR_area, BD %in% all_BioDiv_para)
res <- subset(res, (BD2 %in% c('SR','Redundancy') & random == 'Primer')|(BD2 %in% c('MNTD','GD') & random == 'Barcode') )
res$BD2 <- factor(res$BD2,
                  levels = c('SR','Redundancy','MNTD','GD'),
                  labels = c('Species Richness', 'Functional Redundancy', 'MNTD', 'Genetic Diversity'))
res$type <- factor(res$type,
                   levels = c('reclassified', 'traditional'),
                   labels = c('Reanalyzed','Reported/Conventional'))

bio_factors <- unique(res$BD2)
results <- NULL
for (i in 1:length(bio_factors)) {
  cat('Processing ', as.character(bio_factors[i]), '...\n')
  col = useCors[col_num[i]]
  
  ##calculating slopes
  dt <- subset(res, BD2 == bio_factors[i])
  slopes <- dt %>%
    group_by(DT, random, id, model, x_name, y_name, BD, BD2, type)%>%
    reframe({model_fit <- lm(y ~ log10(x), data = cur_data())
    tidy_fit <- tidy(model_fit)
    slope <- tidy_fit %>% filter(term == "log10(x)") %>% pull(estimate)
    tibble(slope = slope) })
  add <- slopes %>%
    group_by(random, model, x_name, BD, BD2, type) %>%
    mutate(down_median = median(slope[DT == "downsampled"], na.rm = TRUE),
           slope = ifelse(down_median < 0, -slope, slope),
           down_median = abs(down_median) ) %>%
    summarise(down_median = first(down_median),
              reshuffled_prop = mean(slope[DT == "reshuffled"] > down_median, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(p_label = ifelse(reshuffled_prop < 0.001,"p < 0.001",paste0("p = ", sprintf("%.3f", reshuffled_prop)) ))
  annotation_df <- add %>% mutate(x = 1.5, y = Inf,  vjust = 1.5 )
  
  p <- ggplot(slopes, aes(x = DT, y = slope)) +
    facet_wrap(~ type) +
    geom_violin(trim = F, fill = 'transparent', color = col) +  
    geom_boxplot(width =0.45, outliers = T, fill = 'transparent', color = col) +
    geom_point(alpha = 0.5, shape = 21, color = col)+
    geom_text(data = annotation_df, aes(x = x, y = y, label = p_label), inherit.aes = FALSE, vjust = annotation_df$vjust, size = 5) +
    xlab('') + ylab('Area effects on the diversity-human slope') +
    theme_bw() +
    theme(
      legend.position = 'none',
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = 'grey') )
  p
  ##archive
  results$add <- p; names(results)[i] <- bio_factors[i]
}
pdf('../image/Figure_SI_XX_BioDiv_trend_SHR_AreaEffects.pdf', width = 10, height = 10)
ggarrange(results[[1]], results[[2]], results[[3]], results[[4]], ncol = 2, nrow = 2, align = 'hv', labels = letters[1:4])
dev.off()
