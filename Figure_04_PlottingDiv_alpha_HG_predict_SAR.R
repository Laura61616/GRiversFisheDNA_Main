# READNE: plotting species-area relationship, 17.06.2025
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
####Figure_04: plotting ZASR
############################################################
##Latitude trends
data <- read.csv('../result/Res_02_SAR_based_Full_perRiver.csv')
head(data)

info <- read.csv('DataInfo_not_Combine.csv')
info <- unique(info[,c(4:7)])
head(info)

data <- merge.data.frame(info, data, by = 'RivID')
data <- subset(data, BD %in% all_BioDiv_para)
data <- subset(data, (BD2 %in% c('SR','Redundancy') & random == 'Primer')|(BD2 %in% c('MNTD','GD') & random == 'Barcode'))
data$BD2 <- factor(data$BD2,
                   levels = c('SR','Redundancy','MNTD','GD'),
                   labels = c('Species Richness', 'Functional Redundancy', 'MNTD', 'Genetic Diversity'))
data$type <- factor(data$type,
                    levels = c('reclassified', 'traditional'),
                    labels = c('Reanalyzed','Reported/Conventional'))

stats_df <- data %>%
  group_by(BD, type, BD2) %>%
  summarise(model = list(gam(trend ~ abs(Y), data = cur_data()))) %>%
  rowwise() %>%
  mutate(R2 = round(summary(model)$r.sq, 2),
         P = signif(summary(model)$p.pv[2], 3),
         label = ifelse(P < 0.001,
                        paste0("R² = ", R2, ", p < 0.001"),
                        paste0("R² = ", R2, ", p = ", formatC(P, format = "f", digits = 3))) )
positions <- data %>%
  group_by(BD, BD2, type) %>%
  summarise(x = min(abs(Y), na.rm = TRUE) + 1, .groups = "drop") %>%
  group_by(BD2) %>%
  mutate(ymax_overall = max(data$trend[data$BD2 == unique(BD2)], na.rm = TRUE),
         offset = ifelse(type == "Reanalyzed", 0.05, -0.05) * ymax_overall,
         y = ymax_overall + offset)

stats_df <- left_join(stats_df, positions, by = c("BD2", "type"))

p1 <- ggplot(data, aes(x = abs(Y), y = trend, color = BD2, fill = BD2, linetype = type, shape = type)) +
  facet_wrap(~BD2, scales = 'free_y', ncol = 4)+
  geom_point()+
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey')+
  scale_color_manual(values = useCors)+
  geom_smooth(method = 'gam', formula = y ~ x, level = 0.95, lwd = 1, alpha = 0.1)+
  geom_text(data = stats_df, aes(x = x, y = y, label = label), inherit.aes = FALSE, hjust = 0, size = 3.5)+
  xlab('Latitude') +
  ylab('Slope (Biodiversity/catchment area)') +
  theme_bw() +
  theme(
    legend.position = 'right',
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = 'grey') )
p1

##climate verus ZSAR
res <- readRDS('../result/Res_02_SAR_based_Full_all.RDS')
res <- subset(res$SAR_bio01, BD %in% all_BioDiv_para & DT == "downsampled")
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
p_values <- read.csv('../result/Res_02_SAR_based_Full_TempEffects_significance.csv')
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
  mutate(y_reclassified = y_max * 1.05,
         y_traditional  = y_max * 0.95)
p_values2 <- left_join(p_values, y_max_df, by = "BD2")
p_values2 <- p_values2 %>%
  mutate(y_pos = ifelse(type == "Reanalyzed", y_reclassified, y_traditional))

p2 <- ggplot(res, aes(x = x, y = y, color = BD2, fill = BD2, linetype = type, shape = type)) +
  facet_wrap(~BD2, scales = 'free_y', ncol = 4)+
  # geom_point(color = 'grey90', fill = 'transparent')+
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_ribbon(aes(ymin = LCL_95, ymax = UCL_95), alpha = 0.2, data = dt_avg) +
  geom_smooth(method = 'glm', formula = y ~ x, level = 0.95, lwd = 1, alpha = 0.1)+
  geom_text(data = p_values2, aes(x = -Inf, y = y_pos, label = p_label), hjust = -0.1, vjust = 0,inherit.aes = FALSE, size = 4)+
  scale_color_manual(values = useCors)+
  scale_linetype_manual(values = c('solid', 'dashed')) +
  scale_shape_manual(values = 1:2)+
  xlab('Temperature') +
  ylab('Slope (Biodiversity/catchment area)') +
  theme_bw() +
  theme(
    legend.position = 'none',
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = 'grey')
  )
p2

##Human verus ZSAR
res <- readRDS('../result/Res_02_SAR_based_Full_all.RDS')
res <- subset(res$SAR_human, BD %in% all_BioDiv_para & DT == "downsampled")
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
p_values <- read.csv('../result/Res_02_SAR_based_Full_HumanEffects_significance.csv')
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

p3 <- ggplot(res, aes(x = x, y = y, color = BD2, fill = BD2, linetype = type, shape = type)) +
  facet_wrap(~BD2, scales = 'free_y', ncol = 4)+
  # geom_point(color = 'grey90', fill = 'transparent')+
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_ribbon(aes(ymin = LCL_95, ymax = UCL_95), alpha = 0.2, data = dt_avg) +
  geom_smooth(method = 'glm', formula = y ~ x, level = 0.95, lwd = 1, alpha = 0.1)+
  geom_text(data = p_values2, aes(x = -Inf, y = y_pos, label = p_label), hjust = -0.1, vjust = 0,inherit.aes = FALSE, size = 4)+
  scale_color_manual(values = useCors)+
  scale_linetype_manual(values = c('solid', 'dashed')) +
  scale_shape_manual(values = 1:2)+
  xlab('Human activity') +
  ylab('Slope (Biodiversity/catchment area)') +
  theme_bw() +
  theme(
    legend.position = 'none',
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = 'grey')
  )
p3

pdf('../image/Figure_04_SAR2.pdf', height = 12, width = 20)
ggarrange(p1,p2,p3,nrow = 3,align = 'hv',common.legend = T, legend = 'right')
dev.off()
############################################################
####Figure_04: plotting ZASR based on "real data"
############################################################
data <- read.csv('../result/Res_02_SAR_based_Full_perRiver.csv')
head(data)

info <- read.csv('DataInfo_not_Combine.csv')
info <- unique(info[,c(4:7)])
head(info)

data <- merge.data.frame(info, data, by = 'RivID')
data <- subset(data, BD %in% all_BioDiv_para)
data <- subset(data, (BD2 %in% c('SR','Redundancy') & random == 'Primer')|(BD2 %in% c('MNTD','GD') & random == 'Barcode'))
data$BD2 <- factor(data$BD2,
                  levels = c('SR','Redundancy','MNTD','GD'),
                  labels = c('Species Richness', 'Functional Redundancy', 'MNTD', 'Genetic Diversity'))
data$type <- factor(data$type,
                    levels = c('reclassified', 'traditional'),
                    labels = c('Reanalyzed','Reported/Conventional'))

##distribution along temperature
stats_df <- data %>%
  group_by(BD, type, BD2) %>%
  summarise(model = list(gam(trend ~ bio01_unscale, data = cur_data()))) %>%
  rowwise() %>%
  mutate(R2 = round(summary(model)$r.sq, 2),
         P = signif(summary(model)$p.pv[2], 3),
         label = ifelse(P < 0.001,
                        paste0("R² = ", R2, ", p < 0.001"),
                        paste0("R² = ", R2, ", p = ", formatC(P, format = "f", digits = 3))) )
positions <- data %>%
  group_by(BD, BD2, type) %>%
  summarise(x = min(bio01_unscale, na.rm = TRUE) + 1, .groups = "drop") %>%
  group_by(BD2) %>%
  mutate(ymax_overall = max(data$trend[data$BD2 == unique(BD2)], na.rm = TRUE),
         offset = ifelse(type == "Reanalyzed", 0.05, -0.05) * ymax_overall,
         y = ymax_overall + offset)

stats_df <- left_join(stats_df, positions, by = c("BD2", "type"))

p2 <- ggplot(data, aes(x = bio01_unscale, y = trend, color = BD2, fill = BD2, linetype = type, shape = type)) +
  facet_wrap(~BD2, scales = 'free_y', ncol = 4)+
  geom_point()+
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey')+
  scale_color_manual(values = useCors)+
  geom_smooth(method = 'gam', formula = y ~ x, level = 0.95, lwd = 1, alpha = 0.1)+
  geom_text(data = stats_df, aes(x = x, y = y, label = label), inherit.aes = FALSE, hjust = 0, size = 3.5)+
  xlab('Temperature') +
  ylab('Slope (Biodiversity/catchment area)') +
  theme_bw() +
  theme(
    legend.position = 'right',
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = 'grey') )
p2
##distribution along human activity
stats_df <- data %>%
  group_by(BD, type, BD2) %>%
  summarise(model = list(gam(trend ~ HFP_Lights2009_unscale, data = cur_data()))) %>%
  rowwise() %>%
  mutate(R2 = round(summary(model)$r.sq, 2),
         P = signif(summary(model)$p.pv[2], 3),
         label = ifelse(P < 0.001,
                        paste0("R² = ", R2, ", p < 0.001"),
                        paste0("R² = ", R2, ", p = ", formatC(P, format = "f", digits = 3))) )
positions <- data %>%
  group_by(BD, BD2, type) %>%
  summarise(x = min(HFP_Lights2009_unscale, na.rm = TRUE) + 1, .groups = "drop") %>%
  group_by(BD2) %>%
  mutate(ymax_overall = max(data$trend[data$BD2 == unique(BD2)], na.rm = TRUE),
         offset = ifelse(type == "Reanalyzed", 0.05, -0.05) * ymax_overall,
         y = ymax_overall + offset)

stats_df <- left_join(stats_df, positions, by = c("BD2", "type"))

p3 <- ggplot(data, aes(x = HFP_Lights2009_unscale, y = trend, color = BD2, fill = BD2, linetype = type, shape = type)) +
  facet_wrap(~BD2, scales = 'free_y', ncol = 4)+
  geom_point()+
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey')+
  scale_color_manual(values = useCors)+
  geom_smooth(method = 'gam', formula = y ~ x, level = 0.95, lwd = 1, alpha = 0.1)+
  geom_text(data = stats_df, aes(x = x, y = y, label = label), inherit.aes = FALSE, hjust = 0, size = 3.5)+
  xlab('Human activity') +
  ylab('Slope (Biodiversity/catchment area)') +
  theme_bw() +
  theme(
    legend.position = 'right',
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = 'grey') )
p3

pdf('../image/Figure_04_SAR_real.pdf', height = 8, width = 20)
ggarrange(p2,p3,nrow = 2,align = 'hv',common.legend = T, legend = 'right')
dev.off()
############################################################
####Figure_S_XX: plotting ZSAR based on downsampling & reshuffles:temp
############################################################
res <- readRDS('../result/Res_02_SAR_based_Full_all.RDS')
res <- subset(res$SAR_bio01, BD %in% all_BioDiv_para)
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
  
  ##plotting regression curves
  dt <- subset(res, BD2 == bio_factors[i] & DT == 'downsampled')
  dt_avg  <- dt %>%
    group_by(BD, BD2, x_name, x, x_sc, y_name, type) %>%
    summarise(
      y = mean(y, na.rm = TRUE),
      SE = mean(SE, na.rm = TRUE),
      LCL_95 = mean(LCL_95, na.rm = TRUE),
      UCL_95 = mean(UCL_95, na.rm = TRUE),
      LCL_99 = mean(LCL_99, na.rm = TRUE),
      UCL_99 = mean(UCL_99, na.rm = TRUE),
      .groups = "drop"  )
  p1 <- ggplot(dt, aes(x = x, y = y, linetype = type)) +
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
    geom_ribbon(aes(ymin = LCL_95, ymax = UCL_95), alpha = 0.2, data = dt_avg, color = col, fill = col) +
    geom_smooth(method = 'glm', formula = y ~ x, level = 0.95, lwd = 1, alpha = 0.1, color = col, fill = col)+
    scale_linetype_manual(values = c('solid', 'dashed')) +
    xlab('Temperature') +
    ylab('Slope (Biodiversity/catchment area)') +
    theme_bw() +
    theme(
      legend.position = 'none',
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = 'grey')
    )
  p1
  
  ##calculating slopes
  dt <- subset(res, BD2 == bio_factors[i])
  slopes <- dt %>%
    group_by(DT, random, id, model, x_name, y_name, BD, BD2, type)%>%
    reframe({model_fit <- lm(y ~ x, data = cur_data())
    tidy_fit <- tidy(model_fit)
    slope <- tidy_fit %>% filter(term == "x") %>% pull(estimate)
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
  
  p2 <- ggplot(slopes, aes(x = DT, y = slope)) +
    facet_wrap(~ type) +
    geom_violin(trim = F, fill = 'transparent', color = col) +  
    geom_boxplot(width =0.45, outliers = T, fill = 'transparent', color = col) +
    geom_point(alpha = 0.5, shape = 21, color = col)+
    geom_text(data = annotation_df, aes(x = x, y = y, label = p_label), inherit.aes = FALSE, vjust = annotation_df$vjust, size = 5) +
    xlab('') + ylab('Effects of temperature on ZSAR') +
    theme_bw() +
    theme(
      legend.position = 'none',
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = 'grey') )
  p2
  ##archive
  results$add <- list(p1, p2, p1 + p2 + plot_layout(ncol = 2, widths = c(3/5,2/5)) + plot_annotation(title = bio_factors[i])); 
  names(results)[i] <- bio_factors[i]
}
pdf('../image/Figure_SI_XX_BioDiv_trend_SAR_TempEffects.pdf', width = 20, height = 12)
ggarrange(results[[1]][[3]], results[[2]][[3]], results[[3]][[3]], results[[4]][[3]],
          ncol = 2, nrow = 2, align = 'hv')
dev.off()
############################################################
####Figure_S_XX: plotting ZSAR based on downsampling & reshuffles:human
############################################################
res <- readRDS('../result/Res_02_SAR_based_Full_all.RDS')
res <- subset(res$SAR_human, BD %in% all_BioDiv_para)
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
  
  ##plotting regression curves
  dt <- subset(res, BD2 == bio_factors[i] & DT == 'downsampled')
  dt_avg  <- dt %>%
    group_by(BD, BD2, x_name, x, x_sc, y_name, type) %>%
    summarise(
      y = mean(y, na.rm = TRUE),
      SE = mean(SE, na.rm = TRUE),
      LCL_95 = mean(LCL_95, na.rm = TRUE),
      UCL_95 = mean(UCL_95, na.rm = TRUE),
      LCL_99 = mean(LCL_99, na.rm = TRUE),
      UCL_99 = mean(UCL_99, na.rm = TRUE),
      .groups = "drop"  )
  p1 <- ggplot(dt, aes(x = x, y = y, linetype = type)) +
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
    geom_ribbon(aes(ymin = LCL_95, ymax = UCL_95), alpha = 0.2, data = dt_avg, color = col, fill = col) +
    geom_smooth(method = 'glm', formula = y ~ x, level = 0.95, lwd = 1, alpha = 0.1, color = col, fill = col)+
    scale_linetype_manual(values = c('solid', 'dashed')) +
    xlab('Human activity') +
    ylab('Slope (Biodiversity/catchment area)') +
    theme_bw() +
    theme(
      legend.position = 'none',
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = 'grey')
    )
  p1
  
  ##calculating slopes
  dt <- subset(res, BD2 == bio_factors[i])
  slopes <- dt %>%
    group_by(DT, random, id, model, x_name, y_name, BD, BD2, type)%>%
    reframe({model_fit <- lm(y ~ x, data = cur_data())
    tidy_fit <- tidy(model_fit)
    slope <- tidy_fit %>% filter(term == "x") %>% pull(estimate)
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
  
  p2 <- ggplot(slopes, aes(x = DT, y = slope)) +
    facet_wrap(~ type) +
    geom_violin(trim = F, fill = 'transparent', color = col) +  
    geom_boxplot(width =0.45, outliers = T, fill = 'transparent', color = col) +
    geom_point(alpha = 0.5, shape = 21, color = col)+
    geom_text(data = annotation_df, aes(x = x, y = y, label = p_label), inherit.aes = FALSE, vjust = annotation_df$vjust, size = 5) +
    xlab('') + ylab('Effects of human activity on ZSAR') +
    theme_bw() +
    theme(
      legend.position = 'none',
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = 'grey') )
  p2
  ##archive
  results$add <- list(p1, p2, p1 + p2 + plot_layout(ncol = 2, widths = c(3/5,2/5)) + plot_annotation(title = bio_factors[i])); 
  names(results)[i] <- bio_factors[i]
}
pdf('../image/Figure_SI_XX_BioDiv_trend_SAR_HumanEffects.pdf', width = 20, height = 12)
ggarrange(results[[1]][[3]], results[[2]][[3]], results[[3]][[3]], results[[4]][[3]],
          ncol = 2, nrow = 2, align = 'hv')
dev.off()
