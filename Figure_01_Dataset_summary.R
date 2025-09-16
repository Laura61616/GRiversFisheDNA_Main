# READNE: dataset summary, 10.04.2025
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
####A01: Map of datasets
#############################################################
##summary river IDs
info <- read.csv('SiteInfo_not_Combine.csv')
nSites <- info %>% group_by(N_RivID) %>% summarise(nSite = length(unique(N_SiteID))); colnames(nSites)[1] <- 'RivID'

info <- read.csv('DataInfo_not_Combine.csv')
Riv_coord <- unique(info[,c('RivID','X','Y')])

data <- readRDS('Use_All_Sites_modelData_BioDiv.RDS')
BD <- data$BD
tmp <- BD %>% group_by(RivID) %>% summarise(nRe = sum(!is.na(SR_Re)), nOri = sum(!is.na(SR_Ori)))
tmp$nRe <- case_when(tmp$nRe > 0 ~ TRUE, TRUE ~ FALSE); tmp$nOri <- case_when(tmp$nOri > 0 ~ TRUE, TRUE ~ FALSE)

riv_summary <- Reduce(function(x,y) merge.data.frame(x,y,by = 'RivID', all = T), list(Riv_coord, nSites, tmp))
riv_summary <- riv_summary%>%
  mutate(type = case_when(nRe&nOri ~ 'Both', nRe&!nOri ~ 'Reanalyzed only', nOri&!nRe ~ 'Reported only'))

riv_summary%>%group_by(type)%>%summarise(nRiv = n(), nSites = sum(nSite))
write.table(riv_summary, '../result/Res_01_Summary_Rivers.csv', row.names = F, sep = ',', quote = F)
##plotting map
rivers <- read.csv('../result/Res_01_Summary_Rivers.csv')
sf_data <- rivers %>% st_as_sf(coords = c("X", "Y"), crs = 4326)
sf_data$type <- factor(sf_data$type, levels = c('Both',"Reanalyzed only","Reported only")
                       ,labels = c("Original and reanalyzed (51 catchments, 1212 sites)",
                                   "Reanalyzed only (14 catchments, 125 sites)",
                                   "Original only (48 catchments, 481 sites)"))
sites <- read.csv('fish_div_HG_V05_use.csv')
sf_sites <- sites %>% st_as_sf(coords = c("lon_HGV3", "lat_HGV3"), crs = 4326)
AllBasins <- st_read('Basin042017_3119.shp')
UseBasins <- st_read('00_S04_UseBasins_V03.shp')

{
  mapworld <- borders("world",
                      fill = 'grey90',
                      colour = 'grey90')
  mp <- ggplot() + 
    mapworld +
    labs(x ='',y="")+
    theme_bw()+
    theme(legend.position = c(0.2,0.5),
          legend.background = element_blank(),
          legend.title = element_blank(),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.background = element_blank())
}
p <- mp+
  geom_sf(aes(geometry = `geometry`), data = AllBasins, color = 'white',fill = 'transparent') +
  geom_sf(aes(geometry = `geometry`, fill = type ), data = sf_data, alpha = 0.8, color = 'grey50', shape = 21, size = 3) +
  scale_fill_manual(values = c( "#46732EFF", "#FED439FF", "#709AE1FF" )) +
  geom_sf(aes(geometry = `geometry` ),  data = sf_sites, alpha = 0.8, shape = 3, size = 1, stroke = 0.2, color = 'grey40', show.legend = F)
p

pdf('../image/Figure_01_AllSites_map3.pdf', width = 15, height = 5)
print(p)
dev.off()

#############################################################
####A02: distribution of BioDiv Parameters: alpha
#############################################################
BioDiv <- read.csv('Use_All_Sites_Multifacet_Div_all.csv')
bio_factors <- intersect(colnames(BioDiv), all_BioDiv_para)
pp <- NULL
for (i in 1:length(bio_factors)) {
  x <- na.omit(BioDiv[,bio_factors[i]])
  p1 <- ggplot(data.frame(x = x), aes(x = x)) +
    geom_density(fill = useCors[col_num[i]], color =useCors[col_num[i]], alpha = 0.5) +
    labs(x = fac_names[i], y = "Density") +
    theme_bw()+
    theme(legend.position = 'none',
          legend.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.background = element_blank())
  p1
  
  pp$add <- p1; names(pp)[i] <- bio_factors[i]
}

g1 <- ggarrange(pp$SR_Re, pp$Redundancy_Re, pp$MNTD_Re, pp$GD_Re,
                pp$SR_Ori, pp$Redundancy_Ori, pp$MNTD_Tre, pp$GDvalue.COI_Tre,
                ncol = 4, nrow = 2, align = 'hv')
g1

pdf('../image/Figure_01_BioDiv_density.pdf', width = 20, height = 8)
g1
dev.off()
