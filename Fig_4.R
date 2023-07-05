# --- --- --- --- --- --- --- --- --- --- --- ---
#
# Animal-borne sensors as a biologically informed lens on a changing climate
#
# Diego Ellis Soto, diego.ellissoto@yale.edu
#
# Figure 4: 
#
#
# --- --- --- --- --- --- --- --- --- --- --- ---

require(ggplot2);require(tidyverse);require(lubridate);require(plyr);require(sf);require(raster);require(sp);
require(rgdal);require(gridExtra)
library(plyr);library(ggplot2);library(lubridate);require(lubridate);library(dplyr);library(data.table);library(gridExtra);require(raster); require(tidyverse);require(sf);require(sp);require(rgdal);require(raster);require(plyr);require(lubridate);require(dplyr);require(ggplot2);require(tidyverse);require(viridis);require(rasterVis);library(geojsonio);library(sp);library(rgeos);library("ggmap");require(mapedit);require(sf);require(sf);require(maptools);require(shape);require(RColorBrewer);require(GISTools);require(maps);library(mgcv);library(MASS);require(mapview)
cm.cols1=function(x,bias=1) { colorRampPalette(c('grey90','steelblue4','steelblue1','gold','red1','red4'),bias=bias)(x)}
library(rasterVis)

# load('/Users/diegoellis/Desktop/Aqua_tmp.rdata')
load('/Users/diegoellis/Desktop/Submit_NCC_ABI/R4/Data/Tortoise/Aqua_tmp.rdata')
tor_nov_2014_LST <- raster::raster('/Users/diegoellis/Desktop/Submit_NCC_ABI/R4/Data/Tortoise/diego_MODIS_meanLST2014_november_tortuga.tif')
tor_nov_2014_LST = (tor_nov_2014_LST  * 0.02)   - 273.15 # Convert kelvin to celsius
tor_aug_2016_LST <- raster::raster('/Users/diegoellis/Desktop/Submit_NCC_ABI/R4/Data/Tortoise/diego_MODIS_meanLST2016_Aug_for_tortoise.tif')
tor_aug_2016_LST = (tor_aug_2016_LST  * 0.02)   - 273.15 # Convert kelvin to celsius

tor_nov_2016 <- Aqua_tmp %>% filter(Year == 2016, month == 8)
nrow(tor_nov_2016)
range(tor_nov_2016$SRTM)
tor_nov_2014 <- Aqua_tmp %>% filter(Year == 2014, month == 11)
nrow(tor_nov_2014)
range(tor_nov_2014$SRTM)

# Tortoise
tor_nov_2016_sp <- sp::SpatialPointsDataFrame(coords = tor_nov_2016[,c('Long', 'Lat')],
                                          data = tor_nov_2016,
                                          proj4string = sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))


tor_nov_2014_sp <- sp::SpatialPointsDataFrame(coords = tor_nov_2014[,c('Long', 'Lat')],
                                          data = tor_nov_2014,
                                          proj4string = sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# Weather station:

ibut_locs <- read.csv('/Users/diegoellis/Desktop/Submit_NCC_ABI/R4/Data/Tortoise/Ibutton_location.csv')
ibut_locs = ibut_locs[ibut_locs$long <  -90.34264,] # Remote Cerro Fatal
ibutton_sp <-  sp::SpatialPointsDataFrame(sp::coordinates(ibut_locs[, c('long', 'lat')]),data=ibut_locs, proj4string = sp::CRS("+proj=longlat +ellps=WGS84") )

tor_box_2016 <- raster::bbox(raster::extent(tor_nov_2016_sp) * 1.2)
# Look whever the bbox of tortoises is smaller or bigger than the one of weather data:
if(tor_box_2016[1] < raster::bbox(ibutton_sp)[1]){print('Tortoise xmin is smaller, keep it')}else{print('Ibutton xmax is bigger, increase the bbox'); tor_box_2016[1] <- raster::bbox(ibutton_sp)[1] }
if(tor_box_2016[2] < raster::bbox(ibutton_sp)[2]){print('Tortoise xmin is smaller, keep it')}else{print('Ibutton xmax is bigger, increase the bbox'); tor_box_2016[2] <- raster::bbox(ibutton_sp)[2] }
if(tor_box_2016[3] > raster::bbox(ibutton_sp)[3]){print('Tortoise ymin is smaller, keep it')}else{print('Ibutton ymin is bigger, increase the bbox'); tor_box_2016[3] <- raster::bbox(ibutton_sp)[3] }
if(tor_box_2016[4] > raster::bbox(ibutton_sp)[4]){print('Tortoise ymax is smaller, keep it')}else{print('Ibutton ymax is bigger, increase the bbox'); tor_box_2016[4] <- raster::bbox(ibutton_sp)[4] }

tor_box_2014 <- raster::bbox(raster::extent(tor_nov_2014_sp) * 1.2)
# Look whever the bbox of tortoises is smaller or bigger than the one of weather data:
if(tor_box_2014[1] < raster::bbox(ibutton_sp)[1]){print('Tortoise xmin is smaller, keep it')}else{print('Ibutton xmax is bigger, increase the bbox'); tor_box_2014[1] <- raster::bbox(ibutton_sp)[1] }
if(tor_box_2014[2] < raster::bbox(ibutton_sp)[2]){print('Tortoise xmin is smaller, keep it')}else{print('Ibutton xmax is bigger, increase the bbox'); tor_box_2014[2] <- raster::bbox(ibutton_sp)[2] }
if(tor_box_2014[3] > raster::bbox(ibutton_sp)[3]){print('Tortoise ymin is smaller, keep it')}else{print('Ibutton ymin is bigger, increase the bbox'); tor_box_2014[3] <- raster::bbox(ibutton_sp)[3] }
if(tor_box_2014[4] > raster::bbox(ibutton_sp)[4]){print('Tortoise ymax is smaller, keep it')}else{print('Ibutton ymax is bigger, increase the bbox'); tor_box_2014[4] <- raster::bbox(ibutton_sp)[4] }

tor_aug_2016_LST_crop <- raster::crop(tor_aug_2016_LST, tor_box_2016)
tor_nov_2014_LST_crop <- raster::crop(tor_nov_2014_LST, tor_box_2014)

SRTM = raster::raster('/Users/diegoellis/Desktop/Submit_NCC_ABI/R4/Data/Tortoise/SRTM_Santa_Cruz.tif')
# Hasta aca llegue:
plot(SRTM, col = alpha('gray49', 0.8))

pdf(paste0(outdir, '/tortoise_inset_space_2014.pdf'))
par(mfrow=c(1,3))
# MODIS
raster::plot(SRTM, col = alpha('gray49', 0.8), legend=FALSE, ylim = c(-0.75, -0.655), xlim = c(-90.5,-90.388),xlab="",ylab="", axes = F, box = FALSE)
raster::plot(tor_nov_2014_LST_crop, zlim=range(raster::values(tor_nov_2014_LST_crop$LST_Night_1km), na.rm=T),col=cm.cols1(100),las=1,bty='n', add=T)
raster::plot(raster::extent(tor_nov_2014_LST_crop), col = 'black', add=T, lwd = 1.5)
raster::scalebar(2, xy=c(-90.49741, -0.7413617), type='bar', divs=2)
north.arrow(xb=-90.4953, yb=-0.7248967, len=0.0018, lab="N")  
# Tortoise
raster::plot(SRTM, col = alpha('gray49', 0.8), legend=FALSE, ylim = c(-0.75, -0.655), xlim = c(-90.5,-90.388),xlab="",ylab="", axes = F, box = FALSE)

tor_nov_2014_df = data.frame(tor_nov_2014_sp)
tor_1 = tor_nov_2014_df %>% filter(individual_id == unique(tor_nov_2014_df$individual_id)[1]) %>% dplyr::select(Long, Lat)
tor_2 = tor_nov_2014_df %>% filter(individual_id == unique(tor_nov_2014_df$individual_id)[2]) %>% dplyr::select(Long, Lat)
tor_3 = tor_nov_2014_df %>% filter(individual_id == unique(tor_nov_2014_df$individual_id)[3]) %>% dplyr::select(Long, Lat)
lines(tor_1, add=T, lwd = 2)
lines(tor_2, add=T, lwd = 2)
lines(tor_3, add=T, lwd = 2)

# plot(ele_1, add=T, lwd = 2)
# plot(ele_2, add=T, lwd = 2)
# plot(ele_3, add=T, lwd = 2)
raster::plot(raster::extent(tor_nov_2014_LST_crop), col = 'black', add=T, lwd = 1.5)
# IButton:
raster::plot(SRTM, col = alpha('gray49', 0.8), legend=FALSE, ylim = c(-0.75, -0.655), xlim = c(-90.5,-90.388),xlab="",ylab="", axes = F, box = FALSE)
# Remove one Ibutton from there: the one below -0.72 for plotting reasons?
ibutton_sp_tmp <- ibutton_sp[ibutton_sp$lat > -0.72,]
sp::plot(ibutton_sp_tmp, add=T, pch = 1, cex = 1 , col = alpha('black', 0.8))
raster::plot(raster::extent(tor_nov_2014_LST_crop), col = 'black', add=T, lwd = 1.5)

# Now plot the annual temperature profiles:



ibutton_sp$elev = raster::extract(SRTM, ibutton_sp)
tmp_df_2016_temp_sp <- ibutton_sp
# Round to closest 50
tmp_df_2016_temp_sp$elev <- round_any(tmp_df_2016_temp_sp$elev, 50)
tmp_df_2016_temp_sp_df = as.data.frame(tmp_df_2016_temp_sp)
tmp_df_2016_temp_sp_df = tmp_df_2016_temp_sp_df[,c('long', 'lat', 'elev')]
tmp_df_2016_temp_sp_df$elev <- round_any(tmp_df_2016_temp_sp_df$elev, 50)
names(tmp_df_2016_temp_sp_df)[3] <- 'Altitude'


temp_ibutton_long_term <- read.csv('/Users/diegoellis/Desktop/Submit_NCC_ABI/R4/Data/Tortoise/temp_ibutton_long_term.csv') %>%
  filter(year >= range(Aqua_tmp$Year)[1], year <= range(Aqua_tmp$Year)[2]) %>%
  mutate(Date = mdy_hm(dattimhour),
         hour = hour(Date),
         jday = yday(Date),
         month = month(Date)) %>%
  subset(hour >= 0 & hour < 4)

names(temp_ibutton_long_term)
names(temp_ibutton_long_term)[2] <- 'Altitude'
temp_ibutton_long_term <- temp_ibutton_long_term[! temp_ibutton_long_term$site == 'CF',]
ibutton_temperature = left_join(temp_ibutton_long_term, tmp_df_2016_temp_sp_df)
# ibutton_temperature <- ibutton_temperature[ibut_locs$lat > -0.72,]
ibutton_temperature[which(ibutton_temperature$temp == max(ibutton_temperature$temp)),] # Where are the hottest temperature?! 

# Now inset in time:
tmp_aq_n_tor_aq_T = ggplot(Aqua_tmp, aes(x = Jday, y = temperature_calibrated))  + geom_point(alpha = 1/20) + geom_smooth(size = 1.2) + ylab('Tortoise \n collected temperature') + theme_classic() + ylim(0, 30) + theme(text = element_text(size=20))
tmp_aq_n_tor_T_MODIS = ggplot(Aqua_tmp, aes(x = Jday, y = LST)) + geom_point(alpha = 1/20) + geom_smooth(size = 1.2)  +
  theme_classic() + ylab('Aqua Night \n Land Surface temperature') + ylim(0, 30) + theme(text = element_text(size=20))
tmp_aq_n_ibutton_aq_T = ggplot(ibutton_temperature, aes(x = jday, y = temp)) + geom_point(alpha = 1/20) + geom_smooth(size=1.2) +
  theme_classic() + ylab('Ibutton \n Temperature') + ylim(0, 30) + theme(text = element_text(size=20))
pdf(paste0(outdir, '/tortoise_inset_time.pdf'), height = 12, width = 8)
grid.arrange(tmp_aq_n_tor_T_MODIS, tmp_aq_n_ibutton_aq_T, tmp_aq_n_tor_aq_T, ncol = 1)
dev.off()


# Model:

source('/Users/diegoellis/projects/Phd_Chapter_1/ABI_Phd_Chapter_1_code/0_source_packages_ABI.R')
require(move)
require(plyr)
require(dplyr)
require(data.table)
require(lubridate)
require(rpart)
library(randomForest)
require(biomod2)
library(adehabitatHS)
library(pROC)
library(ggplot2)
library(earth)
library(MASS)
library(nnet)
library(gbm)
library(cowplot)
require(gam)
require(Metrics)
require(gridExtra)
cm.cols1=function(x,bias=1) { colorRampPalette(c('grey90','steelblue4','steelblue1','gold','red1','red4'),bias=bias)(x)}

# load('/Users/diegoellis/projects/Phd_Chapter_1/ABI_Phd_Chapter_1_code/Rdata_objects/Biologging_aqua_night.rdata')
# # load('/Users/diegoellis/Desktop/Submit_NCC_ABI/R4/Data/Tortoise/Aqua_tmp.rdata')
# load('/Users/diegoellis/projects/development/Tortoise_data/Biologging_aqua_night.rdata')
# load('/Users/diegoellis/projects/Phd_Chapter_1/Data/Biologging_aqua_night_original_2019.rdata')

LST <- raster::raster('/Users/diegoellis/Desktop/Submit_NCC_ABI/R4/Data/Tortoise/MYD11A2.A2016017_res_compiled_msk_translated_WGS84.tif')
SRTM <- raster::raster('/Users/diegoellis/Desktop/Submit_NCC_ABI/R4/Data/Tortoise/SRTM_30m_Gal_res.tif')
Aspect_sin <- raster::raster('/Users/diegoellis/Desktop/Submit_NCC_ABI/R4/Data/Tortoise/SRTM_30m_Gal_res_aspect_SINUS.tif')
Aspect_cosin <- raster::raster('/Users/diegoellis/Desktop/Submit_NCC_ABI/R4/Data/Tortoise/SRTM_30m_Gal_res_aspect_COSINUS.tif')
roughness <- raster::raster('/Users/diegoellis/Desktop/Submit_NCC_ABI/R4/Data/Tortoise/SRTM_30m_Gal_res_roughness.tif')
TPI <- raster::raster('/Users/diegoellis/Desktop/Submit_NCC_ABI/R4/Data/Tortoise/SRTM_30m_Gal_res_topographic_index.tif')
TRI <- raster::raster('/Users/diegoellis/Desktop/Submit_NCC_ABI/R4/Data/Tortoise/SRTM_30m_Gal_res_terrain_ruggedness_index.tif')
Slope = raster::raster('/Users/diegoellis/Desktop/Submit_NCC_ABI/R4/Data/Tortoise/SRTM_30m_Gal_res_percent_slope.tif')
# LST = raster('/Users/diegoellis/Downloads/drive-download-20210314T202458Z-001/MYD11A2.A2015361_res_compiled_msk_translated_WGS84.tif')
LST = raster::raster('/Users/diegoellis/Desktop/Submit_NCC_ABI/R4/Data/Tortoise/MYD11A2.A2015361_res_compiled_msk_translated_WGS84.tif')
# SRTM2 = raster('/Users/diegoellis/Dropbox/Marius_Galapagos/Inputs/RS_Data_Santa_Cruz/SRTM/SRTM_Santa_Cruz.tif')
SRTM2 = raster::raster('/Users/diegoellis/Desktop/Submit_NCC_ABI/R4/Data/Tortoise/SRTM_Santa_Cruz.tif')
LST = raster::crop(LST, SRTM2)
LST_coarse = LST
LST <- raster::disaggregate(LST, fact=c(30, 30))
names(LST) = 'LST'
names(SRTM) = 'SRTM'
names(Aspect_cosin) = 'Aspect_cosin'
names(Aspect_sin) = 'Aspect_sin'
names(roughness) = 'roughness'
names(TPI) = 'TPI'
names(TRI) = 'TRI'
names(Slope) <- 'Slope'
SRTM = raster::crop(SRTM,LST)
# Crop to Santa Cruz:
Aspect_cosin = raster::crop(Aspect_cosin,LST)
Aspect_sin = raster::crop(Aspect_sin,LST)
roughness = raster::crop(roughness,LST)
TPI = raster::crop(TPI,LST)
TRI = raster::crop(TRI,LST)
Slope = raster::crop(Slope, LST)
LST = (LST  * 0.02)   - 273.15 # Convert kelvin to celsius
senv = raster::stack(LST, SRTM, Aspect_cosin,Aspect_sin,Slope,TPI, TRI,  roughness)
# remove TPI is crap
senv  = raster::dropLayer(senv, 6) 

Aqua_night = Aqua_tmp # Cambie esto para solo cargar las tortugas de arriba ####
names(Aqua_night)
Aqua_night_sub =  Aqua_night[,c('temperature_calibrated', 'LST', 'SRTM',
                                'Ascpect','roughness',  'Percent_slope_degree', 'Long', 'Lat') ]
# Aqua_night_sub = Aqua_night[,c(1:2,5:6 ,8:9)]
# Drop Ascpect, roughness, percent slope degree and SRTM column and just extract everything again
Aqua_night_sub <- Aqua_night_sub[,c(1:2,7:8)]
pdf_tort <- sp::SpatialPointsDataFrame(coords = Aqua_night_sub[,c('Long','Lat')], data = Aqua_night_sub, proj4string = sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
pdf_tort$SRTM <- raster::extract(SRTM , pdf_tort)
pdf_tort$Aspect_cosin <- raster::extract(Aspect_cosin , pdf_tort)
pdf_tort$Aspect_sin <- raster::extract(Aspect_sin , pdf_tort)
pdf_tort$Slope <- raster::extract(Slope , pdf_tort)
pdf_tort$TPI <- raster::extract(TPI , pdf_tort)
pdf_tort$TRI <- raster::extract(TRI , pdf_tort)
pdf_tort$roughness <- raster::extract(roughness , pdf_tort)

Aqua_night_sub = as.data.frame(pdf_tort)
Aqua_night_sub = Aqua_night_sub %>% dplyr::select(-c(Lat.1, Long.1))# [,-11:-12] # Aqua_night_sub = Aqua_night_sub[,-14:-15]
Aqua_night_sub = Aqua_night_sub[complete.cases(Aqua_night_sub),]

# [1] MODEL GAM 1 ####
if(is.element("package:mgcv", search())) detach("package:mgcv") ## make sure the mgcv package is not loaded to avoid conflicts
library(gam)
require(car)

gam2 = gam::gam(temperature_calibrated ~ s(LST,2) + s(SRTM,2) + s(Aspect_sin,2) + s(Aspect_cosin,2) + s(roughness,2)+ s(Slope,2), data=Aqua_night_sub)

pred_model_gam2 = raster::predict(senv, gam2)

cm.cols1=function(x,bias=1) { colorRampPalette(c('grey90','steelblue4','steelblue1','gold','red1','red4'),bias=bias)(x)}

par(mfrow=c(1,2))
# MODIS:
raster::plot(LST,zlim=range(values(LST$LST),na.rm=T), col=cm.cols1(100),las=1,box=FALSE,bty='n', main = 'MODIS Night time Temperature',  ylim = c(-0.8,-0.45))
scalebar(10, type='bar', divs=2, below = 'km')
north.arrow(xb=-90.52827, yb=-0.428331, len=0.008, lab="N")  
# Model:
raster::plot(pred_model_gam2,zlim=range(values(pred_model_gam2$layer),na.rm=T), col=cm.cols1(100),las=1,box=FALSE, main = 'Tortoise microclimate\n at 30m GAM gam2',  ylim = c(-0.8,-0.45))
scalebar(10, type='bar', divs=2, below = 'km')
north.arrow(xb=-90.52827, yb=-0.428331, len=0.008, lab="N")  

# Assume RCP 4.5 so 2 degree celsius increse:
LST_future = LST + 2
GAM_future = pred_model_gam2 + 2

# Threshold for a species that likes 20-28 degrees
#Convert rasters to binary:
LST_future[LST_future$LST < 20,] <- 0
LST_future[LST_future$LST > 28,] <- 0
LST_future[LST_future$LST > 0,] <- 1

GAM_future[GAM_future$layer < 20,] <- 0
GAM_future[GAM_future$layer > 28,] <- 0
GAM_future[GAM_future$layer > 0,] <- 1

cpal <- c('grey59', 'firebrick3')
r2_LST <- ratify(LST_future)
levelplot(r2_LST,col.regions=cpal,att='ID')
r2_GAM <- ratify(GAM_future)
levelplot(r2_GAM,col.regions=cpal,att='ID')
