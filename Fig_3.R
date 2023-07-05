#
# Animal-borne sensors as a biologically informed lens on a changing climate
#
# Diego Ellis Soto, diego.ellissoto@yale.edu
#
# Figure 3: 
#
#
# --- --- --- --- --- --- --- --- --- --- --- ---

# Make elephant inset Figure 3:

cm.cols1=function(x,bias=1) { colorRampPalette(c('grey90','steelblue4','steelblue1','gold','red1','red4'),bias=bias)(x)}
require(GISTools)
require(sp)
require(raster)
require(sf)
require(tidyverse)
outdir = '/Users/diegoellis/Desktop/'
# Read in elephant data and transform to a spatial object

load('/Users/diegoellis/projects/Phd_Chapter_1/Data/Aqua_night_elephant_tortoise_for_plot/Elephant/ele_flux_aq_night.Rdata')
shp_kruger <- readOGR('/Users/diegoellis/projects/Phd_Chapter_1/Data/Shapefile_kruger/WDPA_Mar2020_protected_area_873-shapefile/WDPA_Mar2020_protected_area_873-shapefile-polygons.shp')

# Load elephant data for the speicfic month: April 2009:
Aqua_night_flyovertime <- range(read.csv('/Users/diegoellis/projects/Phd_Chapter_1/Data/ThermochronTracking Elephants Kruger 2007-6340748537733918440/ThermochronTracking Elephants Kruger 2007-6340748537733918440.csv')$MODIS.Land.Surface.Temp...Emissivity.1km.Daily.Aqua.View.Time.Night, na.rm = T)
# Thin elephants by flyover time ####
ele_april_2009 <- read.csv('/Users/diegoellis/projects/Phd_Chapter_1/Data/ThermochronTracking Elephants Kruger 2007-6340748537733918440/ThermochronTracking Elephants Kruger 2007-6340748537733918440.csv') %>%
  mutate(jday = yday(timestamp),month = month(timestamp),
         year = year(timestamp)) %>% dplyr::filter(year == 2009, month == 4) %>% 
  dplyr::rename(Long = location.long,
                Lat = location.lat,
                Temp = external.temperature,
                Animal_id = individual.local.identifier,
                MODIS_Aqua_Daily_1km_night_LST = MODIS.Land.Surface.Temp...Emissivity.1km.Daily.Aqua.Land.Surface.Temperature.Night,
                MODIS_Terra_Daily_1km_night_LST = MODIS.Land.Surface.Temp...Emissivity.1km.Daily.Terra.Land.Surface.Temperature.Night,
                MODIS_Aqua_Daily_1km_day_LST = MODIS.Land.Surface.Temp...Emissivity.1km.Daily.Aqua.Land.Surface.Temperature.Day,
                MODIS_Terra_Daily_1km_day_LST = MODIS.Land.Surface.Temp...Emissivity.1km.Daily.Terra.Land.Surface.Temperature.Day)

ele_april_2009_spdf <- SpatialPointsDataFrame(coords= ele_april_2009[,c('Long', 'Lat')],
                                              data = ele_april_2009,
                                              proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
ele_april_2009$hour <- hour(ele_april_2009$timestamp)
ele_april_2009_aq_night <- ele_april_2009[ele_april_2009$hour >= Aqua_night_flyovertime[1] & ele_april_2009$hour <= Aqua_night_flyovertime[2],] %>% drop_na(Temp, MODIS_Aqua_Daily_1km_night_LST)

ele_april_2009_aq_night_spdf <- SpatialPointsDataFrame(
  coords= ele_april_2009_aq_night[,c('Long', 'Lat')],
  data = ele_april_2009_aq_night,
  proj4string = 
    CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

# --- --- --- --- --- --- --- --- --- ---
# Make a buffer around the flux tower:
# --- --- --- --- --- --- --- --- --- ---

myprj <- paste0("+proj=laea +lat_0=", round(mean(ele_april_2009_spdf$Lat)),' +lon_0=', round(mean(ele_april_2009_spdf$Long)),' +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0')
ele_april_2009_spdf_laea <- spTransform(ele_april_2009_spdf, myprj)
flux_all_sub_aq_night_laea <- spTransform(flux_all_sub_aq_night_spdf, myprj)
flux_all_sub_aq_night_laea_sf <- st_as_sf(flux_all_sub_aq_night_laea)
df_sf_buff <- st_buffer(flux_all_sub_aq_night_laea_sf, 20000)
df_sf_buff_1 <- df_sf_buff[1,] # Buffered flux tower:
buf_sp <- as(df_sf_buff, 'Spatial')
buf_sp <- spTransform(buf_sp, proj4string(shp_kruger))

world.shp=readOGR('/Users/diegoellis/Desktop/Submit_NCC_ABI/R4/Data/TM_WORLD_BORDERS_SIMPL-0.3/TM_WORLD_BORDERS_SIMPL-0.3.shp',verbose=F)
ZAF <- world.shp[world.shp$ISO3 == 'ZAF',]

# img_kruger_mean_april <- raster('/Users/diegoellis/Downloads/diego_MODIS_meanLST2008_may_for_elephant_2.tif')
img_kruger_mean_april <- raster::raster('/Users/diegoellis/Desktop/Submit_NCC_ABI/R4/Data/Elephant/diego_MODIS_meanLST2008_may_for_elephant_2.tif')
img_kruger_mean_april <- raster::crop(img_kruger_mean_april, shp_kruger)
img_kruger_mean_april <- raster::mask(img_kruger_mean_april, shp_kruger)
img_kruger_mean_april <- img_kruger_mean_april * 0.02 - 273.15

# Make bounding box:
bbox_ele <- raster::bbox(raster::extent(ele_april_2009_spdf) * 1.2)
# Look whever the bbox of elephant is smaller or bigger than the one of weather data:
if(bbox_ele[1] < bbox(flux_all_sub_aq_night_spdf)[1]){print('elephant xmin is smaller, keep it')}
if(bbox_ele[2] < bbox(flux_all_sub_aq_night_spdf)[2]){print('elephant xmin is smaller, keep it')}else{print('Ibutton xmax is bigger, increase the bbox'); bbox_ele[2] <- bbox(flux_all_sub_aq_night_spdf)[2] }
if(bbox_ele[3] > bbox(flux_all_sub_aq_night_spdf)[3]){print('elephant ymin is smaller, keep it')}else{print('Ibutton ymin is bigger, increase the bbox'); bbox_ele[3] <- bbox(flux_all_sub_aq_night_spdf)[3] }
if(bbox_ele[4] > bbox(flux_all_sub_aq_night_spdf)[4]){print('elephant ymax is smaller, keep it')}else{print('Ibutton ymax is bigger, increase the bbox'); bbox_ele[4] <- bbox(flux_all_sub_aq_night_spdf)[4] }

img_kruger_mean_april_crop_by_bbox <- raster::crop(img_kruger_mean_april, bbox_ele)


  # MODIS in Kruger:
  pdf(file = paste0(outdir, '_MODIS_ele_april_2008_DES.pdf'))
  raster::plot(shp_kruger, col = alpha('gray49', 0.8) )
  raster::plot(img_kruger_mean_april_crop_by_bbox, zlim=range(raster::values(img_kruger_mean_april_crop_by_bbox$LST_Night_1km),
                                                              na.rm=T),
               col=cm.cols1(100),las=1,bty='n', add=T)
  raster::plot(raster::extent(img_kruger_mean_april_crop_by_bbox), col = 'black', add=T, lwd = 3)
  raster::scalebar(50, xy=c(32.02416, -22.37174), type='bar', divs=4)
  dev.off()
  # plot(extent(ele_april_2009_spdf) * 1.2, col = 'black', add=T, lwd = 3)
  # Flux Tower in Kruger:
  pdf(file = paste0(outdir, '_flux_ele_april_2008_DES.pdf'))
  raster::plot(shp_kruger, col = alpha('gray49', 0.8) )
  sp::plot(buf_sp, add= T, pch = 16,  cex = 1.3, col = alpha('gray49', 0.9))
  sp::plot(flux_all_sub_aq_night_spdf, add= T, pch = 16, col = alpha('black', 0.8))
  raster::plot(raster::extent(img_kruger_mean_april_crop_by_bbox), col = 'black', add=T, lwd = 3)
  raster::scalebar(50, xy=c(32.02416, -22.37174), type='bar', divs=4)
  GISTools::north.arrow(xb=32.20813, yb=-22.74687, len=0.05, lab="N")  
  dev.off()
  
  # Now Elephant
  pdf(file = paste0(outdir, 'tracks_ele_april_2008_points_DES.pdf'))
  raster::plot(shp_kruger, col = alpha('gray49', 0.8) )
  points(ele_april_2009_spdf, pch = 16, cex = 0.3, col = alpha('black', 0.8))
  scalebar(50, xy=c(32.02416, -22.37174), type='bar', divs=4)
  north.arrow(xb=32.20813, yb=-22.74687, len=0.05, lab="N")  
  plot(extent(img_kruger_mean_april_crop_by_bbox), col = 'black', add=T, lwd = 3)
  dev.off()
  
  # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
  pdf(file = paste0(outdir, '_elephant_case.pdf'))
  par(mfrow=c(1,3))
  #MODIS
  raster::plot(shp_kruger, col = alpha('gray49', 0.8) )
  raster::plot(img_kruger_mean_april_crop_by_bbox, zlim=range(raster::values(img_kruger_mean_april_crop_by_bbox$LST_Night_1km), na.rm=T),col=cm.cols1(100),las=1,bty='n', add=T)
  raster::plot(raster::extent(img_kruger_mean_april_crop_by_bbox), col = 'black', add=T, lwd = 3)
  raster::scalebar(50, xy=c(32.0240, -22.37174), type='bar', divs=4)
  GISTools::north.arrow(xb=32.20813, yb=-22.74687, len=0.05, lab="N")  
  # FLUX
  raster::plot(shp_kruger, col = alpha('gray49', 0.8) )
  sp::plot(buf_sp, add= T, pch = 16,  cex = 1.3, col = alpha('gray49', 0.9))
  sp::plot(flux_all_sub_aq_night_spdf, add= T, pch = 16, col = alpha('black', 0.8))
  raster::plot(raster::extent(img_kruger_mean_april_crop_by_bbox), col = 'black', add=T, lwd = 3)
  # Elephant
  raster::plot(shp_kruger, col = alpha('gray49', 0.8) )
  points(ele_april_2009_spdf, pch = 16, cex = 0.3, col = alpha('black', 0.8))
  scalebar(50, xy=c(32.02416, -22.37174), type='bar', divs=4)
  raster::plot(raster::extent(img_kruger_mean_april_crop_by_bbox), col = 'black', add=T, lwd = 3)
  
  north.arrow(xb=32.20813, yb=-22.74687, len=0.05, lab="N")  
  plot(extent(img_kruger_mean_april_crop_by_bbox), col = 'black', add=T, lwd = 3)

  
  # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
  flux_aq_night <- read.csv('/Users/diegoellis/Desktop/Submit_NCC_ABI/R4/Data/Elephant/SkukuzaTemp&Rainfall2007-2009.csv')
  flux_aq_night$tiempo <- mdy_hm(flux_aq_night$Timestamp)
  flux_aq_night$hour <- hour(flux_aq_night$tiempo)
  flux_aq_night$Air_temp <- flux_aq_night$Temperature.Ambient..deg.C.
  # elephant_thermal_env_dat <- read.csv('/Users/diegoellis/Downloads/ThermochronTracking Elephants Kruger 2007-1269689165289715956/ThermochronTracking Elephants Kruger 2007-1269689165289715956.csv')
  elephant_thermal_env_dat <- read.csv('/Users/diegoellis/Desktop/Submit_NCC_ABI/R4/Data/Elephant/ThermochronTracking Elephants Kruger 2007-6340748537733918440.csv')
  Aqua_night_flyovertime <- range(elephant_thermal_env_dat$MODIS.Land.Surface.Temp...Emissivity.1km.Daily.Aqua.View.Time.Night, na.rm = T)
  flux_aq_night_sub <- flux_aq_night[flux_aq_night$hour >= Aqua_night_flyovertime[1] & flux_aq_night$hour <= Aqua_night_flyovertime[2],] %>% tidyr::drop_na(Temperature.Ambient..deg.C.)
  flux_aq_night_sub$jday <- yday(flux_aq_night_sub$tiempo)
  
  # Subset   
  tmp_aq_n_ele_T = ggplot(ele_aq_night, aes(x = jday, y = Temp)) + geom_smooth() + geom_point(alpha = 1/20) + ylab('Elephant \n collected temperature') + theme_classic() + ylim(0, 30) + theme(text = element_text(size=20))
  tmp_aq_n_ele_aq_T = ggplot(ele_aq_night, aes(x = jday, y = MODIS_Aqua_Daily_1km_night_LST)) + geom_smooth() + geom_point(alpha = 1/20) +
    theme_classic() + ylab('Aqua Night \n Land Surface temperature') + ylim(0, 30) + theme(text = element_text(size=20))
  tmp_aq_night_flux_T <- ggplot(flux_aq_night_sub,aes(x = jday, y = Temperature.Ambient..deg.C.))  + geom_point(alpha = 1/20)  + geom_smooth() + theme_classic()+ylab('Flux tower \n collected temperature') + theme_classic() + ylim(0, 30) + theme(text = element_text(size=20))
  pdf(paste0(outdir, '/elephant_inset_time.pdf'), height = 12, width = 8)
  grid.arrange(tmp_aq_n_ele_aq_T, tmp_aq_night_flux_T, tmp_aq_n_ele_T, ncol = 1)
  dev.off()
  
  
  
  
  ele_aq_night$year <- year(ele_aq_night$timestamp)
  ele_aq_night_averaged = plyr::ddply(ele_aq_night, 'year', function(x){
    plyr::ddply(x, 'jday', function(y){
      data.frame(
        u_T = mean(y$Temp, na.rm=T),
        u_LST = mean(y$MODIS_Aqua_Daily_1km_night_LST, na.rm=T),
        jday = unique(y$jday),
        year = unique(y$year),
        jday_year = paste0(unique(y$year), '-', unique(y$jday))
      )
    })
  })
  
  moving_window_T_LST = ele_aq_night_averaged  %>% 
    dplyr::mutate(T_ele_03da = zoo::rollmean(u_T, k = 3, fill = NA),
                  T_ele_05da = zoo::rollmean(u_T, k = 5, fill = NA),
                  T_ele_07da = zoo::rollmean(u_T, k = 7, fill = NA),
                  T_ele_15da = zoo::rollmean(u_T, k = 15, fill = NA),
                  T_ele_21da = zoo::rollmean(u_T, k = 21, fill = NA),
                  
                  LST_03da = zoo::rollmean(u_LST, k = 3, fill = NA),
                  LST_05da = zoo::rollmean(u_LST, k = 5, fill = NA),
                  LST_07da = zoo::rollmean(u_LST, k = 7, fill = NA),
                  LST_15da = zoo::rollmean(u_LST, k = 15, fill = NA),
                  LST_21da = zoo::rollmean(u_LST, k = 21, fill = NA)
    ) %>% 
    dplyr::ungroup()
  
  
  
  moving_window_T_LST_2009 = moving_window_T_LST %>% filter(year == 2009)
  
  moving_window_T_LST_2009 %>% filter(T_ele_21da == min(moving_window_T_LST_2009$T_ele_21da, na.rm=T)) %>% dplyr::select(jday_year)
  
  moving_window_T_LST_2009 %>% filter(LST_21da == min(moving_window_T_LST_2009$LST_21da, na.rm=T)) %>% dplyr::select(jday_year)
  
  moving_window_T_LST_2009 %>% filter(T_ele_21da == max(moving_window_T_LST_2009$T_ele_21da, na.rm=T)) %>% dplyr::select(jday_year)
  
  moving_window_T_LST_2009 %>% filter(LST_21da == max(moving_window_T_LST_2009$LST_21da, na.rm=T)) %>% dplyr::select(jday_year)
  
  
  # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
  # Fluxtower ####
  # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
  
  flux_all <- read.table('/Users/diegoellis/projects/Phd_Chapter_1/Data/Fluxtower_Skukuza/FLX_ZA-Kru_FLUXNET2015_FULLSET_2000-2013_1-4/FLX_ZA-Kru_FLUXNET2015_ERAI_HH_1989-2014_1-4.csv', sep =',', header=T)
  
  flux_all$TIMESTAMP_START <- as.character(flux_all$TIMESTAMP_START)
  flux_all$TIMESTAMP_START <- ymd_hm(flux_all$TIMESTAMP_START)
  flux_all$TIMESTAMP_END <- as.character(flux_all$TIMESTAMP_END)
  flux_all$TIMESTAMP_END <- ymd_hm(flux_all$TIMESTAMP_END)
  flux_all$year <- lubridate::year(flux_all$TIMESTAMP_END)
  flux_all_sub <- flux_all[flux_all$year >= 2007 & flux_all$year <2010,  ]
  flux_all_sub$hour <- lubridate::hour(flux_all_sub$TIMESTAMP_START)
  flux_all_sub$jday <- lubridate::yday(flux_all_sub$TIMESTAMP_START)
  # Subset to flyover time of aqua night:
  flux_all_sub_aq_night <- flux_all_sub[flux_all_sub$hour >= Aqua_night_flyovertime[1] & flux_all_sub$hour <= Aqua_night_flyovertime[2],] %>% drop_na(TA_ERA)
  
  flux_aq_night_averaged = plyr::ddply(flux_all_sub_aq_night, 'year', function(x){
    plyr::ddply(x, 'jday', function(y){
      data.frame(
        u_T = mean(y$TA_ERA, na.rm=T),
        jday = unique(y$jday),
        year = unique(y$year),
        jday_year = paste0(unique(y$year), '-', unique(y$jday))
      )
    })
  })
  
  moving_window_T_flux = flux_aq_night_averaged  %>% 
    dplyr::mutate(T_flux_03da = zoo::rollmean(u_T, k = 3, fill = NA),
                  T_flux_05da = zoo::rollmean(u_T, k = 5, fill = NA),
                  T_flux_07da = zoo::rollmean(u_T, k = 7, fill = NA),
                  T_flux_15da = zoo::rollmean(u_T, k = 15, fill = NA),
                  T_flux_21da = zoo::rollmean(u_T, k = 21, fill = NA)
                  
    ) %>% 
    dplyr::ungroup()
  
  moving_window_T_flux %>%  dplyr::filter(year == 2009)  %>% 
    tidyr::pivot_longer(names_to = "rolling_mean_key", 
                        values_to = "rolling_mean_value", 
                        cols = c(u_T, 
                                 T_flux_03da, 
                                 T_flux_21da)) %>% 
    ggplot2::ggplot(aes(x = jday, 
                        y = rolling_mean_value, 
                        color = rolling_mean_key)) +
    ggplot2::geom_line() +   
    ggplot2::labs(title = "Fluxtower temperature rolling average", 
                  subtitle = "For the year 2009",
                  y = "Fluxtower temperature", 
                  color = "Metric",
                  x = "jday_year") + 
    theme_classic()
  
  moving_window_T_flux_2009 = moving_window_T_flux %>% filter(year == 2009)
  
  # 21 days
  moving_window_T_flux_2009 %>% filter(T_flux_21da == min(moving_window_T_flux_2009$T_flux_21da, na.rm=T)) %>% dplyr::select(jday_year)
  # 208
  moving_window_T_flux_2009 %>% filter(T_flux_21da == max(moving_window_T_flux_2009$T_flux_21da, na.rm=T)) %>% dplyr::select(jday_year)
  # 18
  