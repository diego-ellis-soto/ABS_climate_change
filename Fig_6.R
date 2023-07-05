# --- --- --- --- --- --- --- --- --- --- --- ---
#
# Animal-borne sensors as a biologically informed lens on a changing climate
#
# Diego Ellis Soto, diego.ellissoto@yale.edu
#
# Figure 6: 
#
# Clean code for FPT and recurse for Albatross and white stork: 
#
#
#
# --- --- --- --- --- --- --- --- --- --- --- ---

require(hddtools)
require(scales)
require(lubridate)
require(move)
library(adehabitatLT)
library(raster)
library(spatial)
library(rgdal)
library(dplyr)
library(ctmm)
library(recurse)
require(ggplot2)
require(viridis)
require(ggmap)
require(recurse)  
require(moveHMM)
require(sf)
require(ggspatial)
library(raster)
require(adehabitatHR)
library(GISTools)  

getContinuousPalette = function(n, alpha = 1)
{
  cols = alpha(brewer_pal(palette = "RdYlBu", direction = -1)(9), alpha)
  return( gradient_n_pal(cols)(seq(0, 1, length = n)) )
}


# Albatross: ####
albatross <- move('/Users/diegoellis/projects/development/Animove_2019/ENSO/Movement_data/Galapagos Albatrosses.csv')
albatross_df = as.data.frame(albatross)
move_object = albatross[[3]]
move_object = albatross[['X4261.2228']]
move_object$year <- year(move_object$timestamp)
move_object$hour <- hour(move_object$timestamp)
move_object$id_year <- paste0(move_object@idData$individual.local.identifier, '_', year(move_object$timestamp))
move_object_df = as.data.frame(move_object)
PROJ = paste0("+proj=laea +lat_0=", round( mean(move_object@coords[,2]) ), ' +lon_0=', round( mean(move_object@coords[,1]) ), ' +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0')

move_obj_laea <- spTransform(move_object, PROJ)

move::plot(move_object, col = viridis_pal()(nrow(move_object)), pch = 20, 
           xlab = "Longitude", ylab = "Latitude", asp = 1, main = paste0('Animal track individual ', unique(move_object@idData$individual.local.identifier), ' ',
                                                                         unique(move_object@idData$individual.taxon.canonical.name)))
move_object_df_sp <- SpatialPointsDataFrame(coords = move_object_df[,c('location.long', 'location.lat')], data = move_object_df, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
proj4string(move_object_df_sp)

# Calculate are in km2:
st_as_sf(move_obj_laea) %>% st_bbox %>% st_as_sfc %>% st_area()
area_albatros = st_as_sf(move_obj_laea) %>% st_bbox %>% st_as_sfc %>% st_area()
units::set_units(area_albatros, km^2)

# In lambert equal area for FPT:
move_obj_ltraj <- as.ltraj(xy=move_object@coords, date=move_object$timestamp, id=move_object$id_year,typeII = TRUE) #bursts are animal yearsalbatross_1$id_year 
move_obj_ltraj_laea <- as.ltraj(xy=move_obj_laea@coords, date=move_obj_laea$timestamp, id=move_obj_laea$id_year,typeII = TRUE) #bursts are animal yearsalbatross_1$id_year 
#determine fix schedule
# First passage time:
hist(move_obj_ltraj_laea[[1]]$dt/3600, xlim=c(0,24), breaks=100) #schedule changes seasonally and among ids, 
#fix rates typically 1,2,4 hours
refda <- min(move_obj_ltraj_laea[[1]]$date)
wost_NA <- setNA(move_obj_ltraj_laea,refda,1,units="hour")
# Set variogram:
wost_demo <- sett0(wost_NA,refda,1,units="hour")
is.regular(wost_demo)
wost_fpt1 <- fpt(wost_demo, radii=1:100, units="hours")
varlogfpt(wost_fpt1, graph=TRUE)
title('In lambert equal area', outer=TRUE)
refda <- min(move_obj_ltraj[[1]]$date)
wost_NA <- setNA(move_obj_ltraj,refda,1,units="hour")
# Set variogram:
wost_demo <- sett0(wost_NA,refda,1,units="hour")
is.regular(wost_demo)
wost_fpt1 <- fpt(wost_demo, radii=1:12, units="hours")

varlogfpt(wost_fpt1, graph=TRUE)
title('Long-Lat', outer=TRUE)
# plot(wost_fpt1, scale = 500, warn=FALSE)

# Recurse:  
move_object_df = as.data.frame(move_object) %>% dplyr::mutate(x =location.long,
                                                              y = location.lat,
                                                              t = timestamp,
                                                              id = individual.local.identifier) %>%   dplyr::select(x, y, t, id)

move_object_df_revisit = getRecursions(move_object_df, 0.5)  # Defined by first pssage time
albatross_df_revisit = move_object_df_revisit


albatross_df = as.data.frame(move_object) %>% dplyr::mutate(x =location.long,
                                                            y = location.lat,
                                                            t = timestamp,
                                                            id = individual.local.identifier) %>%   dplyr::select(x, y, t, id)
# Plot the revisitations:
par(mfrow=c(1,1))
plot(albatross_df_revisit, albatross_df, legendPos = c(-84, -8))
drawCircle(-89, -8, 0.5)

colours = unique(getContinuousPalette(max(move_object_df_revisit$revisits), 0.5)[move_object_df_revisit$revisits])

hist(move_object_df_revisit$revisits, breaks = 20, main = "", xlab = "Revisits (radius = 2)", 
     col = colours)
# Average and var and sd of revisitation:
print('Calculate average, variance, st of revisitation')
print(paste0('Average revisitation is ', round(mean(move_object_df_revisit$revisits), 2),
             '; sd is ', round(sd(move_object_df_revisit$revisits),2),
             '; variance is ', round(var(move_object_df_revisit$revisits),2)))

# examine first passage time
# the first visits passes through the center of the circle, thus representing a first passage time
# the first visit at each location equals the first passage time
# this only works in teh recurse package if you use single individuals, not available for multiple individuals
move_object_df_revisit$firstPassageTime = move_object_df_revisit$revisitStats$timeInside[move_object_df_revisit$revisitStats$visitIdx==1]

# Make a map #####
move_object_sp = as(move_object, 'Spatial')

map.albatross <- qmap( ( bbox((move_object_sp))  * 1.01 ), zoom = 9, maptype = 'hybrid', legend="topright") # try a smaller zoom
map.albatross <- qmap( ( bbox((move_object_sp))   ), zoom = 9, maptype = 'hybrid', legend="topright") # try a smaller zoom
# map.albatross <- qmap( bbox((move_object_sp)) , zoom = 9, maptype = 'hybrid', legend="topright", extent = 'normal') # try a smaller


print(map.albatross + 
        geom_path(data = albatross_df, aes(x = x, y = y), color = "white", size = 0.3) +
        geom_spatial_point(data = albatross_df, aes(x = x, y = y), 
                           color = getContinuousPalette(max(albatross_df_revisit$revisits), 0.5)[albatross_df_revisit$revisits],
                           crs = 4326) + 
        coord_sf(crs = 4326) +
        annotation_scale() 
      
)

cutOff = 100
# blue larger less than 100 hours, these areas are very revisited. 
# Areas used a lot by animals.
print(map.albatross + 
        geom_path(data = move_object_df, aes(x = x, y = y),color = "white", size = 0.3) + 
        geom_point(data = move_object_df, aes(x = x, y = y), 
                   color=alpha(ifelse(move_object_df_revisit$firstPassageTime >cutOff, "blue", "grey"),.5)))

temp_NA<-setNA(ltraj=move_obj_ltraj, date.ref=min(move_obj_ltraj[[1]]$date), dt=1, units="hour")
temp_NA<-sett0(temp_NA, date.ref=min(move_obj_ltraj[[1]]$date), dt=1, units="hour")

# Calculate HR size of an albatross: 
radii<-round(seq(1,100, length.out = 100)) #make upper limit equal to HR size
# Needs to be aseries of radius so you get first passage time from 10 meters to 3000m 
# Look for maximum semivariance over those different size radii, scale on x axis and variance of log of FPT is on y axis

tempFPT<-fpt(temp_NA, radii=1:100, units = 'hours')

var<-varlogfpt(tempFPT, graph=T)

# adfFPT<-cbind(unique(infolocs(move_obj_ltraj)[[1]]$id_year), var)

# head(infolocs(move_obj_ltraj)[[1]])
# }

PROJ = paste0("+proj=laea +lat_0=", round( mean(move_object@coords[,2]) ), ' +lon_0=', round( mean(move_object@coords[,1]) ), ' +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0')

df_albatross_laea <- spTransform(move_object, PROJ)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# White storks ####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

loburg <- read.csv('/Users/diegoellis/Downloads/loburg_2014.csv') %>% 
  dplyr::mutate(timestamp = ymd_hms(timestamp),
         year = year(timestamp),
         hour = hour(timestamp),
         id_year = niche_name)

df_stork<-   move::move(x=loburg$lon,
                        y=loburg$lat,
                        time=loburg$timestamp,
                        animal=loburg$niche_name,
                        removeDuplicatedTimestamps=T)

proj4string(df_stork) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

PROJ = paste0("+proj=laea +lat_0=", round( mean(df_stork@coords[,2]) ), ' +lon_0=', round( mean(df_stork@coords[,1]) ), ' +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0')

df_stork_laea <- spTransform(df_stork, PROJ)

# Calculate are in km2:
st_as_sf(df_stork_laea[['Magnus.2014']]) %>% st_bbox %>% st_as_sfc %>% st_area()
area_stork = st_as_sf(df_stork_laea[['Magnus.2014']]) %>% st_bbox %>% st_as_sfc %>% st_area()
units::set_units(area_stork , km^2)

# In lambert equal area for FPT:
move_obj_ltraj_laea <- as.ltraj(xy=df_stork_laea@coords, date=df_stork_laea$time, id=df_stork_laea@trackId,typeII = TRUE) #bursts are animal yearsalbatross_1$id_year 

magnus_2014 <- loburg %>% filter(niche_name == 'Magnus-2014') %>% dplyr::select(lon, lat, niche_name, dist2nest) %>% 
  dplyr::rename(ID = niche_name)

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# First passage time:
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

hist(move_obj_ltraj_laea[[1]]$dt/3600, xlim=c(0,24), breaks=100) #schedule changes seasonally and among ids, 
#fix rates typically 1,2,4 hours
refda <- min(move_obj_ltraj_laea[[1]]$date)
wost_NA <- setNA(move_obj_ltraj_laea,refda,1,units="hour")
# Set variogram:
wost_demo <- sett0(wost_NA,refda,1,units="hour")
is.regular(wost_demo)
wost_fpt1 <- fpt(wost_demo, radii=1:10000, units="hours")
varlogfpt(wost_fpt1, graph=TRUE)
title('In lambert equal area', outer=TRUE)
refda <- min(move_obj_ltraj[[1]]$date)
wost_NA <- setNA(move_obj_ltraj,refda,1,units="hour")

# Set variogram:
wost_demo <- sett0(wost_NA,refda,1,units="hour")
is.regular(wost_demo)
wost_fpt1 <- fpt(wost_demo, radii=1:12, units="hours")

varlogfpt(wost_fpt1, graph=TRUE)
title('Long-Lat', outer=TRUE)
plot(wost_fpt1, scale = 500, warn=FALSE)
# Just Magnus:

refda_magnus <- min(move_obj_ltraj_laea[[8]]$date)
wost_NA_magnus <- setNA(move_obj_ltraj_laea[8],refda,1,units="hour")
# Set variogram:
wost_demo_magnus <- sett0(wost_NA_magnus,refda,1,units="hour")
is.regular(wost_demo_magnus)
wost_fpt1_magnus <- fpt(wost_demo_magnus, radii=1:10000, units="hours")
varlogfpt(wost_fpt1_magnus, graph=TRUE)
title('In lambert equal area', outer=TRUE)

# Recurse 

move_object_df = loburg %>% dplyr::filter(niche_name == 'Magnus-2014') %>% 
  dplyr::select(lon, lat, niche_name, timestamp) %>% 
  dplyr::mutate(x =lon,
                y = lat,
                t = timestamp,
                id = niche_name) %>% 
  dplyr::select(x, y, t, id)


move_object_df_revisit = getRecursions(move_object_df, 0.05)  # Defined by first pssage time


par(mfrow = c(1, 2), mar = c(4, 4, 1, 1))
plot(move_object_df_revisit, move_object_df, legendPos = c(mean(move_object_df$x), min(move_object_df$y)))
drawCircle(mean(move_object_df$x), min(move_object_df$y), 0.05)

hist(move_object_df_revisit$revisits, breaks = 20, main = "", xlab = "Revisits (radius = 0.05)")
summary(move_object_df_revisit$revisits)

print(paste0('Average revisitation is ', round(mean(move_object_df_revisit$revisits), 2),
             '; sd is ', round(sd(move_object_df_revisit$revisits),2),
             '; variance is ', round(var(move_object_df_revisit$revisits),2)))

move_object_df_revisit$firstPassageTime = move_object_df_revisit$revisitStats$timeInside[move_object_df_revisit$revisitStats$visitIdx==1]

hist(as.numeric(move_object_df_revisit$firstPassageTime), breaks = 20, col = "darkgrey", border = NA, main = "", xlab = "First passage (hrs)") # histogram of the first 

# Make a map #####

magnus_2014_sp <- SpatialPointsDataFrame(coords= magnus_2014[,c('lon', 'lat')],data = magnus_2014,
                                         proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))



# map.move <- qmap(bbox(extent(move_object_sp)), zoom = 9, maptype = 'hybrid', legend="topright") # try a smaller zoom

map.move <- qmap(bbox((magnus_2014_sp)), zoom = 13, maptype = 'hybrid', legend="topright") # try a smaller zoom

print(map.move + 
        geom_path(data = magnus_2014, aes(x = lon, y = lat), color = "white", size = 0.3) + 
        geom_point(data = magnus_2014, aes(x = lon, y = lat), 
                   color = getContinuousPalette(max(move_object_df_revisit$revisits), 0.5)[move_object_df_revisit$revisits]))



print(map.move + 
        geom_path(data = magnus_2014, aes(x = lon, y = lat), color = "white", size = 0.3) + 
        geom_spatial_point(data = magnus_2014, aes(x = lon, y = lat), 
                           color = getContinuousPalette(max(move_object_df_revisit$revisits), 0.5)[move_object_df_revisit$revisits],
                           crs = 4326)) +  coord_sf(crs = 4326) + annotation_scale()


cutOff = 200


print(map.move + 
        geom_path(data = magnus_2014, aes(x = lon, y = lat),color = "white", size = 0.3) + 
        geom_point(data = magnus_2014, aes(x = lon, y = lat), 
                   color=alpha(ifelse(move_object_df_revisit$firstPassageTime >cutOff, "blue", "grey"),.5)))


# Recurse
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Caribou                                                   #####
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

caribou <- move('/Users/diegoellis/Desktop/Submit_NCC_ABI/R4/Data/Albatross_caribou_stork/Mountain caribou in British Columbia-gps.csv')
caribou_df = as.data.frame(caribou)

# sort(table(caribou@trackId))
move_object = caribou[['BP_car145']]
move_object$year <- year(move_object$timestamp)
move_object$hour <- hour(move_object$timestamp)
move_object$id_year <- paste0(move_object@idData$individual.local.identifier, '_', year(move_object$timestamp))
move_object_df = caribou_df


PROJ = paste0("+proj=laea +lat_0=", round( mean(move_object@coords[,2]) ), ' +lon_0=', round( mean(move_object@coords[,1]) ), ' +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0')

move_obj_laea <- spTransform(move_object, PROJ)

move::plot(move_object, col = viridis_pal()(nrow(move_object)), pch = 20, 
           xlab = "Longitude", ylab = "Latitude", asp = 1, main = paste0('Animal track individual ', unique(move_object@idData$individual.local.identifier), ' ',
                                                                         unique(move_object@idData$individual.taxon.canonical.name)))

move_object_df_sp <- SpatialPointsDataFrame(coords = move_object_df[,c('location.long', 'location.lat')], data = move_object_df, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
proj4string(move_object_df_sp)

# Calculate are in km2:
st_as_sf(move_obj_laea) %>% st_bbox %>% st_as_sfc %>% st_area()
area_caribou = st_as_sf(move_obj_laea) %>% st_bbox %>% st_as_sfc %>% st_area()
units::set_units(area_caribou, km^2) # 4204.297 [km^2]

# In lambert equal area for FPT:
move_obj_ltraj <- as.ltraj(xy=move_object@coords, date=move_object$timestamp, id=move_object$id_year,typeII = TRUE) #bursts are animal yearsalbatross_1$id_year 
move_obj_ltraj_laea <- as.ltraj(xy=move_obj_laea@coords, date=move_obj_laea$timestamp, id=move_obj_laea$id_year,typeII = TRUE) #bursts are animal yearsalbatross_1$id_year 
#determine fix schedule
# First passage time:
hist(move_obj_ltraj_laea[[1]]$dt/3600, xlim=c(0,24), breaks=100) #schedule changes seasonally and among ids, 
summary(move_obj_ltraj_laea[[1]]$dt/3600)
#fix rates typically 1,2,4 hours
refda <- min(move_obj_ltraj_laea[[1]]$date)
wost_NA <- setNA(move_obj_ltraj_laea,refda,7,units="hour")
# Set variogram:
wost_demo <- sett0(wost_NA,refda,7,units="hour")
is.regular(wost_demo)
wost_fpt1 <- fpt(wost_demo, radii=1:100, units="hours")
wost_fpt1 <- fpt(wost_demo, radii=0.1:2, units="hours")
varlogfpt(wost_fpt1, graph=TRUE)
title('In lambert equal area', outer=TRUE)
refda <- min(move_obj_ltraj[[1]]$date)
wost_NA <- setNA(move_obj_ltraj,refda,1,units="hour")
# Set variogram:
wost_demo <- sett0(wost_NA,refda,1,units="hour")
is.regular(wost_demo)
wost_fpt1 <- fpt(wost_demo, radii=1:12, units="hours")

varlogfpt(wost_fpt1, graph=TRUE)
title('Long-Lat', outer=TRUE)
plot(wost_fpt1, scale = 500, warn=FALSE)

# Recurse:  
move_object_df = as.data.frame(move_object) %>% 
  dplyr::mutate(x =location.long,
                y = location.lat,
                t = timestamp,
                id = individual.local.identifier) %>% 
  dplyr::select(x, y, t, id)

move_object_df_revisit = getRecursions(move_object_df, 0.01)  # Defined by first pssage time

caribou_df_revisit = move_object_df_revisit

caribou_df = as.data.frame(move_object) %>% dplyr::mutate(x =location.long,
                                                          y = location.lat,
                                                          t = timestamp,
                                                          id = individual.local.identifier) %>%   select(x, y, t, id)


# par(mfrow = c(1, 2), mar = c(4, 4, 1, 1))
par(mfrow=c(1,1))
plot(caribou_df_revisit, caribou_df, legendPos = c(-123.0, 55.2))
drawCircle(-122.6, 55.2, 0.01)

colours = unique(getContinuousPalette(max(move_object_df_revisit$revisits), 0.5)[move_object_df_revisit$revisits])

hist(move_object_df_revisit$revisits, breaks = 20, main = "", xlab = "Revisits (radius = 2)", 
     col = colours)
# Average and var and sd of revisitation:
print('Calculate average, variance, st of revisitation')
summary(move_object_df_revisit$revisits)
var(move_object_df_revisit$revisits)
# Caribou info ####
# Caribou: Temporal grain: 7 hours, Scale of seleciton 2 degrees? Domain: 4204.297 [km^2], revisitation: mean 13.04, variance:  146

print(paste0('Average revisitation is ', round(mean(move_object_df_revisit$revisits), 2),
             '; sd is ', round(sd(move_object_df_revisit$revisits),2),
             '; variance is ', round(var(move_object_df_revisit$revisits),2)))

# the first visits passes through the center of the circle, thus representing a first passage time
# the first visit at each location equals the first passage time
# this only works in teh recurse package if you use single individuals, not available for multiple individuals
move_object_df_revisit$firstPassageTime = move_object_df_revisit$revisitStats$timeInside[move_object_df_revisit$revisitStats$visitIdx==1]

# Make a map #####
move_object_sp = as(move_object, 'Spatial')


Area_Map <- ggmap(get_map(location = bbox((move_object_sp)),
                          zoom = 4,
                          maptype = "satellite",
                          color = "color"
)
)


# map.move <- qmap(bbox(extent(move_object_sp)), zoom = 9, maptype = 'hybrid', legend="topright") # try a smaller zoom
move_object$longitude = move_object@coords[,1]
move_object$latitude = move_object@coords[,2]

map.caribou <- qmap( ( bbox((move_object_sp))  * 1.01 ), zoom = 9, maptype = 'hybrid', legend="topright") # try a smaller zoom

# map.albatross <- qmap( bbox((move_object_sp)) , zoom = 9, maptype = 'hybrid', legend="topright", extent = 'normal') # try a smaller

# print(map.move + 
#         geom_path(data = albatross_df, aes(x = x, y = y), color = "white", size = 0.3) + 
#         geom_point(data = albatross_df, aes(x = x, y = y), 
#                    color = getContinuousPalette(max(albatross_df_revisit$revisits), 0.5)[albatross_df_revisit$revisits]))
# 

print(map.caribou + 
        geom_path(data = caribou_df, aes(x = x, y = y), color = "white", size = 0.3) +
        geom_spatial_point(data = caribou_df, aes(x = x, y = y), 
                           color = getContinuousPalette(max(caribou_df_revisit$revisits), 0.5)[caribou_df_revisit$revisits],
                           crs = 4326) + 
        coord_sf(crs = 4326) +
        annotation_scale() 
      
)
# This one works: 
qmap(( bbox(as(move_object, 'Spatial')) ) , zoom = 12, maptype = 'hybrid', legend="topright") + geom_path(aes(x = longitude, y = latitude), data = as.data.frame(move_object), color = "white", size = 0.3 ) + 
  geom_spatial_point(data = as.data.frame(move_object), aes(x = longitude, y = latitude), 
                     color = getContinuousPalette(max(caribou_df_revisit$revisits), 0.5)[caribou_df_revisit$revisits],
                     crs = 4326) + 
  coord_sf(crs = 4326) +
  annotation_scale() 


#
cutOff = 100
# blue larger less than 100 hours, these areas are very revisited. 
# Areas used a lot by animals.
print(map.caribou + 
        geom_path(data = move_object_df, aes(x = x, y = y),color = "white", size = 0.3) + 
        geom_point(data = move_object_df, aes(x = x, y = y), 
                   color=alpha(ifelse(move_object_df_revisit$firstPassageTime >cutOff, "blue", "grey"),.5)))

temp_NA<-setNA(ltraj=move_obj_ltraj, date.ref=min(move_obj_ltraj[[1]]$date), dt=1, units="hour")
temp_NA<-sett0(temp_NA, date.ref=min(move_obj_ltraj[[1]]$date), dt=1, units="hour")

# Calculate HR size of an albatross: 
radii<-round(seq(1,100, length.out = 100)) #make upper limit equal to HR size
# Needs to be aseries of radius so you get first passage time from 10 meters to 3000m 
# Look for maximum semivariance over those different size radii, scale on x axis and variance of log of FPT is on y axis

tempFPT<-fpt(temp_NA, radii=1:100, units = 'hours')

var<-varlogfpt(tempFPT, graph=T)

adfFPT<-cbind(unique(infolocs(move_obj_ltraj)[[1]]$id_year), var)

head(infolocs(move_obj_ltraj)[[1]])

#################################
