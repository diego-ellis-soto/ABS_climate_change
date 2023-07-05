# --- --- --- --- --- --- --- --- --- --- --- ---
#
# Animal-borne sensors as a biologically informed lens on a changing climate
#
# Diego Ellis Soto, diego.ellissoto@yale.edu
#
#
# Generate Figure 1:
# Create a map showcasing a collection of studies of Animal Borne Sensors across the world in land, air and sea
# --- --- --- --- --- --- --- --- --- --- --- ---



require(rnaturalearthhires);require(palettetown);require(plyr);library(ggplot2);library(rgdal);require(tidyverse);library(mapproj);require(sp);require(ggthemes);library(sf);library(raster);library(dplyr);library(spData);library(spDataLarge);library(rworldmap);library(geosphere);library(gpclib)


load('/Users/diegoellis/Desktop/Submit_NCC_ABI/R4/Data/2023_centroid_locs_biolog_studies.Rdata')

load(url("https://github.com/valentinitnelav/RandomScripts/blob/master/NaturalEarth.RData?raw=true")) # This dude has a nice worldmap I think


crs = "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
world <- fortify(spTransform(getMap(), CRS(paste0(crs))))
NE_box_proj        <- spTransform(NE_box, CRSobj = paste0(crs))
NE_graticules_proj <- spTransform(NE_graticules, CRSobj = paste0(crs))
centroid_locations <- SpatialPointsDataFrame(coords = centroid_locations[,c('Longitude','Latitude')], data = centroid_locations, proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +towgs84=0,0,0"))
centroid_locations_proj <- spTransform(centroid_locations, CRSobj = paste0(crs))
centroid_locations_proj_df <- centroid_locations_proj %>% as.data.frame() %>% mutate(Longitude = Longitude.1,Latitude = Latitude.1)  %>% dplyr::select(Longitude, Latitude, Locomotion)
prj.coord <- project(cbind(lbl.Y$lon, lbl.Y$lat), proj=paste0(crs))
lbl.Y.prj <- cbind(prj.coord, lbl.Y)
names(lbl.Y.prj)[1:2] <- c("X.prj","Y.prj")
prj.coord <- project(cbind(lbl.X$lon, lbl.X$lat), proj=paste0(crs))
lbl.X.prj <- cbind(prj.coord, lbl.X)
names(lbl.X.prj)[1:2] <- c("X.prj","Y.prj")
centroid_locations_proj_df$Locomotion <- as.factor(centroid_locations_proj_df$Locomotion)


centroid_locations_proj_df$Locomotion <- ordered(centroid_locations_proj_df$Locomotion, levels = c("Terrestrial", "Aerial", "Aquatic"))

ggplot() + 
  geom_map(data=world, map=world,
                    aes(x=long, y=lat, map_id=id),
                    color="gray50", fill="gray50", size=0.25)+
  geom_polygon(data=NE_box_proj, aes(x=long, y=lat),colour="black",fill="transparent", size = 0.25)+
  coord_equal() + 
  theme_map() +
  geom_point(data = centroid_locations_proj_df, aes(x = Longitude, y = Latitude, colour = factor(Locomotion)) , alpha = 0.5, size = 4)+
  geom_path(data=NE_graticules_proj, aes(long, lat, group=group), linetype="dotted", color="grey50", size = 0.25)+ scale_color_manual(values = c('Aerial' = "#E69F00",'Aquatic' =  "#56B4E9",'Terrestrial' =  "#009E73") )
