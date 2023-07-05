# --- --- --- --- --- --- --- --- --- --- --- ---
#
# Animal-borne sensors as a biologically informed lens on a changing climate
#
# Diego Ellis Soto, diego.ellissoto@yale.edu
#
# --- --- --- --- --- --- --- --- --- --- --- ---

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Figure 2 (a) Elephant seal collecting oceanographic measurements
#
# Code adapted from: https://github.com/aodn/imos-user-code-library/wiki/Using-the-IMOS-User-Code-Library-with-R
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

source('/Users/diegoellis/projects/Phd_Chapter_1/ABI_Phd_Chapter_1_code/0_source_packages_ABI.R')
source('/Users/diegoellis/projects/Phd_Chapter_1/Older_stuff/Meop/ncParse.R')

file_URL <- 'http://thredds.aodn.org.au/thredds/dodsC/IMOS/eMII/demos/AATAMS/marine_mammal_ctd-tag/2009_2011_ct64_Casey_Macquarie/ct64-M746-09/IMOS_AATAMS-SATTAG_TSP_20100205T043000Z_ct64-M746-09_END-20101029T071000Z_FV00.nc'
dataset <- ncParse( file_URL)

map <- maps::map( 'world', plot = F)
map2 <- map
map2$x <- map2$x +360
map$x <- c(map$x,map2$x)
map$y <- c(map$y,map2$y)

## Extract data from variables
temp <- dataset$variables$TEMP$data
psal <- dataset$variables$PSAL$data
pres <- dataset$variables$PRES$data
index <- dataset$variables$parentIndex$data
lat <- round( dataset$variables$LATITUDE$data, 1)
lon <- round( dataset$variables$LONGITUDE$data, 1)
date <- dataset$variables$TIME$data


## Extract data from global and variable attributes
title <- dataset$metadata$title
templab <- gsub( '_', ' ', dataset$variables$TEMP$long_name)
datelab <- gsub( '_', ' ', dataset$variables$TIME$standard_name)
psallab <- gsub( '_', ' ', dataset$variables$PSAL$long_name)
preslab <- gsub( '_', ' ', dataset$variables$PRES$long_name)
latlab <- dataset$variables$LATITUDE$long_name
lonlab <- dataset$variables$LONGITUDE$long_name
tempunit <- gsub( '_', ' ', dataset$variables$TEMP$units)
psalunit <- gsub( '_', ' ', dataset$variables$PSAL$units)
presunit <- gsub( '_', ' ', dataset$variables$PRES$units)

## Creation of a temperature vector
tempsq <- seq( round( min( temp), 0), round( max( temp), 0), 0.1)
tsq <- matrix( ncol = 1, nrow = length( temp))
for (i in 1:length( temp)){
  tsq[i] <- which.min( abs( tempsq- round( temp[i], 1)))
}

## Creation of a date vector
dates <- c()
for (i in as.numeric( rownames( table( index)))){
  dates <- c( dates, as.numeric(rep( date[i], length( which(( index) == i)))))
}

## Creation of a sea surface salinity vector
sssal <- c()
for (i in as.numeric( rownames( table( index)))){
  sssal <- c( sssal, head( psal[which(( index) ==i )], n = 1L))
}

## Creation of a salinity vector
salsq <- seq( round( min( sssal), 2), round( max( sssal), 2), 0.01)
ssq <- matrix( ncol = 1, nrow = length( sssal))
for (i in 1:length( sssal)){
  ssq[i] <- which.min( abs( salsq- round( sssal[i], 1)))
}

## The longitude in the original dataset ranges from [-180 to 180].
## To deal with the international date line issue, we change the range of longitude values to [0 to 360].
for (i in 1:length( lon)){
  if ( lon[i] < 0) {lon[i] = lon[i] + 360}
}

## To deal with the international date line issue, we change the range of the map longitude values to [0 to 360].

# Plot all the profiles of a time series
par( mar = c( 4.5, 4.5, 3, 4.5), cex=  1.2)
split.screen( c( 2, 1))
screen( 1)

plot( dates, -pres, col = colorRampPalette( viridis(5))( length( tempsq))[tsq], type='n',xaxt = 'n',
      xlab = datelab, ylab = paste( preslab, ' (negative dbar)', sep = ''), main = paste( dataset$metadata$species_name, ' - released in ', dataset$metadata$release_site, ' / animal reference number : ', dataset$metadata$unique_reference_code, sep = ''))

segments(
  x0 = dates[1:( length( dates) - 1)], y0 = -pres[1:( length( pres)-1)], x1 = dates[2: length( dates)], y1 = -pres[2:length( pres)], col = colorRampPalette( magma(5))( length( tempsq))[tsq]
)
axis.POSIXct( 1, at = seq( date[1], date[length(date)], by = 'months'), seq( date[1], date[length(date)], by = 'months'), format = '%b %Y')
mtext( paste( templab, ' ', '(', tempunit, ')', sep = ''), 4, line = 2.5)
vertical.image.legend( zlim = range( tempsq), colorRampPalette( magma(5))( length( tempsq)))

screen( 2)
par( mar = c( 4.5, 4.5, .5, 4.5), cex = 1.2)

plot( lon, lat, col = colorRampPalette( magma(5))( length(salsq))[ssq], pch = 19, type = 'p', cex = .5, ylim= c(-80,-40),
      xlab = lonlab, ylab = latlab)
segments(x0 = lon[1:( length( lon)-1)], y0 = lat[1:( length( lat)-1)], x1 = lon[2:length( lon)], y1 = lat[2:length(lat)], col = colorRampPalette( magma(5))( length( salsq))[ssq], lwd = 2)


lines( map$x, map$y)
mtext( paste( psallab, ' ', '(', psalunit, ')', sep = ''), 4, line = 2.75)
vertical.image.legend (zlim = range( salsq), colorRampPalette(  magma(5))( length( salsq)))
close.screen( all = TRUE)

plot( dates, -pres, col = colorRampPalette( viridis(5))( length( tempsq))[tsq], type='n',xaxt = 'n',
      xlab = datelab, ylab = paste( preslab, ' (negative dbar)', sep = ''), main = paste( dataset$metadata$species_name, ' - released in ', dataset$metadata$release_site, ' / animal reference number : ', dataset$metadata$unique_reference_code, sep = ''))


# --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Fig 2 (b) White stork                                 #
# --- --- --- --- --- --- --- --- --- --- --- --- --- ---

# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# Wind estimation by birds: 
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# This code was adapted from Weinzerl et al. 2016 'Wind estimation based on thermal soaring of birds' Ecology & Evolution

require(moveWindSpeed)

set.seed(3425)
data <- storks[[1]]

# --- --- --- --- --- --- ---# --- --- --- --- --- --- ---# --- --- --- --- --- --- ---# --- --- --- --- --- --- ---
#### Top panel:
# --- --- --- --- --- --- ---# --- --- --- --- --- --- ---# --- --- --- --- --- --- ---# --- --- --- --- --- --- ---

# Wind estimation results can only be trusted if the bird is actually flying and if it is making turns, usually while circling in a thermal. We will only consider track segments in which the bird turns at least 360 degrees and has an average speed of at least 4 m/s:

windowSize <- 29
# In order to determine an appropriate sampling rate we check the sampling rate in the stork track:
tdiffs <- diff(as.numeric(timestamps(data)))
quantile(tdiffs, c(0.01, 0.1, 0.5, 0.9, 0.99))
# The majority of data points were sampled in 1 second intervals, so we choose 1 second as our sampling interval:
isSamplingRegular <- 1

isThermallingFunction <- getDefaultIsThermallingFunction(360, 4)
# The wind estimation is based on the assumption that, within a given track segment, the air speed of the bird fluctuates randomly around a fixed mean. For high sampling rates, however, the deviations from the mean are not statistically independent, but tend to be correlated between subsequent points. We therefore need to estimate the temporal autocorrelation parameter phi (see AR(1) processes). Phi is not estimated for each track segment separately, but across all segments of a track. For performance reasons the number of segments used in the estimation of phi can be restricted:

maxPointsToUseInEstimate <- 25

# The move package has not updated to account for raster now requiring terra. I downloaded the move package source from CRAN (move_4.0.6.tar.gz) and ran their code for the distance function after loading the move package in to R. The code for their distance method is:
setMethod("distance", 
          signature=c(x=".MoveTrackSingle",y="missing"),
          definition=function(x){ 
            Dists<-raster::pointDistance(x[-sum(n.locs(x)),], x[-1,], longlat=isLonLat(x))
            return(Dists)
          })
setMethod("distance", 
          signature=c(".MoveTrackStack",y="missing"), 
          definition=function(x){
            lst <- lapply(split(x), distance)
            return(lst)
          })

# With the parameters defined so far we can run the function for estimating temporal autocorrelation parameter phi:
estimationResult <- estimatePhi(data, isSamplingRegular=isSamplingRegular, windowSize=windowSize, isThermallingFunction=isThermallingFunction, maxPointsToUseInEstimate=maxPointsToUseInEstimate, phiInitialEstimate=0, returnPointsUsedInEstimate=T)
# The estimation function returns the estimated value of phi and the number of points used in the estimate.

# We are now equipped to run the actual wind estimation along the original track. Technically wind can be estimated for every point in the track, but then the track segments used in the estimation will overlap for adjacent points. In order to avoid that we restrict the estimation to every r windowSizeth point:

isFocalPoint<-function(i, ts) i%%windowSize==0
# Now we run the estimation function:

phi <- estimationResult$phi
windEst <- getWindEstimates(data, isSamplingRegular=isSamplingRegular, windowSize=windowSize, isThermallingFunction=isThermallingFunction, phi=phi, isFocalPoint=isFocalPoint)
names(windEst)
# We only want to retain points for which the estimation was successful and where the bird was thermalling:

windEst2 <- windEst[!is.na(windEst$estimationSuccessful),]
# Here is a simple plot of the wind speeds obtained:

windSpeed <-sqrt(windEst2$windX^2+windEst2$windY^2)

#  We focus on 2 minutes where we know that high resolution trajectories were recorded while the bird was in a thermal and select the first individual. We transform the trajectory to a local projection for better plotting.

storksSub <-
  storks[strftime(timestamps(storks), format = "%H%M", tz="UTC") %in% c("1303", "1304"), ]

storksWind <- getWindEstimates(storksSub)
storksWind <- spTransform(storksWind, center = T)
indiv1 <- storksWind[[1]]

# Lets first assure that the individual is sampled with one hertz. This is crucial because some of the calculations we do below for plotting assume this.

unique(as.numeric(diff(timestamps(indiv1)), units='secs'))
# Now lets have a look at the track in one individual within the thermal. We add an arrow to each 10th location to indicate the estimated wind speed.

# We can also select one rotation and plot how head wind speed vector combined with the airspeed vector results in the ground speed.

indivSub <- indiv1[45:72,]
plot(indivSub, pch = 19)
arrows(
  coordinates(indivSub)[-n.locs(indivSub), 1],
  coordinates(indivSub)[-n.locs(indivSub), 2],
  coordinates(indivSub)[-1, 1],
  coordinates(indivSub)[-1, 2],
  length = .1,
  lwd = 1.5
)
arrows(
  coordinates(indivSub)[-n.locs(indivSub), 1],
  coordinates(indivSub)[-n.locs(indivSub), 2],
  coordinates(indivSub)[-n.locs(indivSub), 1] + head(indivSub$windX,-1),
  coordinates(indivSub)[-n.locs(indivSub), 2] + head(indivSub$windY,-1),
  col = "darkgoldenrod3",
  length = .1,
  lwd = 1.5
)
arrows(
  coordinates(indivSub)[-n.locs(indivSub), 1] + head(indivSub$windX,-1),
  coordinates(indivSub)[-n.locs(indivSub), 2] + head(indivSub$windY,-1),
  coordinates(indivSub)[-1, 1],
  coordinates(indivSub)[-1, 2],
  col = "deepskyblue",
  length = .1,
  lwd = 1.5
)
legend(
  "topright",
  legend = c("Ground speed", "Wind Speed", "Air Speed"),
  bty = "n",
  col = c("black", "darkgoldenrod3", "deepskyblue"),
  lty = 1,
  lwd = 3
)

# --- --- --- --- --- --- --- --- --- ---
#   Multiple storks in a thermal        #
# --- --- --- --- --- --- --- --- --- ---


wind_indivSub <- getWindEstimates(indivSub, windowSize = 27)# 31)

storksSub <- getWindEstimates(storksSub, windowSize = 31)

# We can also create a vertical wind speed profile through a thermal. To do so we select a thermal where most birds are present. We see that there seems to be quite a consistent pattern.

windData <- getWindEstimates(storks)
windData <- windData[strftime(timestamps(windData), format = "%H", tz="UTC") == "12" &
                       !is.na(windData$windX), ]


# windData$windX_sqrt <- windData$windX^ 2
windData$vertical_wind_profile <- sqrt(windData$windX ^ 2 + windData$windY ^ 2)
windData_df = as.data.frame(windData)
windData_df$trackId <- as.factor(windData_df$trackId)
str(windData_df$trackId)

ggplot(windData_df, aes(x = vertical_wind_profile, y = height_above_ellipsoid))+ geom_point(aes(colour = trackId)) + scale_colour_viridis_d() + theme_classic()+ggtitle('Vertical wind profile')+xlab('Wind speed')+ylab('Altitude')


ggplot(windData_df[1:1020,], aes(x = vertical_wind_profile, y = height_above_ellipsoid))+ geom_point(aes(colour = trackId)) + scale_colour_viridis_d() + theme_classic()+ggtitle('Vertical wind profile')+xlab('Wind speed')+ylab('Altitude')+ theme(legend.position="none")
