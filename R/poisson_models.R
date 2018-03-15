# helper methods for formatting a point pattern analysis

# find most commonly occuring year in tick dataset
find_most_common_year <-function(obs,day){

  day = as.character(day)
  years <- sapply(obs$layer, function(date) substring(date,2,5))
  factor(years)

  year_counts = sort(table(years), decreasing=TRUE)
  year_string = names(year_counts[1])
  common_year_string = paste0("A", year_string, day , sep="",collapse = NULL)

  return(common_year_string)
}


# convert raster to image format
# @param raster_load raster stack
# @param year to load
# @return an image representation of a raster
format.im <-function(raster_load,year){

  year_list = names(raster_load)
  layer_to_load = match(year,year_list)

  read_img = as.im(raster(raster_load, layer=layer_to_load))
  return(read_img)

}

# returns sorted predictors as images
# creates a list of images to return
create_point_process_predictors <- function(day,year, predictor_stack){

  year <-  paste0("A", year, day , sep="",collapse = NULL)
  predictors.im = lapply(predictor_stack,function(d) format.im(d,year))

  return(predictors.im)

}

# converts tick dataframe to a ppp format, only using observations from the given year
ticks_to_ppp <- function(obs, year, proj=PROJ ){

  obs$layer = sapply(obs$layer, function(date) substring(date,2,5))
  obs = filter(obs,layer == year)

  cat("  number of points in observations set:  ")
  cat(nrow(obs))
  cat(" \n")

  obs = obs[,c("x","y")]
  ticks.sp = SpatialPoints(obs,CRS(proj))
  ticks.ppp  <- as(ticks.sp, "ppp")

  return(ticks.ppp)
}

# return dataset of tick obs and extracted environmental covariates
set_up_point_process <- function(SS,day, window,params,year){
  # load/format tick observations
  obs = trim_prescence_points(tickdata::read_obs(),
                              SS,
                              day, window)

  predictor_stack <- process_predictor_stack(SS, params, window, day)
  prescence_points <- extract_precsence_points(predictor_stack,obs)

  dataset <- cbind(obs,prescence_points)

  dataset<-filter(dataset, grepl(year, dataset$layer))

  return(dataset)

}

# create multi histogram of data
explore_covariates <-function(dataset, params,year){

  cpredictrs = dataset[,params]

  # do some plot stuff
  fpath = paste("/home/vchisto/epx2/150/",year, sep="")
  jpeg(fpath)
  d <- melt(cpredictrs)
   ggplot(d,aes(x = value)) +
         facet_wrap(~variable,scales = "free_x") +
         geom_histogram()
  dev.off()



}


convert_obs <-function(dataset){

  # grab lat long
  current_obs = dataset[,c("x","y")]
  # reformat into spatstat formats
  ticks.sp = SpatialPoints(current_obs,CRS(PROJ))
  ticks.ppp  <- as(ticks.sp, "ppp")
}

######## CONSTANTS ########
PROJ <- c(longlat =  "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0",
          lcc = "+proj=lcc +lat_1=25 +lat_0=25 +lon_0=-95 +k_0=1 +x_0=0 +y_0=0 +a=6367470.21484375 +b=6367470.21484375 +units=km +no_defs")
