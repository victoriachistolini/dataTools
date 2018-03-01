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
format.im <-function(raster_load,year){
  
  year_list = names(raster_load)
  layer_to_load = match(year,year_list)
  
  read_img = as.im(raster(raster_load, layer=layer_to_load))
  return(read_img)
  
}

# returns sorted predictors as images
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



######## CONSTANTS ########
PROJ <- c(longlat =  "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0",
          lcc = "+proj=lcc +lat_1=25 +lat_0=25 +lon_0=-95 +k_0=1 +x_0=0 +y_0=0 +a=6367470.21484375 +b=6367470.21484375 +units=km +no_defs")