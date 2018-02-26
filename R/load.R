# loading 
#
# help to load data / extract points 
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

library(dplyr)
library(tibble)
library(sp)
library(namtools)
library(raster)
library(dismo)
library(dismotools)
library(tickdata)
library(rscripting)
library(binaryLogic)
library(spatstat)  
library(rgdal)     
library(maptools)  
library(raster)


#' Compute the doy window
#' @param x the day of year (should be 1-366)
#' @param w the window as 2 elements [days before, days after]
#' @param MAX the highest possible day number
#' @return numeric doy vector
doy_window <- function(x = 1, w = c(-5,5), MAX = 366){
  newday <- x + seq(from = w[1], to = w[2])
  ix <- newday < 1
  if (any(ix)) newday[ix] = MAX + newday[ix]
  ix <- newday > MAX
  if (any(ix)) newday[ix] = newday[ix] - MAX
  newday
}

#' Compute the union of a list of vectors
#' 
#' @param x a list of zero or more vectors
#' @return a vector of the union or NULL is the input is empty
munion <- function(x){
  if (!is.list(x)) stop("input must be a list")
  
  n <- length(x)
  if (n == 0) return(NULL)
  if (n == 1) return(x[[1]])
  
  s <- x[[1]]
  for (i in 2:n) s <- union(s, x[[i]])
  s
}

#' Compute the intersection of a list of vectors
#' 
#' @param x a list of zero or more vectors
#' @return a vector of the intersection or NULL is the input is empty
mintersect <- function(x){
  n <- length(x)
  if (n == 0) return(NULL)
  if (n == 1) return(x[[1]])
  
  s <- x[[1]]
  for (i in 2:n) s <- intersect(s, x[[i]])
  s
}

#' Project tick obseravtion coordinates
#' 
#' @param x a tibble with 'lon' and 'lat' columns
#' @param from_proj character source projection
#' @param to_proj character destiantion projection
#' @return the original input with updated x and y columns
project_ticks <- function(x,
                          from_proj = PROJ[['longlat']], 
                          to_proj = PROJ[['lcc']]){
  
  ll <- as.data.frame(dplyr::select(x, lon, lat))
  sp::coordinates(ll) <- ~lon+lat
  proj4string(ll) <- from_proj
  ll <- sp::coordinates(sp::spTransform(ll, to_proj))
  x$x <- ll[,1]
  x$y <- ll[,2]
  x
}

#trim tick obs
trim_prescence_points <- function(
  obs,
  ss,
  doy,
  daywindow){
  
  
  
  # trim the observations to the user specified date range
  TICKS <- obs %>%
    dplyr::slice(which(findInterval(date,ss) == 1))
  # project to match the NAM predictiors
  TICKS <- project_ticks(TICKS)
  # add a layer id so we can navigate the predictor stacks/bricks/layers
  TICKS$layer <- format(TICKS$date, "A%Y%j")
  
  # the the collection of days around 'this' day
  days <- doy_window(doy, w = daywindow)
  
  # trim the obs to our dayswindow
  cat("  selecting presence points\n")
  obs <- TICKS %>% 
    dplyr::filter(doy %in% days) %>%
    dplyr::select(x, y, layer) %>%
    as.data.frame()
  
  # total number of points 
  return(obs)
  
}

process_predictor_stack <- function(ss, param, window,doy){
  cat("  selecting predictor stack\n")
  #crop region  #NOt correct region 
  DAY <- namtools::read_db_nwaday() %>%
    dplyr::slice(which(dplyr::between(D,ss[1],ss[2]) ))
  
  DAY$doy <- as.numeric(format(DAY$D, "%j"))
  
  pp <- split(DAY, paste0(DAY$trt,'_', DAY$n))
  pp <- pp[param]
  
  PP <- namtools::read_nam218_daybrick(params = names(pp))
  
  JD <- sapply(PP,
               function(P){
                 format(getZ(P), "%j")
               })
  jd <- mintersect(JD)
  
  #get window
  days <- doy_window(doy, w = window)
  jdays <- sprintf("%0.3i", days)
  
  
  # trimming predictors stack
  jd <- sapply(names(PP), function(n) which(JD[[n]] %in% jdays), simplify = FALSE )
  pp <- sapply(names(PP), function(n) {
    p <- PP[[n]][[jd[[n]]]]
    raster::setZ(p, as.POSIXct(paste(names(p), "00:00:00"), format = 'A%Y%j %H:%M:%S', tz = 'UTC'))
  })
  
  
  return(pp)
  
  
}

extract_precsence_points <- function(pp,obs){
  cat("  extracting points\n")
  yy <- sapply(names(pp),
               function(pname, x = NULL){
                 dismotools::layers_extractPoints(pp[[pname]], x)
               },
               x = obs %>% dplyr::select(x,y,layer) %>% as.data.frame(), simplify = FALSE)
  prs <- as.data.frame(yy)
  
  return(prs)
  
}


extract_back_points <- function(pp,obs,nback){
  
  bkgpts <- dismotools::layers_randomPoints(pp[[1]], N = nback,
                                            pts = obs %>% dplyr::select(x,y,layer) %>% as.data.frame())
  bkgpts <- tibble::as.tibble(bkgpts)
  
  bkg <- sapply(names(pp),
                function(pname, x = NULL){
                  dismotools::layers_extractPoints(pp[[pname]], x)
                },
                x = bkgpts %>% dplyr::select(x = lon, y = lat, layer) %>% as.data.frame(),
                simplify = FALSE)
  bkg <- as.data.frame(bkg)
  
  return(bkg)
  
}

# return auc / num pres points used 
create_model <- function(day, window,ss,params){
  
  obs = trim_prescence_points(tickdata::read_obs(), ss, day, window)
  
  # download prescence/background data
  predictor_stack <- process_predictor_stack(ss, params, window, day)
  prescence_points <- extract_precsence_points(predictor_stack,obs)
  background_points <- extract_back_points(predictor_stack,obs, NBACK)
  
  
  # run MAXENT model
  flag <- c(rep(1, nrow(prescence_points)), rep(0, nrow(background_points)))
  input <- rbind(prescence_points, background_points)
  mpath <- "/dev/shm/model"
  model <- try(dismo::maxent(input, flag, path = mpath))
  
  # store model
  if (!inherits(model, 'try-error') ) {
    # extract the AUC
    auc <- dismotools::maxent_get_results(model, 'auc')
    
  }
  
  #delete model
  if (dir.exists(mpath)) unlink(mpath, force = TRUE, recursive = TRUE)
  
  out <- c(nrow(prescence_points), auc)
  return(out)
}


get_precsence_points_obs <- function(day, window,ss,params) {
  
  obs = trim_prescence_points(tickdata::read_obs(), ss, day, window)
  predictor_stack <- process_predictor_stack(ss, params, window, day)
  prescence_points <- extract_precsence_points(predictor_stack,obs)
  
  pres_points_obs <- cbind(obs, prescence_points)
  
  return(pres_points_obs)
  
}


#### CONSTANTS ####
PROJ <- c(longlat =  "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0",
          lcc = "+proj=lcc +lat_1=25 +lat_0=25 +lon_0=-95 +k_0=1 +x_0=0 +y_0=0 +a=6367470.21484375 +b=6367470.21484375 +units=km +no_defs")



