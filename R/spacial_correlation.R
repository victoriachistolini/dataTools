
#sort, order dataframe 
create_dataset <- function(day, window, params){
  
  pres_obs = get_precsence_points_obs(day,
                                      window,
                                      as.POSIXct(c("2006-06-01 00:00:00", "2013-12-31 00:00:00"), tz = 'UTC'),
                                      params)
  
  pres_obs$layer<- sapply(pres_obs$layer, function(date) paste0(substring(date,2,5), as.character(day), sep="", collapse=NULL))
  
  pres_obs = pres_obs[order(pres_obs$layer),] 
  pres_obs = na.omit(pres_obs)
  pres_obs = aggregate(pres_obs[, 4:16], list(pres_obs$layer), mean)
  
  return(pres_obs)
  
  
  
}



scatter <- function(out_dir,params, pres_obs){
  
  num_year <- as.numeric(pres_obs$Group.1)
  
  for (i in 2:length(params)+1){
    cat("  evaluating predictor:  ")
    cat(i)
    cat(" \n")
    
    ext = paste0(params[i],".jpg" , sep="",collapse = NULL)
    fpath = file.path(out_dir, ext)
    jpeg(fpath)
    
    # Generate covariate plot
    plot(num_year,pres_obs[, i] )
    dev.off
    
  }
  
}



# create another function for pairwise covariance between covariates
corr_matrx <- function(pres_obs, num_param, outLoc){
  data_matrx = matrix(nrow=14,ncol=14)
  total = num_param+1
  
  for (i in 2:num_param){
    cat(" i *** : ")
    cat(i)
    cat(" \n")
    start = i+1
    for (j in start:total){
      cat(j)
      cat(" \n")
      data_matrx[i,j] = cor(pres_obs[, i], pres_obs[, j]) 
      
    }
    
  }
  write.table(data_matrx, file = outLoc, row.names=FALSE, col.names=FALSE, sep = ",")
}


create_multi_week_set <-function(days, window, params){
  
  data=NULL
  for (i in 1:length(days)){
    partial_data = create_dataset(days[i], window, params)
    
    if(i>1){
      data = rbind(data,partial_data)
    }
    else{
      data = partial_data
    }
    
  }
  
  pres_data = data[order(data$Group.1),] 
  return(pres_data)
}



process_corr_matrix <- function(data_matrx, params){
  
  data_matrx = data_matrx[,-c(1:2)]
  chr_vector = character(13*13)
  num_vector = numeric(13*13)
  countr=0
  
  rownames(data_matrx) <- params
  colnames(data_matrx)<- params[2:length(params)]
  
  for (i in 1:length(params)){
    
    for (j in 1:(length(params)-1)){
      
      vals = data_matrx[i,j]
      #vals = data_matrx[1,1]
      if (is.na(vals)){
        next
      }
      
      if(vals > 0.4){
        chr_vector[countr] = paste(rownames(data_matrx)[i],colnames(data_matrx)[j],sep=" ~ ")
        num_vector[countr] = vals
        countr=countr+1
        
      }
      
    }
  }
  
  
  output = data.frame(num_vector,chr_vector)
  colnames(output) = c("Correlation_Coef", "Params")
  
  
  return(output)
  
  
  
}