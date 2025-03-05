## function to estimate full posterior of abundance i.e., regional population sizes

post_abund <- function(x, #data frame with summed relative abundance within a given region.
                       draws = cali_use$calibration, # posterior draws of caslibration
                       fun = "mean", # summary statistic to return alternates are "median", "quantile" or "hpdi"
                       p = 0.025,
                       use_log = FALSE){

  y <- vector("numeric",length(x))

  if(use_log){
    xl <- log(x)
  }
  if(fun == "mean"){
    for(i in 1:length(x)){
      if(use_log){
        y[i] <- exp(mean(xl[i]+as.numeric(draws)))
      }else{
        y[i] <- mean(x[i]*as.numeric(draws))
      }
    }
  }

  if(fun == "median"){
    for(i in 1:length(x)){
      if(use_log){
        y[i] <- exp(median(xl[i]+as.numeric(draws)))
      }else{
        y[i] <- median(x[i]*as.numeric(draws))
      }
    }
  }

  if(fun == "hpdi"){
    for(i in 1:length(x)){
      if(use_log){
        y[i] <- exp(interval_function_hpdi(xl[i]+as.numeric(draws),p, names = FALSE))

      }else{
        y[i] <- interval_function_hpdi(x[i]*as.numeric(draws),p, names = FALSE)
      }
    }
  }


  if(fun == "quantile"){
    for(i in 1:length(x)){
      if(use_log){
        y[i] <- exp(quantile(xl[i]+as.numeric(draws),p, names = FALSE))

      }else{
        y[i] <- quantile(x[i]*as.numeric(draws),p, names = FALSE)
      }
    }
  }

  return(y)
}

