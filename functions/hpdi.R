
# function to calculate the highest posterior density interval for the quantiles
# these intervals are often a better descriptor of skewed posterior distributions
interval_function_hpdi <- function(x,probs, names = FALSE){
  y <- vector("numeric",length = length(probs))
  if(names){
    names(y) <- paste0(probs*100,"%")
  }
  for(j in 1:length(probs)){
    prob <- probs[j]
    if(prob > 0.67 | prob < 0.33){
      if(prob < 0.33){
        q2 <- 1-(prob*2)
        y[j] <- HDInterval::hdi(x,q2)[1]

      }else{
        q2 <- 1-((1-prob)*2)
        y[j] <- HDInterval::hdi(x,q2)[2]

      }
    }else{
      y[j] <- stats::quantile(x,prob)
    }
  }
  return(y)
}
