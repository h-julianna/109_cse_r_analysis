Bf<- function(sd, obtained, dfdata = 1, likelihood = c("normal", "t"), 
              modeloftheory = c("normal", "t", "cauchy", "uniform"), uniform, 
              lower=0, upper=1, meanoftheory=0, sdtheory=1, modeoftheory =0, scaleoftheory = 1, 
              dftheory = 1, tail=1, method = "old"){
  if (method == "old") {
    area <- 0
    if(identical(uniform, 1)){
      theta <- lower
      range <- upper - lower
      incr <- range / 2000
      for (A in -1000:1000){
        theta <- theta + incr
        dist_theta <- 1 / range
        height <- dist_theta * dnorm(obtained, theta, sd)
        area <- area + height * incr
      }
    }else
    {theta <- meanoftheory - 5 * sdtheory
    incr <- sdtheory / 200
    for (A in -1000:1000){
      theta <- theta + incr
      dist_theta <- dnorm(theta, meanoftheory, sdtheory)
      if(identical(tail, 1)){
        if (theta <= 0){
          dist_theta <- 0
        } else {
          dist_theta <- dist_theta * 2
        }
      }
      height <- dist_theta * dnorm(obtained, theta, sd)
      area <- area + height * incr
    }
    }
    LikelihoodTheory <- area
    Likelihoodnull <- dnorm(obtained, 0, sd)
    BayesFactor <- LikelihoodTheory / Likelihoodnull
    ret <- list("LikelihoodTheory" = LikelihoodTheory,"Likelihoodnull" = Likelihoodnull, "BayesFactor" = BayesFactor)
    ret
  } else if (method == "new") {
    if(likelihood=="normal"){
      dfdata=10^10
    }
    if(modeloftheory=="normal"){
      dftheory = 10^10
    } else if(modeloftheory=="cauchy"){
      dftheory = 1
    }
    area <- 0
    normarea <- 0
    if(modeloftheory=="uniform"){
      theta <- lower
      range <- upper - lower
      incr <- range / 2000
      for (A in -1000:1000){
        theta <- theta + incr
        dist_theta <- 1 / range
        height <- dist_theta * dt((obtained-theta)/sd, df=dfdata)
        area <- area + height * incr
      }
      LikelihoodTheory <- area
    }else{
      theta <- modeoftheory - 10 * scaleoftheory
      incr <- scaleoftheory/200
      for (A in -2000:2000){
        theta <- theta + incr
        dist_theta <- dt((theta-modeoftheory)/scaleoftheory, df=dftheory)
        if(identical(tail, 1)){
          if (theta <= modeoftheory){
            dist_theta <- 0
          } else {
            dist_theta <- dist_theta * 2
          }
        }
        height <- dist_theta * dt((obtained-theta)/sd, df = dfdata)
        area <- area + height * incr
        normarea <- normarea + dist_theta*incr
      }
      LikelihoodTheory <- area/normarea
    }
    Likelihoodnull <- dt(obtained/sd, df = dfdata)
    BayesFactor <- LikelihoodTheory/Likelihoodnull
    BayesFactor
  }
}

Bf_range <- function(sd, obtained, dfdata = 1, likelihood = c("normal", "t"), 
                     modeloftheory= c("normal","t","cauchy") , meanoftheory = 0,
                     modeoftheory = 0, sdtheoryrange, dftheory = 1, tail = 1, method = "old")
{
  if (method == "old") {
    
    x = c(0)
    y = c(0)
    
    for(sdi in sdtheoryrange)
    {
      #sdi = sdtheoryrange[1]
      # uses old Bf method
      B = as.numeric(Bf(sd, obtained, meanoftheory=0, uniform = 0, sdtheory=sdi, tail=tail)[3])
      
      #following line corrects for the fact that the calcuator does not correctly compute BF when sdtheory==0; this code ensures that if sdtheory ==0, BF=1
      
      if (sdi ==0 ) {B=1}
      
      x= append(x,sdi)  
      y= append(y,B)
      output = cbind(x,y)
      
    }
    
  }
  else if (method == "new") {
    
    x = c(0)
    y = c(0)
    
    for(sdi in sdtheoryrange)
    {
      #sdi = sdtheoryrange[1]
      # uses new Bf method
      B = as.numeric(Bf(sd = sd, obtained = obtained, dfdata = dfdata, likelihood = likelihood, 
                        modeloftheory = modeloftheory, modeoftheory=modeoftheory, scaleoftheory=sdi, 
                        dftheory = dftheory, tail = tail, method="new"))
      
      #following line corrects for the fact that the calcuator does not correctly compute BF when sdtheory==0; this code ensures that if sdtheory ==0, BF=1
      
      if (sdi ==0 ) {B=1}
      
      x= append(x,sdi)  
      y= append(y,B)
      output = cbind(x,y)
      
    }
    
  }
  output = output[-1,] 
  colnames(output) = c("sdtheory", "BF")
  return(output) 
}

old_Bf <- function(sd, obtained, dfdata, meanoftheory, sdtheory, dftheory, tail = 2)
  
{
  
  area <- 0
  
  normarea <- 0
  
  theta <- meanoftheory - 10 * sdtheory
  
  incr <- sdtheory/200
  
  for (A in -2000:2000){
    
    theta <- theta + incr
    
    dist_theta <- dt((theta-meanoftheory)/sdtheory, df=dftheory)
    
    if(identical(tail, 1)){
      
      if (theta <= 0){
        
        dist_theta <- 0
        
      } else {
        
        dist_theta <- dist_theta * 2
        
      }
      
    }
    
    height <- dist_theta * dt((obtained-theta)/sd, df = dfdata)
    
    area <- area + height * incr
    
    normarea <- normarea + dist_theta*incr
    
  }
  
  LikelihoodTheory <- area/normarea
  
  Likelihoodnull <- dt(obtained/sd, df = dfdata)
  
  BayesFactor <- LikelihoodTheory/Likelihoodnull
  
  BayesFactor
  
}