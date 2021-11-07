GenerateData_Deng <- function(m, pi, par1, par2, distribution, time){
  
  pp <- 1 - pi
  for (run in 1:1000){
    
    #### Generating data
    set.seed(run)
    obs_time = rep(NA,m)
    for (i in 1:m){
      select = rbinom(1,1,pp)
      if (select==0){
        if(distribution == "gamma"){
          obs_time[i] = rgamma(1,par1,par2)
        }
        else if(distribution == "weibull"){
          obs_time[i] = rweibull(1,par1,par2) 
        }
        else if(distribution == "lognormal"){
          obs_time[i] = exp(rnorm(1,par1,par2)) 
        }
      }
      else{
        while (TRUE){
          C = runif(1,0,30)		#departure time
          if(distribution == "gamma"){
            Y = rgamma(1,par1,par2)
          }
          else if(distribution == "weibull"){
            Y = rweibull(1,par1,par2) 
          }
          else if(distribution == "lognormal"){
            Y = exp(rnorm(1,par1,par2))
          }
          if (Y>C) break
        }
        if(time == "forward"){ 
          obs_time[i] <- round(Y-C,4) #observed forward time
        }
        else if(time == "backward"){  
          obs_time[i] = round(Y-(Y-C),4)		#observed backward time
        }  
      }
    }
  }
  return(obs_time)
}  

GenerateData_Deng(1200, 0.2, 5, 0.8, "gamma", "forward")