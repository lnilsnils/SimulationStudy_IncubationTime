simulate_outbreak_wuhan <- function(total_cases, leaveenter_per_day_prior19th, leave_per_day_from19thonwards, growth_rate, distribution, par1, par2, time){
  
  infected_list <- list()
  incubation_periods <- list()
  infected_travel <- list()
  time_obs <- list()
  left_wuhan <- list()
  population_total <- c(1:11000000) #population Wuhan
  population_minus_infected <- c(1:11000000) # same as population_total at the start - to ensure that not > 1 infections per person
  
  day_of_inf <- vector(mode = "list", length = 19)
  
  for (i in 1:19){
    #for infected_list- draw a sample from people who are still in Wuhan and have not been infected before -> population_minus_infected
    infected_list[[i]] <- sample(x = population_minus_infected, size = total_cases, replace = FALSE)
    if (i < 15){
      #add people who enter city
      population_total <- c(population_total, (tail(population_total, 1) + 1):(tail(population_total, 1) + leaveenter_per_day_prior19th))
      #for left_wuhan - draw a sample from everyone that is in the city on that day
      left_wuhan[[i]] <- sample(x = population_total, size = leaveenter_per_day_prior19th, replace = FALSE)
      population_total <- population_total[!(population_total %in% left_wuhan[[i]])]
      population_minus_infected <- population_total[!(population_total %in% unlist(infected_list[1:i]))] 
    }
    else{
      #draw people who travel out Wuhan from population
      left_wuhan[[i]] <- sample(x = population_total, size = leave_per_day_from19thonwards, replace = FALSE)
      population_total <- population_total[!(population_total %in% left_wuhan[[i]])]
      population_minus_infected <- population_total[!(population_total %in% unlist(infected_list[1:i]))] 
      infected_travel[[i]] <- intersect(left_wuhan[[i]], unlist(infected_list))
      
      #look up time of infection for each infected traveler in List of infected
      for (k in 1:length(infected_travel[[i]])){
        day_of_inf[[i]][k] <- which(sapply(X = infected_list, FUN = function(X) infected_travel[[i]][k] %in% X))[1]
      }
      #for those the compute incubation periods are
      if(distribution == "gamma"){
        incubation_periods[[i]] <- rgamma(n = length(infected_travel[[i]]), shape = par1, rate = par2)
      } else if(distribution == "weibull"){
        incubation_periods[[i]] <- rweibull(n = length(infected_travel[[i]]), shape = par1, scale = par2)
      } else if(distribution == "lognormal"){
        incubation_periods[[i]] <- rlnorm(n = length(infected_travel[[i]]), meanlog = par1, sdlog = par2)
      }
      
      if(time == "forward"){ #forward time
        #the forward times when symptom onset is not before day of travel are:
        time_obs[[i]] <- (day_of_inf[[i]][day_of_inf[[i]] + incubation_periods[[i]] > i] +
                            incubation_periods[[i]][day_of_inf[[i]] + incubation_periods[[i]] > i]) - i
      }
      else if(time == "backward"){ #backward time
        time_obs[[i]] <- i - (day_of_inf[[i]][day_of_inf[[i]] + incubation_periods[[i]] > i])
      }
    }
    #number of cases per day
    total_cases <- total_cases + round(total_cases * (exp(growth_rate) - 1))
  }
  return(unlist(time_obs))
}


