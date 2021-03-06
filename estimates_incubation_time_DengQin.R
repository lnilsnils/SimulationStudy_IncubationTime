incubationPeriod_Estimates_Deng <- function(runs, distribution, par1, par2){
  
  result = NULL
  
  aa <- par1
  bb <- par2
  
  
  if(distribution == "gamma"){
    ##########################################gamma incubation
    f <- function(y,a,b){
      return(dgamma(y,shape=a,rate=b))
    }
    
    #pdf of forward time
    h <- function(y,a,b){
      RE = b/a*(1-pgamma(y,shape=a,rate=b))
      return(RE)
    }
    FI <- function(y,a,b) return(pgamma(y,shape=a,rate=b))
    
    F <- function(k,a,b,p){
      RE = rep(0,length(k))
      posit = which(k>0)
      REh = pgamma(k[posit],shape=a+1,rate=b)+k[posit]*b/a*(1-pgamma(k[posit],shape=a,rate=b))
      if (p>0.999) RE[posit] = REh
      else{
        REf = pgamma(k[posit],shape=a,rate=b)
        RE[posit] = REf*(1-p)+REh*p
      }
      return(RE)
    }
    
    
    
  }else if(distribution == "weibull"){
    ##########################################weibull incubation
    f <- function(y,a,b){
      return(dweibull(y,shape=a,scale=b))
    }
    h <- function(y,a,b){
      RE = (1-pweibull(y,shape=a,scale=b))/gamma(1+1/a)/b
      return(RE)
    }
    FI <- function(y,a,b) return(pweibull(y,shape=a,scale=b))
    F <- function(k,a,b,p){
      RE = rep(0,length(k))
      posit = which(k>0)
      REh = pgamma((k[posit]/b)^a,shape=1/a,rate=1)
      if (p>0.999) RE[posit] = REh
      else{
        REf = pweibull(k[posit],shape=a,scale=b)
        RE[posit] = REf*(1-p)+REh*p
      }
      return(RE)
    }
    
  }else if(distribution == "lognormal"){
    ##########################################lognormal incubation
    f <- function(y,u,s){
      RE = rep(0,length(y))
      posit = which(y>0)
      RE[posit]=dnorm(log(y[posit]),u,s)/y[posit]
      return(RE)
    }
    h <- function(y,u,s){
      RE = exp(-u-s^2/2)*(1-pnorm(log(y),u,s))
      return(RE)
    }
    FI <- function(y,u,s){
      RE = rep(0,length(y))
      posit = which(y>0)
      RE[posit]=pnorm(log(y[posit]),u,s)
      return(RE)
    }
    F <- function(k,u,s,p){
      RE = rep(0,length(k))
      posit = which(k>0)
      REf = pnorm(log(k[posit]),u,s)
      REh = pnorm(log(k[posit]),u+s^2,s)+exp(-u-s^2/2)*k[posit]*(1-pnorm(log(k[posit]),u,s))
      RE[posit] = REf*(1-p)+REh*p  # equation 2 of paper Deng
      return(RE)
    }
    
  }
  
  loglik <- function(pa){
    a = pa[1]
    b = pa[2]
    p = pa[3]
    if (p>0.999) P=h(x,a,b)
    else P = p*h(x,a,b)+(1-p)*f(x,a,b)
    P[P<0.00001]=0.00001
    RE = - sum( log(P) )
    if (RE>10000) RE=10000
    return(RE)
  }
  
  
  Loglik <- function(pa){
    a = pa[1]
    b = pa[2]
    p = pa[3]
    P = F(x+0.5,a,b,p)-F(x-0.5,a,b,p)
    P[P<0.00001]=0.00001
    RE = - sum( log(P) )
    if (RE>10000) RE=10000
    return(RE)
  }
  
  
  
  for (w in 1:9){
    if (w%%3==1) m=600
    if (w%%3==2) m=1200
    if (w%%3==0) m=1800
    if (w<=3) pp=1
    if (w>=4 & w <=6) pp=0.9
    if (w>=7) pp=0.8
    
    
    a.vec = numeric(runs)
    b.vec = numeric(runs)
    q.vec = numeric(runs)
    ave.vec = numeric(runs)
    Q1.vec = numeric(runs)
    Q2.vec = numeric(runs)
    Q3.vec = numeric(runs)
    Q4.vec = numeric(runs)
    Q5.vec = numeric(runs)
    Q6.vec = numeric(runs)
    
    a.vec.qin = numeric(runs)
    b.vec.qin = numeric(runs)
    q.vec.qin = numeric(runs)
    ave.vec.qin = numeric(runs)
    Q1.vec.qin = numeric(runs)
    Q2.vec.qin = numeric(runs)
    Q3.vec.qin = numeric(runs)
    Q4.vec.qin = numeric(runs)
    Q5.vec.qin = numeric(runs)
    Q6.vec.qin = numeric(runs)
    
    bias_par1_Deng = numeric(runs)
    bias_par2_Deng = numeric(runs)
    bias_p_Deng = numeric(runs)
    
    bias_par1_Qin = numeric(runs)
    bias_par2_Qin = numeric(runs)
    bias_p_Qin = numeric(runs)
    
    
    for (run in 1:runs){
      #### Generating data
      set.seed(run)
      x = rep(NA,m)
      for (i in 1:m){
        select = rbinom(1,1,pp) #additional infected at day of travel
        if (select==0){
          if(distribution == "gamma"){
            x[i] = rgamma(1,aa,bb)
          }else if(distribution == "weibull"){
            x[i] = rweibull(1,aa,bb)
          }else if(distribution == "lognormal"){
            x[i] = rlnorm(1,aa,bb)
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
          x[i] <- (Y-C) #observed forward time
        }
      }
      #(suggestion for correction: select <- size 1000 0/1, x-1000 from gamma (if gamma),y-1000 simulate_outbreakWuhan, x*select+y*(1-select), ignore for now)
      
      #### p: censor probability
      
      x = round(x,0)
      x[x>25] = 25 
      
      
      
      if(distribution == "gamma"){
        par_1 <- optim(par=c(4,0.5,0.8),Loglik,method='L-BFGS-B',
                       lower=c(1.1,0.1,0),upper=c(10,2,1))$par
      }
      else if(distribution == "weibull"){
        par_1 <- optim(par=c(2,10,0.8),Loglik,method='L-BFGS-B',
                       lower=c(1.1,2,0),upper=c(5,15,1))$par
        
      }
      else if(distribution == "lognormal"){
        par_1 <- optim(par=c(2,0.4,0.8),Loglik,method='L-BFGS-B',
                       lower=c(1,0.1,0),upper=c(5,1,1))$par
      }
      
      #model parameter and pi
      
      a = par_1[1]
      b = par_1[2]
      p = par_1[3]
      
      a.vec[run] = a
      b.vec[run] = b
      q.vec[run] = 1-p
      
      
      
      
      #mean, median, quartiles
      if(distribution == "gamma"){
        #########################gamma
        est <- c(
          c(a,b,1-p) ,  # estimates for parameters alpha and beta, and pi
          a/b,         # mean
          round(qgamma(c(0.25,0.5,0.75,0.9,0.95,0.99),a,b),2))  # estimated quantiles
        # -round(loglik(c(a,b)),2),
        
      } else if(distribution == "weibull"){
        #########################weibull
        est <- c(
          c(a,b,1-p) ,   # estimates for parameters k and lambda, and pi
          b*gamma(1+1/a) , # mean
          round(qweibull(c(0.25,0.5,0.75,0.9,0.95,0.99),shape=a,scale=b),2))  # estimated quantiles
        #-round(loglik(c(a,b)),2),
        
      } else if(distribution == "lognormal"){
        ########################lognormal
        est <- c(
          c(a,b,1-p),  # estimates for parameters mu and sigma, and pi
          exp(a+b^2/2), # mean
          round(exp(qnorm(c(0.25,0.5,0.75,0.9,0.95,0.99),a,b)),2))  # estimated quantiles
        #-round(loglik(c(a,b)),2)
        
      }
      
      #mean, median, quartiles, bias
      
      ave.vec[run] = est[4]
      Q1.vec[run] = est[5]
      Q2.vec[run] = est[6]
      Q3.vec[run] = est[7]
      Q4.vec[run] = est[8]
      Q5.vec[run] = est[9]
      Q6.vec[run] = est[10]
      
      
      bias_par1_Deng[run] <- est[1] - aa
      bias_par2_Deng[run] <- est[2] - bb
      bias_p_Deng[run] <- est[3] - (1 - pp)
      
      
      
      if(distribution == "gamma"){
        par_2 <- optim(par=c(4,0.5,0.8),loglik,method='L-BFGS-B',
                       lower=c(1.1,0.1,0),upper=c(10,2,1))$par
      }
      else if(distribution == "weibull"){
        par_2 <- optim(par=c(2,10,0.8),loglik,method='L-BFGS-B',
                       lower=c(1.1,2,0),upper=c(5,15,1))$par
      }
      else if(distribution == "lognormal"){
        par_2 <- optim(par=c(2,0.4,0.8),loglik,method='L-BFGS-B',
                       lower=c(1,0.1,0),upper=c(5,1,1))$par
      }
      
      #model parameter and pi
      
      a = par_2[1]
      b = par_2[2]
      p = par_2[3]
      
      a.vec.qin[run] = a
      b.vec.qin[run] = b
      q.vec.qin[run] = 1-p
      
      
      #mean, median, quartiles
      if(distribution == "gamma"){
        #########################gamma
        est <- c(
          c(a,b,1-p),  # estimates for parameters alpha and beta, and pi
          a/b,         # mean
          round(qgamma(c(0.25,0.5,0.75,0.9,0.95,0.99),a,b),2))  # estimated quantiles
        # -round(loglik(c(a,b)),2),
        
      } else if(distribution == "weibull"){
        #########################weibull
        est <- c(
          c(a,b,1-p),   # estimates for parameters k and lambda, and pi
          b*gamma(1+1/a), # mean
          round(qweibull(c(0.25,0.5,0.75,0.9,0.95,0.99),shape=a,scale=b),2))  # estimated quantiles
        #-round(loglik(c(a,b)),2),
        
      } else if(distribution == "lognormal"){
        ########################lognormal
        est <- c(
          c(a,b,1-p),  # estimates for parameters mu and sigma, and pi
          exp(a+b^2/2), # mean
          round(exp(qnorm(c(0.25,0.5,0.75,0.9,0.95,0.99),a,b)),2))  # estimated quantiles
        #-round(loglik(c(a,b)),2)
        
      }
      
      
      #mean, median, quartiles, bias
      
      ave.vec.qin[run] = est[4]
      Q1.vec.qin[run] = est[5]
      Q2.vec.qin[run] = est[6]
      Q3.vec.qin[run] = est[7]
      Q4.vec.qin[run] = est[8]
      Q5.vec.qin[run] = est[9]
      Q6.vec.qin[run] = est[10]
      
      
      bias_par1_Qin[run] <- est[1] - aa
      bias_par2_Qin[run] <- est[2] - bb
      bias_p_Qin[run] <- est[3] - (1 - pp)
      
      
    }
    
    #MSE and coverage
    MSE_par1_Deng <- sum(bias_par1_Deng^2)/runs
    MSE_par2_Deng <- sum(bias_par2_Deng^2)/runs
    MSE_p_Deng <- sum(bias_p_Deng^2)/runs
    
    MSE_par1_Qin <- sum(bias_par1_Qin^2)/runs
    MSE_par2_Qin <- sum(bias_par2_Qin^2)/runs
    MSE_p_Qin <- sum(bias_p_Qin^2)/runs
    
    
    
    res <- rbind(c(round(mean(a.vec),2), round(t.test(a.vec)$conf.int[1], 2), round(t.test(a.vec)$conf.int[2], 2)),
                 c(round(mean(b.vec),2), round(t.test(b.vec)$conf.int[1], 2), round(t.test(b.vec)$conf.int[2], 2)),
                 c(round(mean(q.vec),2), round(t.test(q.vec)$conf.int[1], 2), round(t.test(q.vec)$conf.int[2], 2)),  
                 c(round(mean(ave.vec),2), round(t.test(ave.vec)$conf.int[1], 2), round(t.test(ave.vec)$conf.int[2], 2)),
                 c(round(mean(Q1.vec),2), round(t.test(Q1.vec)$conf.int[1], 2), round(t.test(Q1.vec)$conf.int[2], 2)),
                 c(round(mean(Q2.vec),2), round(t.test(Q2.vec)$conf.int[1], 2), round(t.test(Q2.vec)$conf.int[2], 2)),
                 c(round(mean(Q3.vec),2), round(t.test(Q3.vec)$conf.int[1], 2), round(t.test(Q3.vec)$conf.int[2], 2)),
                 c(round(mean(Q4.vec),2), round(t.test(Q4.vec)$conf.int[1], 2), round(t.test(Q4.vec)$conf.int[2], 2)),
                 c(round(mean(Q5.vec),2), round(t.test(Q5.vec)$conf.int[1], 2), round(t.test(Q5.vec)$conf.int[2], 2)),
                 c(round(mean(Q6.vec),2), round(t.test(Q6.vec)$conf.int[1], 2), round(t.test(Q6.vec)$conf.int[2], 2)),
                 c(round(MSE_par1_Deng, 2), NULL, NULL),
                 c(round(MSE_par2_Deng, 2), NULL, NULL),
                 c(round(MSE_p_Deng, 2), NULL, NULL),
                 
                 c(round(mean(a.vec.qin),2), round(t.test(a.vec.qin)$conf.int[1], 2), round(t.test(a.vec.qin)$conf.int[2], 2)),
                 c(round(mean(b.vec.qin),2), round(t.test(b.vec.qin)$conf.int[1], 2), round(t.test(b.vec.qin)$conf.int[2], 2)),
                 c(round(mean(q.vec.qin),2), round(t.test(q.vec.qin)$conf.int[1], 2), round(t.test(q.vec.qin)$conf.int[2], 2)),
                 c(round(mean(ave.vec),2), round(t.test(ave.vec)$conf.int[1], 2), round(t.test(ave.vec)$conf.int[2], 2)),
                 c(round(mean(Q1.vec.qin),2), round(t.test(Q1.vec.qin)$conf.int[1], 2), round(t.test(Q1.vec.qin)$conf.int[2], 2)),
                 c(round(mean(Q2.vec.qin),2), round(t.test(Q2.vec.qin)$conf.int[1], 2), round(t.test(Q2.vec.qin)$conf.int[2], 2)),
                 c(round(mean(Q3.vec.qin),2), round(t.test(Q3.vec.qin)$conf.int[1], 2), round(t.test(Q3.vec.qin)$conf.int[2], 2)),
                 c(round(mean(Q4.vec.qin),2), round(t.test(Q4.vec.qin)$conf.int[1], 2), round(t.test(Q4.vec.qin)$conf.int[2], 2)),
                 c(round(mean(Q5.vec.qin),2), round(t.test(Q5.vec.qin)$conf.int[1], 2), round(t.test(Q5.vec.qin)$conf.int[2], 2)),
                 c(round(mean(Q6.vec.qin),2), round(t.test(Q6.vec.qin)$conf.int[1], 2), round(t.test(Q6.vec.qin)$conf.int[2], 2)),
                 c(round(MSE_par1_Qin, 2), NULL, NULL),
                 c(round(MSE_par2_Qin, 2), NULL, NULL),
                 c(round(MSE_p_Qin, 2), NULL, NULL))
    
    
    
    parameter <- c('alpha','beta', 'pi', 'mean', '.25Q', 'median', '.75Q', '.90Q', '.95Q', '.99Q', 'MSE alpha', 'MSE beta', 'MSE p', 
                   'alpha','beta', 'pi', 'mean', '.25Q', 'median', '.75Q', '.90Q', '.95Q', '.99Q', 'MSE alpha', 'MSE beta', 'MSE p')
    method <- c(rep('Deng', nrow(res) / 2), rep('Qin', nrow(res) / 2)) 
    
    colnames(res) <- c("estimate", "L_95%CI", "U_95%CI")
    
    res <- data.frame(parameter, method, res)
    
    result = rbind(result, res)
    print(cat("sample size:", m, "pi:", 1-pp))
    print(res)
  }
  return(result)
}

Gamma_DengQin_estimates <- incubationPeriod_Estimates_Deng(1000, distribution = "gamma", par1 = 5, par2 = 0.8)

Weibull_DengQin_estimates <- incubationPeriod_ParEstimates(1000, distribution = "weibull", par1 = 2, par2 = 8)

Lognormal_DengQin_estimates <- incubationPeriod_ParEstimates(1000, distribution = "lognormal", par1 = 1.8, par2 = 0.4)

