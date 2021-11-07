###Comparison foward and backward time

#proposed data generation method

sample_my_forwardtime <- NULL
for (i in 1:1000){
  sample_my_forwardtime <- append(sample_my_forwardtime, simulate_outbreak_wuhan(125, 245000, 490000, 0.14, "gamma", 5, 0.8, "forward"))
}

sample_my_backwarddtime <- NULL
for (i in 1:1000){
  sample_my_backwarddtime <- append(sample_my_backwarddtime, simulate_outbreak_wuhan(125, 245000, 490000, 0.14, "gamma", 5, 0.8, "backward"))
}

hist(unlist(sample_my_forwardtime), main = 'Forward time proposed generation method', xlim = c(0,20) , xlab = 'Forward time', breaks = 20, probability = TRUE)

hist(unlist(sample_my_backwarddtime), main = 'Backward time proposed generation method', xlim = c(0,20) , xlab = 'Backward time',  breaks = 20, probability = TRUE)

#Deng's and Qin's data generation method with Pi = 0.2

sample_Deng_forwardtime_0.2 <- NULL
for (i in 1:1000){
  sample_Deng_forwardtime_0.2 <- append(sample_Deng_forwardtime_0.2, GenerateData_Deng(1200, 0.2, 5, 0.8, "gamma", "forward"))
}

sample_Deng_backwardtime_0.2 <- NULL
for (i in 1:1000){
  sample_Deng_backwardtime_0.2 <- append(sample_Deng_backwardtime_0.2, GenerateData_Deng(1200, 0.2, 5, 0.8, "gamma", "backward"))
}

hist(unlist(sample_Deng_forwardtime_0.2), main = 'Forward time Deng and Qin, Pi = 0.2', xlim = c(0,20) , xlab = 'Forward time', breaks = 20, probability = TRUE)

hist(unlist(sample_Deng_backwardtime_0.2), main = 'Backward time Deng and Qin, Pi = 0.2', xlim = c(0,20) , xlab = 'Backward time', breaks = 20, probability = TRUE)


#Deng's and Qin's data generation method with Pi = 0.0

sample_Deng_forwardtime_0.0 <- NULL
for (i in 1:1000){
  sample_Deng_forwardtime_0.0 <- append(sample_Deng_forwardtime_0.0, GenerateData_Deng(1200, 0, 5, 0.8, "gamma", "forward"))
}

sample_Deng_backwardtime_0.0 <- NULL
for (i in 1:1000){
  sample_Deng_backwardtime_0.0 <- append(sample_Deng_backwardtime_0.0, GenerateData_Deng(1200, 0, 5, 0.8, "gamma", "backward"))
}

hist(unlist(sample_Deng_forwardtime_0.0), main = 'Forward time Deng and Qin, Pi = 0.0', xlim = c(0,20) , xlab = 'Forward time', breaks = 20, probability = TRUE)

hist(unlist(sample_Deng_backwardtime_0.0), main = 'Backward time Deng and Qin, Pi = 0.0', xlim = c(0,20) , xlab = 'Backward time', breaks = 20, probability = TRUE)


