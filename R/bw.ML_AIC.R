#'To Find Optimal Band Width Using K-Fold Maximum Likelihood Method
#'
#'@param data a numeric vector
#'@param k an integer greater than 1, by default `k = 2`
#'@param mplot by default `FALSE`, if TRUE then returns the mean log-likelihood vs. h plot
#'@param mkernel by default `"gaussian"`. It can be any of the kernel supported by `density` function
#'@return Optimal bandwidth with scale adjustment
#'@examples
#'install.packages("nor1mix")
#'library(nor1mix)  #need to install
#'ex1 <- norMix(mu = c(-2, 2),sigma = c(1, 0.8),w = c(.4, .6))
#'data1 <- rnorMix(300,ex1)  #generated data from Mixture Normal
#'h1 <- Optimal.KDE_ML(data1)
#'plot(density(data1,h1),col = "cornflowerblue")
#'@export
Optimal.KDE_ML <- function(data,k = 2,mplot = FALSE,mkernel = "gaussian"){
  
  if(k < 2){
    stop("k must be a finite integer greater than 1")
  }
  if(!is.integer(k)){
    stop("k must be a finite integer greater than 1")
  } 
  eta <- switch(mkernel,gaussian = 1,rectangular = 1/sqrt(3),
                triangular = 1/sqrt(6),epanechnikov = 1/sqrt(5),
                biweight = 1/sqrt(7),optcosine = sqrt(1-8/pi^2),
                cosine = sqrt(1/3 - 2/pi^2))
  
  m <- 1
  m <- ifelse(mkernel == "gaussian",3.5,m)
  
  split.data <- split(1:length(data),sample(1:length(data),length(data),replace = F)%%k)
  r <- extendrange(range(data))
  
  Estimated_Likelihood <- function(h){
    
    lik.vals <- lapply(split.data,FUN = function(test.data.index){  
      
      training.data <- data[-test.data.index]
      test.data <- data[test.data.index]
      fitted.density <- density(training.data,kernel = mkernel,bw = h*eta,
                                from = min(training.data) - m*h,to = max(training.data) + m*h,
                                n = max(2^(log2(ceiling(31*(max(training.data) - min(training.data) + 2*m*h)/(2*m*h) + 1))),
                                        2^10))
      approx.fitted.density <- approxfun(fitted.density$x,fitted.density$y)
      lik.val <- approx.fitted.density(test.data)
      lik.hh <- ifelse(is.na(lik.val) | (lik.val == 0),
                       min(lik.val[!(is.na(lik.val) | (lik.val == 0))]),lik.val)
      sum(log(lik.hh))
      
    })
    
    sum(unlist(lik.vals))/length(data)
  }  
  
  n <- mean(unlist(lapply(split.data,FUN = length)))
  hhrange <- diff(r) / n^c(0.9, 0.4)
  hh <- seq(hhrange[1], hhrange[2], length.out = 1500)
  
  mean_est.loglikelihood <-
    vapply(hh, FUN = Estimated_Likelihood,
           FUN.VALUE = 2)
  
  if(mplot){
    plot(hh, mean_est.loglikelihood,type = 'l',col = 'red',xlab = 'h',
         ylab = 'Mean Estimated log Likelihood',main = paste("K = ",k))
    
  }
  
  h.optimal <- hh[which.max(mean_est.loglikelihood)]
  
  return(h.optimal*eta)
}
#'To Find Optimal Band Width Using AIC
#'
#'@param data a numeric vector
#'@param mplot by default `FALSE`, if TRUE then returns the AIC vs. h plot
#'@param mkernel by default `"gaussian"`. It can be any of the kernel supported by `density` function
#'@return Optimal bandwidth with scale adjustment
#'@examples
#'install.packages("nor1mix")
#'library(nor1mix) # Need to install
#'ex2 <- norMix(mu = c(-2, 0,2),sigma = c(0.6,0.4,0.6),w = c(0.3,0.4,0.3))
#'data2 <- rnorMix(300,ex2)  #generated data from Mixture Normal
#'h2 <- Optimal.KDE_AIC(data2)
#'plot(density(data2,h2),col = "cornflowerblue")
#'@export
Optimal.KDE_AIC <- function(data,mplot = FALSE,mkernel = "gaussian"){ 
  
  eta <- switch(mkernel,gaussian = 1,rectangular = 1/sqrt(3),
                triangular = 1/sqrt(6),epanechnikov = 1/sqrt(5),
                biweight = 1/sqrt(7),optcosine = sqrt(1-8/pi^2),
                cosine = sqrt(1/3 - 2/pi^2))
  
  m <- 1
  m <- ifelse(mkernel == "gaussian",3.5,m)
  
  r <- extendrange(range(data))
  n <- length(data)
  
  h.AIC <- function(h) {
    
    fitted.density <- density(data,kernel = mkernel,bw = h*eta,
                              from = min(data) - m*h,to = max(data) + m*h,
                              n = max(2^(log2(ceiling(31*(max(data) - min(data) + 2*m*h)/(2*m*h) + 1))),
                                      2^10))
    approx.fitted.density <- approxfun(fitted.density$x,fitted.density$y)
    lik.val <- approx.fitted.density(data)
    lik.hh <- ifelse(is.na(lik.val) | (lik.val == 0),
                     min(lik.val[!(is.na(lik.val) | (lik.val == 0))]),lik.val)
    
    return(- 2*sum(log(lik.hh)) + 2*(diff(range(data))/h))
  }
  
  hhrange <- diff(r) / n^c(0.9, 0.4)
  hh <- seq(hhrange[1], hhrange[2], length.out = 1500)
  
  AICVal <-
    vapply(hh, FUN = h.AIC,FUN.VALUE = 2)
  
  if(mplot){
    plot(hh,AICVal,type = 'l',col = 'red',xlab = 'h',
         ylab = 'AIC')
    
  }
  
  h.optimal <- hh[which.min(AICVal)]
  return(h.optimal*eta)
}
