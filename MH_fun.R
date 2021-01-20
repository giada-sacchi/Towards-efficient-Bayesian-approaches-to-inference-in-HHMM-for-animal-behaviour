#### Packages & Scripts ####
if(!require(boot)) install.packages("boot")
library(boot)
if(!require(Rcpp)) install.packages("Rcpp")
library(Rcpp)
if(!require(RcppArmadillo)) install.packages("RcppArmadillo")
library(RcppArmadillo)   
if(!require(invgamma)) install.packages("invgamma")
library(invgamma) 
if(!require(numbers)) install.packages("numbers")
library(numbers) # mod
sourceCpp("13253_2017_282_MOESM4_ESM.cpp")


##### UPDATEPAR #####
# Function for single-update MH
updatepar_MH <- function(par.list, npar, data, 
                      likhood,           # current likelihood
                      mu, sigma,        
                      alpha, beta,
                      shape, rate,
                      delta){            # matrix with sd of random walk proposal distribution 
  # counts of the accepted moves for each parameter (rows=parameters, cols=prod states)
  counts <- matrix(rep(0,33), ncol=3)
  for (j in 5:11){ # positions corresponding to target variable parameters
    for (i in 1:3) { # for the three levels of each variable (ul, ll)
      # Keep a record of the current parameter value being updated 
      oldpar <- par.list[[j]][i]
      # Propose a new candidate value using a random walk from normal proposal 
      par.list[[j]][i] <- rnorm(1, par.list[[j]][i], delta[j,i])
      if (j==11) { # dw.pi
        if (par.list[[j]][i]>=0 & par.list[[j]][i]<=1){ # dw.pi is a probability
          # Compute the log-likelihood (of the move)
          newlikhood <- -mllk(pn=par.list, data=data, ll.N=2, ul.N=3, fit=TRUE)
          if(is.nan(newlikhood)){
            A <- 0
          } else {
            # Include likelihood and (log) prior contributions (symmetric proposal distributions cancel out)
            num <- newlikhood + log(dbeta(par.list[[j]][i], alpha[i], beta[i]))
            den <- likhood + log(dbeta(oldpar, alpha[i], beta[i]))
            A <- min(1, exp(num-den))
          }
        }  else { # if parameter is outside range set the acceptance probability A to 0
          A <- 0 
        }
      } else if (is.element(j,c(5,7,9))){ # parameter is a mean
        if (par.list[[j]][i]>=0){
          # Compute the log-likelihood (of the move)
          newlikhood <- -mllk(pn=par.list, data=data, ll.N=2, ul.N=3, fit=TRUE)
          if(is.nan(newlikhood)){
            A <- 0
          } else {
          # Include likelihood and (log) prior contributions (symmetric proposal distributions cancel out)
          num <- newlikhood + log(dlnorm(par.list[[j]][i], mu[i], sd[i]))
          den <- likhood + log(dlnorm(oldpar, mu[i], sd[i]))
          A <- min(1, exp(num-den))
          }
        }  else { # Otherwise set the acceptance probability A to 0
          A <- 0 
        }
      } else { # parameter is a sd
        if (par.list[[j]][i]>=0){
          # Compute the log-likelihood (of the move)
          newlikhood <- -mllk(pn=par.list, data=data, ll.N=2, ul.N=3, fit=TRUE)
          if(is.nan(newlikhood)){
            A <- 0
          } else {
          # Include likelihood and (log) prior contributions (symmetric proposal distributions cancel out)
          num <- newlikhood + log(dinvgamma(par.list[[j]][i], shape[i], rate[i]))
          den <- likhood + log(dinvgamma(oldpar, shape[i], rate[i]))
          A <- min(1, exp(num-den))
          }
        }  else { # Otherwise set the acceptance probability A to 0
          A <- 0 
        }
      }
      if (runif(1) <= A) {  # If the generated random number is smaller than the acceptance prob A
        # Accept the move with probability A and store its likelihood
        likhood <- newlikhood
        counts[j,i] <- 1
      } else {     
        # Reject move and return parameter value to previous value
        par.list[[j]][i] <- oldpar
      }
    }
  }
  return(list(par.list=par.list, likhood=likhood, counts=counts))
}
