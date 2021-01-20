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
# Function for the update of the whole parameter set
updatepar_tot <- function(par.list, npar, data, 
                      likhood,           # current likelihood
                      mu, sigma, alpha, beta, shape, rate,
                      acc.count=rapply(par.list, function(x)ifelse(x!=0,0,x), how="replace"),
                      delta){            # matrix with sd of random walk proposal distribution 
  for (j in 5:11){ # positions corresponding to target variable parameters
    for (i in 1:3) { # for the three levels of each variable (ul, ll)
      # Keep a record of the current parameter value being updated 
      oldpar <- par.list[[j]][i]
      # Propose a new candidate value using a random walk from normal proposal 
      par.list[[j]][i] <- rnorm(1, par.list[[j]][i], delta[[j]][i])
      if (j==11) { # dw.pi
        if (par.list[[j]][i]>=0 & par.list[[j]][i]<=1){ # dw.pi is a probability
          # Compute the log-likelihood (of the move)
          newlikhood <- -mllk(pn=par.list, data=data, ll.N=2, ul.N=3, fit=TRUE)
          if(is.nan(newlikhood)){
            A <- 0
          } else {
            # Include likelihood and (log) prior contributions (symmetric proposal distributions cancel out)
            num <- newlikhood + dbeta(par.list[[j]][i], alpha[i], beta[i], log=TRUE)
            den <- likhood + dbeta(oldpar, alpha[i], beta[i], log=TRUE)
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
          num <- newlikhood + dlnorm(par.list[[j]][i], mu[i], sd[i], log=TRUE)
          den <- likhood + dlnorm(oldpar, mu[i], sd[i], log=TRUE)
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
          num <- newlikhood + dinvgamma(par.list[[j]][i], shape[i], rate[i], log=TRUE)
          den <- likhood + dinvgamma(oldpar, shape[i], rate[i], log=TRUE)
          A <- min(1, exp(num-den))
          }
        }  else { # Otherwise set the acceptance probability A to 0
          A <- 0 
        }
      }
      if (runif(1) <= A) {  # If the gengerated random number is smaller than the acceptance prob A
        # Accept the move with probability A and store its likelihood
        likhood <- newlikhood
        acc.count[[j]][i] <-  acc.count[[j]][i] + 1
      } else {     
        # Reject move and return parameter value to previous value
        par.list[[j]][i] <- oldpar
      }
    }
  }
  for (j in 1:4){
    if (j==1){ # positions corresponding to target variable parameters
      # Keep a record of the current parameter value being updated 
      oldpar <- par.list[[j]]
      # Propose a new candidate value using a random walk from normal proposal 
      par.list[[j]][1] <- rnorm(1, par.list[[j]][1], delta[[j]][1])
      if ( par.list[[j]][1]>=0 &  par.list[[j]][1]<=1){ # probability
        par.list[[j]][2] <- 1 - par.list[[j]][1]
        # Compute the log-likelihood (of the move)
        newlikhood <- -mllk(pn=par.list, data=data, ll.N=2, ul.N=3, fit=TRUE)
        if(is.nan(newlikhood)){
          A <- 0
        } else {
          # Include likelihood and (log) prior contributions (symmetric proposal distributions cancel out)
          num <- newlikhood + sum(dbeta(par.list[[j]],2,2, log=TRUE))
          den <- likhood + sum(dbeta(oldpar,2,2, log=TRUE))
          A <- min(1, exp(num-den))
          if(is.nan(A)){A <- 0}
        }
      }  else { # if parameter is outside range set the acceptance probability A to 0
        A <- 0 
      }
      if (runif(1) <= A) {  # If the gengerated random number is smaller than the acceptance prob A
        # Accept the move with probability A and store its likelihood
        likhood <- newlikhood
        acc.count[[j]] <-  acc.count[[j]] + 1
      } else {     
        # Reject move and return parameter value to previous value
        par.list[[j]] <- oldpar
      }
    }
    if (j==2){ 
      for (k in 1:2){ # rows of the transition matrix
        # Keep a record of the current parameter value being updated 
        oldpar <- par.list[[j]][k,]
        # Propose a new candidate value using a random walk from normal proposal 
        par.list[[j]][k,1] <- rnorm(1, par.list[[j]][k,1], delta[[j]][k,1])
        if (par.list[[j]][k,1]>=0 & par.list[[j]][k,1]<=1){ # probability
          par.list[[j]][k,2] <- 1 - par.list[[j]][k,1]
          # Compute the log-likelihood (of the move)
          newlikhood <- -mllk(pn=par.list, data=data, ll.N=2, ul.N=3, fit=TRUE)
          if(is.nan(newlikhood)){
            A <- 0
          } else {
            # Include likelihood and (log) prior contributions (symmetric proposal distributions cancel out)
            num <- newlikhood + sum(dbeta(par.list[[j]][k,], 0.5, 0.5, log=TRUE))
            den <- likhood + sum(dbeta(oldpar, 0.5, 0.5, log=TRUE))
            A <- min(1, exp(num-den))
            if(is.nan(A)){A <- 0}
          }
        }  else { # if parameter is outside range set the acceptance probability A to 0
          A <- 0 
        }
        if (runif(1) <= A) {  # If the gengerated random number is smaller than the acceptance prob A
          # Accept the move with probability A and store its likelihood
          likhood <- newlikhood
          acc.count[[j]][k,] <- acc.count[[j]][k,] + 1
        } else {     
          # Reject move and return parameter value to previous value
          par.list[[j]][k,] <- oldpar
        }
      }
    }
    if (j==3){ # positions corresponding to target variable parameters
      for (k in 1:2){
        # Keep a record of the current parameter value being updated 
        oldpar <- par.list[[j]][[k]]
        # Propose a new candidate value using a random walk from normal proposal 
        par.list[[j]][[k]][1] <- rnorm(1, par.list[[j]][[k]][1], delta[[j]][[k]][1])
        if (par.list[[j]][[k]][1]>=0 & par.list[[j]][[k]][1]<=1){ # probability
          par.list[[j]][[k]][2] <- rnorm(1, par.list[[j]][[k]][2], delta[[j]][[k]][2])
          if ( par.list[[j]][[k]][2]>=0 &  par.list[[j]][[k]][2]<=1 & 
              (par.list[[j]][[k]][2] + par.list[[j]][[k]][1]) <=1){
            par.list[[j]][[k]][3] <- 1 - (par.list[[j]][[k]][1] + par.list[[j]][[k]][2])
            # Compute the log-likelihood (of the move)
            newlikhood <- -mllk(pn=par.list, data=data, ll.N=2, ul.N=3, fit=TRUE)
            if(is.nan(newlikhood)){
              A <- 0
            } else {
              # Include likelihood and (log) prior contributions (symmetric proposal distributions cancel out)
              num <- newlikhood + sum(dbeta(par.list[[j]][[k]], 0.5, 0.5, log=TRUE))
              den <- likhood + sum(dbeta(oldpar, 0.5, 0.5, log=TRUE))
              A <- min(1, exp(num-den))
              if(is.nan(A)){A <- 0}
            }
          }
        }  else { # if parameter is outside range set the acceptance probability A to 0
          A <- 0 
        }
        if (runif(1) <= A) {  # If the gengerated random number is smaller than the acceptance prob A
          # Accept the move with probability A and store its likelihood
          likhood <- newlikhood
          acc.count[[j]][[k]] <-  acc.count[[j]][[k]] + 1
        } else {     
          # Reject move and return parameter value to previous value
          par.list[[j]][[k]] <- oldpar
        }
      }
    }
    if (j==4){ 
      for (k in 1:2){ #length(par.list[[j]])
        for (i in 1:3){ # nrow(par.list[[j]][[k]]) # rows of the transition matrix
          # Keep a record of the current parameter value being updated 
          oldpar <- par.list[[j]][[k]][i,]
          # Propose a new candidate value using a random walk from normal proposal 
          par.list[[j]][[k]][i,1] <- rnorm(1, par.list[[j]][[k]][i,1], delta[[j]][[k]][i,1])
          if ( par.list[[j]][[k]][i,1]>=0 &  par.list[[j]][[k]][i,1]<=1){ # probability
            par.list[[j]][[k]][i,2] <- rnorm(1, par.list[[j]][[k]][i,2], delta[[j]][[k]][i,2])
            if ( par.list[[j]][[k]][i,2]>=0 &  par.list[[j]][[k]][i,2]<=1 & 
                (par.list[[j]][[k]][i,2] + par.list[[j]][[k]][i,1])<=1){
              par.list[[j]][[k]][i,3] <- 1 - (par.list[[j]][[k]][i,1] + par.list[[j]][[k]][i,2])
              # Compute the log-likelihood (of the move)
              newlikhood <- -mllk(pn=par.list, data=data, ll.N=2, ul.N=3, fit=TRUE)
              if(is.nan(newlikhood)){
                A <- 0
              } else {
                # Include likelihood and (log) prior contributions (symmetric proposal distributions cancel out)
                num <- newlikhood + sum(dbeta(par.list[[j]][[k]][i,], 0.5, 0.5, log=TRUE))
                den <- likhood + sum(dbeta(oldpar, 0.5, 0.5, log=TRUE))
                A <- min(1, exp(num-den))
                if(is.nan(A)){A <- 0}
              }
            }
          }  else { # if parameter is outside range set the acceptance probability A to 0
            A <- 0 
          }
          if (runif(1) <= A) {  # If the gengerated random number is smaller than the acceptance prob A
            # Accept the move with probability A and store its likelihood
            likhood <- newlikhood
            acc.count[[j]][[k]][i,] <- acc.count[[j]][[k]][i,] + 1
          } else {     
            # Reject move and return parameter value to previous value
            par.list[[j]][[k]][i,] <- oldpar
          }
        }
      }
    }
  } 
  return(list(par.list=par.list, likhood=likhood, acc.count=acc.count))
}
