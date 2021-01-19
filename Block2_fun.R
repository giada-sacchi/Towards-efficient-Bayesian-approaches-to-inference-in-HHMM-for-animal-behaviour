#### Packages & Scripts ####
if(!require(boot)) install.packages("boot")
library(boot)
if(!require(Rcpp)) install.packages("Rcpp")
library(Rcpp)
if(!require(RcppArmadillo)) install.packages("RcppArmadillo")
library(RcppArmadillo)   
if(!require(invgamma)) install.packages("invgamma")
library(invgamma)   # ingamma
if(!require(numbers)) install.packages("numbers")
library(numbers) # mod
sourceCpp("13253_2017_282_MOESM4_ESM.cpp")


##### MLLK #####
## MINUS LOG-LIKELIHOOD
## defining a function that computes the negative log-likelihood
## if third line is not "active", then par.vec should have natural parameters

mllk <- function(pn, data, ll.N, ul.N, fit=TRUE){
  # create a matrix of 1's where to write the mllk for each data point in each ul state
  ul.mllk = matrix(1, nrow=length(data), ncol=ll.N);
  # initialize lists for ul probs at state1 and state2
  ul.all.probs1 = ul.all.probs2 = vector("list");
  # retrieve natural parameters from the working ones
  #pn = pw2pn(par.vec=par.vec, ll.N=ll.N, ul.N=ul.N);
  
  for(i in 1:ll.N){  ## for all the production states (ll)
    for(j in 1:length(data)){ ## for all the data points
      
      # initialize a prob matrix of 1's, with as many rows as data and columns as upper levels
      all.probs = matrix(1, nrow=nrow(data[[j]]), ncol=ul.N);
      # store the indices of non-missing data
      dd.ind = which(!is.na(data[[j]]$dive_duration));
      md.ind = which(!is.na(data[[j]]$maximum_depth));
      dw.ind = which(!is.na(data[[j]]$dive_wiggliness));
      # "zero true": store the indices of data for which wiggliness value is null but not missing
      zt.ind = which(data[[j]]$dive_wiggliness==0&!is.na(data[[j]]$dive_wiggliness));
      # "zero false": store the indices of data for which wiggliness value is neither null nor missing
      zf.ind = which(data[[j]]$dive_wiggliness!=0&!is.na(data[[j]]$dive_wiggliness));
      
      for(k in 1:ul.N){ ## for all the upper level states
        # initialize a prob matrix for the variables "duration" and "max depth" with all 1's and as many rows as data
        dd.probs = md.probs = dw.probs = rep(1, nrow(data[[j]]));         
        # get the prob for non-missing values from a gamma distribution
        dd.probs[dd.ind] = dgamma(data[[j]]$dive_duration[dd.ind], # data for which the variable value is not missing
                                  shape=pn$dd.mu[k]^2/pn$dd.sigma[k]^2, 
                                  scale=pn$dd.sigma[k]^2/pn$dd.mu[k]);
        # same steps as for dd.probs
        md.probs[md.ind] = dgamma(data[[j]]$maximum_depth[md.ind], 
                                  shape=pn$md.mu[k]^2/pn$md.sigma[k]^2, 
                                  scale=pn$md.sigma[k]^2/pn$md.mu[k]);
        # deal with variable "dw" ("wiggliness")
        if(length(zt.ind)>0){ ## if there are some null (BUT not missing) values for "wiggliness"
          # set the prob of the null values to the corresponding estimated proportion
          dw.probs[zt.ind] = pn$dw.pi[k];
          # set the prob for non-null values to the prob of the intersection of independent events
          # (prob of being non-null)*(prob of being non-missing)
          dw.probs[zf.ind] = (1-pn$dw.pi[k])*dgamma(data[[j]]$dive_wiggliness[zf.ind], 
                                                    shape=pn$dw.mu[k]^2/pn$dw.sigma[k]^2, 
                                                    scale=pn$dw.sigma[k]^2/pn$dw.mu[k]);
        }else{ # if there are NO null values for "wiggliness"
          # prob equal to the intersection(/product) of the independent events of being non-null AND non-missing
          dw.probs[dw.ind] = (1-pn$dw.pi[k])*dgamma(data[[j]]$dive_wiggliness[dw.ind], 
                                                    shape=pn$dw.mu[k]^2/pn$dw.sigma[k]^2, 
                                                    scale=pn$dw.sigma[k]^2/pn$dw.mu[k]);
        }
        # set the total probs for state "k" to the product of the probs of all the variables (prob of having "that" observation)
        all.probs[,k] = dd.probs*md.probs*dw.probs;
      } ## for all the upper level states
      # compute mllk for observation "j" when in production state (ul) "i" 
      ul.mllk[j,i] = nLogLike_Rcpp(ap=all.probs, 
                                   gamma=pn$ul.gamma[[i]], 
                                   foo=pn$ul.delta[[i]]*all.probs[1,],
                                   n=nrow(all.probs));
      if(i==1){ # for the first internal state
        ul.all.probs1[[j]] = all.probs;
      }
      if(i==2){ # for the second internal state
        ul.all.probs2[[j]] = all.probs; 
      }
    } ## for all the data points
  } ## for all the production states (ll)
  # initialize the prob matrix for internal (ll) states to all 1's
  ll.all.probs = matrix(1, nrow=length(data), ncol=ll.N);
  for(i in 1:ll.N){ # for all the internal states (ll)
    c = 700;    # consider a large costant (otherwise values too close to 0)
    # calculate the prob for internal state "i" as the exponential of the negative mllk
    # NOTE: exp(-(-log(likelihood))) = likelihood ~ cond prob
    ll.all.probs[,i] = exp(-ul.mllk[,i]+c); 
  }
  if(fit==TRUE){ 
    return(nLogLike_Rcpp(ap=ll.all.probs, 
                         gamma=pn$ll.gamma, 
                         foo=pn$ll.delta*ll.all.probs[1,], 
                         n=nrow(ll.all.probs))+length(data)*c);
  }else{
    foo = vector("list");     # create an empty list for ul probs
    foo[[1]] = ul.all.probs1; # first prod state
    foo[[2]] = ul.all.probs2; # second prod state
    list(ll.all.probs=ll.all.probs, ul.all.probs=foo);
  }
}


##### UPDATEPAR.BLOCK #####
# Function for updating the parameter vector through block updating (Version 2)
updatepar.block2 <- function(par.list, npar, data, n.iter, tuning=F, thr=c(0.2,0.4),
                             likhood,       # current likelihood
                             mu, sigma, alpha, beta, shape, rate, n.tune=n.iter,
                             delta){        # matrix with sd of random walk proposal distribution 
  # Define an array to store sample from posterior distribution for each chain
  itns <- array(0, dim=c(npar, n.iter))
  # counts of the accepted moves for each parameter (rows=parameters, cols=prod states)
  counts <- array(0, dim=c(11,3, n.iter))
  for (t in 1:n.iter){ ## MH iterations
      # Keep a record of the current parameter value being updated 
      oldpar <- par.list[[11]]
      # Propose a new candidate value using a random walk from normal proposal 
      par.list[[11]] <- rnorm(3, par.list[[11]], delta[11,])
      if (all(par.list[[11]]>=0) & all(par.list[[11]]<=1)){ # dw.pi is a probability
        # Compute the log-likelihood (of the move)
        newlikhood <- -mllk(pn=par.list, data=data, ll.N=2, ul.N=3, fit=TRUE)
        if(is.nan(newlikhood)){
          A <- 0
        } else {
          # Include likelihood and (log) prior contributions (symmetric proposal distributions cancel out)
          num <- newlikhood + sum(dbeta(par.list[[11]], alpha, beta,log=TRUE))
          den <- likhood + sum(dbeta(oldpar, alpha, beta,log=TRUE))
          A <- min(1, exp(num-den))
        }
      }  else { # if parameter is outside range set the acceptance probability A to 0
        A <- 0 
      }
      if (runif(1) <= A) {  # If the gengerated random number is smaller than the acceptance prob A
        # Accept the move with probability A and store its likelihood
        likhood <- newlikhood
        counts[11,,t] <- 1
      } else {     
        # Reject move and return parameter value to previous value
        par.list[[11]] <- oldpar
      }
    for (j in c(5,7,9)){ # parameter is a mean
        # Keep a record of the current parameter value being updated 
        oldpar <- par.list[[j]]
        # Propose a new candidate value using a random walk from normal proposal 
        par.list[[j]] <- rnorm(3, par.list[[j]], delta[j,])
        if (par.list[[j]][1]>=0 & par.list[[j]][2]>=0 & par.list[[j]][3]>=0){
          # Compute the log-likelihood (of the move)
          newlikhood <- -mllk(pn=par.list, data=data, ll.N=2, ul.N=3, fit=TRUE)
          if(is.nan(newlikhood)){
            A <- 0
          } else {
          # Include likelihood and (log) prior contributions (symmetric proposal distributions cancel out)
          num <- newlikhood + sum(dlnorm(par.list[[j]], mu, sd,log=TRUE))
          den <- likhood + sum(dlnorm(oldpar, mu, sd,log=TRUE))
          A <- min(1, exp(num-den))
          } 
        }  else { # Otherwise set the acceptance probability A to 0
          A <- 0 
        }
        if (runif(1) <= A) {  # If the gengerated random number is smaller than the acceptance prob A
          # Accept the move with probability A and store its likelihood
          likhood <- newlikhood
          counts[j,,t] <- 1
        } else {     
          # Reject move and return parameter value to previous value
          par.list[[j]] <- oldpar
        }
    }
    for (j in c(6,8,10)) { # parameter is a sd
        # Keep a record of the current parameter value being updated 
        oldpar <- par.list[[j]]
        # Propose a new candidate value using a random walk from normal proposal 
        par.list[[j]] <- rnorm(3, par.list[[j]], delta[j,])
        if (par.list[[j]][1]>=0 & par.list[[j]][2]>=0 & par.list[[j]][3]>=0){
          # Compute the log-likelihood (of the move)
          newlikhood <- -mllk(pn=par.list, data=data, ll.N=2, ul.N=3, fit=TRUE)
          if(is.nan(newlikhood)){
            A <- 0
          } else {
            # Include likelihood and (log) prior contributions (symmetric proposal distributions cancel out)
            num <- newlikhood + sum(dinvgamma(par.list[[j]], shape, rate,log=TRUE))
            den <- likhood + sum(dinvgamma(oldpar, shape, rate,log=TRUE))
            A <- min(1, exp(num-den))
          }
        }  else { # Otherwise set the acceptance probability A to 0
          A <- 0 
        }
        if (runif(1) <= A) {  # If the gengerated random number is smaller than the acceptance prob A
          # Accept the move with probability A and store its likelihood
          likhood <- newlikhood
          counts[j,,t] <- 1
        } else {     
          # Reject move and return parameter value to previous value
          par.list[[j]] <- oldpar
        }
    }
      itns[,t] <- as.vector(unlist(par.list))
      if (tuning==T){
        if (t<n.tune){
          if (mod(t,100)==0){ # every 100 iterations
            #look at accepted moves for the previous t iterations
            foo <- (apply(counts, c(1,2), sum)/t)
            # if less than thr[1]% of the moves were accepted, reduce delta (only for that parameter)
            delta[foo < thr[1]] <- 0.9*delta[foo < thr[1]]
            # if more than thr[2]% of the moves were accepted, increase delta (only for that parameter)
            delta[foo > thr[2]] <- 1.1*delta[foo > thr[2]]
          }
        }
      }
  } ## MH iterations
  return(list(itns=itns, counts=counts, delta=delta))
}


