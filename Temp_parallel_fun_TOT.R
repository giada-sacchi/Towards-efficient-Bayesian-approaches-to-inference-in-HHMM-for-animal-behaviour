##### UPDATEPAR.TEMP #####
# Function used to update parameters (acceptance prob computed considering chain temperature)

updatepar.temp.tot <- function(par.list, npar, data, 
                      likhood,                     # current likelihood
                      B.temp=1,                    # temperature of the chain
                      E,                          # initial energy of the chain
                      mu, sigma, alpha, beta, shape, rate,
                      acc.count=rapply(par.list, function(x)ifelse(x!=0,0,x), how="replace"),
                      delta){            # matrix with sd of random walk proposal distribution 
  for (j in 5:11){ # positions corresponding to target variable parameters
    for (i in 1:3) { # for the three prod states
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
            num <- newlikhood + log(dbeta(par.list[[j]][i], alpha[i], beta[i]))
            den <- likhood + log(dbeta(oldpar, alpha[i], beta[i]))
            E <- den-num
            A <- min(1, exp(-B.temp*E))
          }
        }  else { # if parameter is outside range set the acceptance probability A to 0
          A <- 0 
        }
        A
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
          E <- den-num
          A <- min(1, exp(-B.temp*E))
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
          E <- den-num
          A <- min(1, exp(-B.temp*E))
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
          num <- newlikhood + sum(dbeta(par.list[[j]],2,2,log=TRUE))
          den <- likhood + sum(dbeta(oldpar,2,2,log=TRUE))
          E <- den-num
          A <- min(1, exp(-B.temp*E))
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
            num <- newlikhood + sum(dbeta(par.list[[j]][k,],0.5,0.5,log=TRUE))
            den <- likhood + sum(dbeta(oldpar,0.5,0.5,log=TRUE))
            E <- den-num
            A <- min(1, exp(-B.temp*E))
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
              num <- newlikhood + sum(dbeta(par.list[[j]][[k]],0.5,0.5,log=TRUE))
              den <- likhood + sum(dbeta(oldpar,0.5,0.5,log=TRUE))
              E <- den-num
              A <- min(1, exp(-B.temp*E))
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
                num <- newlikhood + sum(dbeta(par.list[[j]][[k]][i,],0.5,0.5,log=TRUE))
                den <- likhood + sum(dbeta(oldpar,0.5,0.5,log=TRUE))
                E <- den-num
                A <- min(1, exp(-B.temp*E))
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
  # Output the parameter values and log-likelihood value (and the vector of counts for accepted moves)
  output <- list(par.list=par.list, likhood=likhood, energy=E, counts=acc.count) 
  output
}


##### PRIOR #####
# Function to compute log-prior probability of a parameter vector

prior <- function(par.list, mu, sigma, alpha, beta, shape, rate){
  prior <- 0
  for (j in 5:11){ # positions corresponding to target variable parameters
    for (i in 1:3) {
      if (j==11) { 
        prior <- prior + log(dbeta(par.list[[j]][i], alpha[i], beta[i]))
      } else if (is.element(j, c(5,7,9))){
        prior <- prior + log(dlnorm(par.list[[j]][i], mu[i], sigma[i]))
      } else {  
        prior <- prior + log(dinvgamma(par.list[[j]][i], shape[i], rate[i]))
      }
    }
  }
  for (j in 1:4){
    if (j==1){
      prior <- prior + sum(log(dbeta(par.list[[j]],2,2)))
    }
    if (j==2){ 
      for (k in 1:2){ 
        prior <- prior + sum(log(dbeta(par.list[[j]][k,],.5,.5)))
      }
    }
    if (j==3){ 
      for (k in 1:2){
        prior <- prior + sum(log(dbeta(par.list[[j]][[k]],0.5,0.5)))
      }
    }
    if (j==4){ 
      for (k in 1:2){ #length(par.list[[j]])
        for (i in 1:3){ # nrow(par.list[[j]][[k]])
          prior <- prior + sum(log(dbeta(par.list[[j]][[k]][i,],0.5,0.5)))
        }
      }
    }
  }
  return(prior)
}


##### TEMP.TOT.UPDATE #####
temp.tot.update <- function(j){
  sourceCpp("13253_2017_282_MOESM4_ESM.cpp")
  B <- B.temp[j]
  E <- E[j]
  for (t in 1:n.iter){
    output <- updatepar.temp.tot(par.list=pn, npar=n.par, data=data, likhood=likelihood, 
                             shape=shape, rate=rate, alpha=alpha, beta=beta, delta=delta,
                             mu=mu, sigma=sd, B.temp=B, E=E, acc.count=acc.count)
    pn <- output$par.list 
    likelihood <- output$likhood
    acc.count <- output$counts
    E <- output$energy
    # save temporarily iterations for each chain
    itns.tmp[,t] <- as.vector(unlist(pn))
    if (t<n.tune){
      if (mod(t,100)==0){ # every 100 iterations
        foo <- as.vector(unlist(acc.count))/t
        delta <- as.vector(unlist(delta))
        delta_tmp <- ifelse(foo < thr[1], 0.9*delta,
                            ifelse(foo > thr[2], 1.1*delta, delta))
        delta <- list(ll.delta=delta[1:2],
                      ll.gamma=matrix(delta[3:6], nrow=2),
                      ul.delta=list(delta[7:9],delta[10:12]),
                      ul.gamma=list(matrix(delta[13:21],nrow=3),
                                    matrix(delta[22:30],nrow=3)),
                      dd.mu=delta_tmp[31:33],
                      dd.sigma=delta_tmp[34:36],
                      md.mu=delta_tmp[37:39],
                      md.sigma=delta_tmp[40:42],
                      dw.mu=delta_tmp[43:45],
                      dw.sigma=delta_tmp[46:48],
                      dw.pi=delta_tmp[49:51])
      }
    }
  }
  return(list(pn=pn, likelihood=likelihood, energy=E, itns=itns.tmp, B.temp=B,counts=acc.count))
}
