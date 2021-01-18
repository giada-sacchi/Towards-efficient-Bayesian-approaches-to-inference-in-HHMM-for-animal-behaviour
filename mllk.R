## MINUS LOG-LIKELIHOOD
## Defining a function that computes the negative log-likelihood
## if third line is not "active", then par.vec should have natural parameters

mllk <- function(pn, data, ll.N, ul.N, fit=TRUE){
  
  # create a matrix of 1's where to write the mllk for each data point in each ul state
  ul.mllk = matrix(1, nrow=length(data), ncol=ll.N)
  # initialize lists for ul probs at state1 and state2
  ul.all.probs1 = ul.all.probs2 = vector("list")

  for(i in 1:ll.N){ ## for all the production states (ll)
    for(j in 1:length(data)){ ## for all the data points
      
      # initialize a prob matrix of 1's, with as many rows as data and columns as upper levels
      all.probs = matrix(1, nrow=nrow(data[[j]]), ncol=ul.N)
      # store the indices of non-missing data
      dd.ind = which(!is.na(data[[j]]$dive_duration))
      md.ind = which(!is.na(data[[j]]$maximum_depth))
      dw.ind = which(!is.na(data[[j]]$dive_wiggliness))
      # "zero true": store the indices of data for which wiggliness value is null but not missing
      zt.ind = which(data[[j]]$dive_wiggliness==0&!is.na(data[[j]]$dive_wiggliness))
      # "zero false": store the indices of data for which wiggliness value is neither null nor missing
      zf.ind = which(data[[j]]$dive_wiggliness!=0&!is.na(data[[j]]$dive_wiggliness))
      
      for(k in 1:ul.N){ ## for all the upper level states
        # initialize a prob matrix for the variables "duration" and "max depth" with all 1's and as many rows as data
        dd.probs = md.probs = dw.probs = rep(1, nrow(data[[j]]))         
        # get the prob for non-missing values from a gamma distribution
        dd.probs[dd.ind] = dgamma(data[[j]]$dive_duration[dd.ind], # data for which the variable value is not missing
                                  shape=pn$dd.mu[k]^2/pn$dd.sigma[k]^2, 
                                  scale=pn$dd.sigma[k]^2/pn$dd.mu[k])
        md.probs[md.ind] = dgamma(data[[j]]$maximum_depth[md.ind], 
                                  shape=pn$md.mu[k]^2/pn$md.sigma[k]^2, 
                                  scale=pn$md.sigma[k]^2/pn$md.mu[k])
        
        # deal with variable "dw" ("wiggliness")
        if(length(zt.ind)>0){ # if there are some null (BUT not missing) values for "wiggliness"
          # set the prob of the null values to the corresponding estimated proportion
          dw.probs[zt.ind] = pn$dw.pi[k]
          # set the prob for non-null values to the prob of the intersection of independent events
          # (prob of being non-null)*(prob of being non-missing)
          dw.probs[zf.ind] = (1-pn$dw.pi[k])*dgamma(data[[j]]$dive_wiggliness[zf.ind], 
                                                    shape=pn$dw.mu[k]^2/pn$dw.sigma[k]^2, 
                                                    scale=pn$dw.sigma[k]^2/pn$dw.mu[k])
        }else{ # if there are NO null values for "wiggliness"
          # prob equal to the intersection(/product) of the independent events of being non-null AND non-missing
          dw.probs[dw.ind] = (1-pn$dw.pi[k])*dgamma(data[[j]]$dive_wiggliness[dw.ind], 
                                                    shape=pn$dw.mu[k]^2/pn$dw.sigma[k]^2, 
                                                    scale=pn$dw.sigma[k]^2/pn$dw.mu[k])
        }
        # set the total probs for state "k" to the product of the probs of all the variables (prob of having "that" observation)
        all.probs[,k] = dd.probs*md.probs*dw.probs
      } ## for all the upper level states
      
      # compute mllk for observation "j" when in production state (ul) "i" 
      ul.mllk[j,i] = nLogLike_Rcpp(ap=all.probs, 
                                   gamma=pn$ul.gamma[[i]], 
                                   foo=pn$ul.delta[[i]]*all.probs[1,],
                                   n=nrow(all.probs))
      if(i==1){ # for the first internal state
        ul.all.probs1[[j]] = all.probs
      }
      if(i==2){ # for the second internal state
        ul.all.probs2[[j]] = all.probs
      }
      
    } ## for all the data points
  } # for all the production states (ll)
  
  # initialize the prob matrix for internal (ll) states to all 1's
  ll.all.probs = matrix(1, nrow=length(data), ncol=ll.N)
  
  for(i in 1:ll.N){ # for all the internal states (ll)
    c = 700;    # consider a large constant (otherwise values too close to 0)
    # calculate the prob for internal state "i" as the exponential of the negative mllk
    # NOTE: exp(-(-log(likelihood))) = likelihood ~ cond prob
    ll.all.probs[,i] = exp(-ul.mllk[,i]+c); 
  }
  
  if(fit==TRUE){ 
    return(nLogLike_Rcpp(ap=ll.all.probs, 
                         gamma=pn$ll.gamma, 
                         foo=pn$ll.delta*ll.all.probs[1,], 
                         n=nrow(ll.all.probs))+length(data)*c)
  }else{
    foo = vector("list");     ## create an empty list for ul probs
    foo[[1]] = ul.all.probs1; ## first prod state
    foo[[2]] = ul.all.probs2; ## second prod state
    list(ll.all.probs=ll.all.probs, ul.all.probs=foo) ## output all ll and ul probs
  }
}
