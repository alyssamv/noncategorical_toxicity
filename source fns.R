######################################################
############### Functions ############################
######################################################
#Function to obtain average toxicity scores at each dose level, which is the weighted sum of toxicity weights and the probability 
# of toxicity at each dose level. This computes equation #1.23 in the manuscript. 
get.thresh <- function(ntox, weights, TOX){
  tox.type.names <- paste("tox.type", seq(1:ntox), sep="")
  tox_type <- seq(from=0, to=4,by=1)
  tox.type <- matrix(rep(tox_type, each=ntox), ncol=ntox, byrow=T)
  tox.type <- as.data.frame(tox.type)
  possible_outcomes <- expand.grid(tox.type[,1:ntox]) #Permutes all possible toxicity profiles into a data set
  #map toxicity profiles to weights
  mapped.weight <- NA
  v <- NULL
  for(i in 1:ntox){
    for(k in 1:nrow(possible_outcomes)){
      mapped.weight[k] <- ifelse(possible_outcomes[k,i] > 0, weights[i,possible_outcomes[k,i]], 0)
    }
    v <- cbind(v, mapped.weight)
  }
  mapped_data <- matrix(v, ncol=3, byrow=F)
  scores <- apply(mapped_data, 1, sum) #Calulate toxicity scores for each profile in data set
  prob.data <- array(NA, c(nrow(possible_outcomes), ntox, d))
  for(j in 1:d){
    for(i in 1:nrow(possible_outcomes)){
      for(k in 1:ntox){
        tox.prob <- ifelse(possible_outcomes[i,k] > 0, TOX[j, possible_outcomes[i,k]+1, k], TOX[j,1,k])
        prob.data[i,k,j] <- tox.prob
      }
    }
  }
  prob.tox <- apply(prob.data[,,],c(1,3),prod)
  prob.scores <- scores*prob.tox
  thresh <- apply(prob.scores, 2, sum)
  return(thresh)
}


#Function to caculate the probabiity of a DLT at each dose level using 
#the individual toxicity probabilities across grades and doses
dlt.prob <- function(TOX, ntox, tox.dlt){
  cp <- function(p) 
  {
    ev <- do.call(expand.grid,replicate(length(p),0:1,simplify=FALSE)) # all permutations (rows) of DLT for each tox type
    pe <- apply(ev,1,function(x) prod(p*(x==1)+(1-p)*(x==0))) # Law of total prob ; DLT occurs with probability p 
    tapply(pe,rowSums(ev),sum) # sums probabilities for each group (total # DLTs across tox types)
  }
  ## probability of DLT for each dose at each tox type
  dlt.matrix <- matrix(0, nrow=nrow(TOX), ncol=d)
  for(i in 1:ntox){
    xxx <- TOX[,-(1:tox.dlt[i]),i] # for the dose and tox type, give only the cols corresponding to the grades that qualify as DLT
    xxx <- cbind(xxx, 0) # ensures that xxx is a matrix (instead of vector)
    dlt.matrix[i,] <- apply(xxx, 1, sum) # why sum?
  }
  ## probability of DLT for each dose across tox types
  ptox <- NULL
  for(i in 1:d){
    ptox[i] <- 1-cp(dlt.matrix[1:ntox,i])[1] #Ptox is probability of DLT at each dose level
  }
  return(ptox)
}

tox.weights <- function(W){
  thetamax <- sum(W[,4])
  W2 <- W/thetamax
  return(W2)
}