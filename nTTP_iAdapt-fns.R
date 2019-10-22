###############################################################
############### iAdapt for truncated Normal LR ################
###############################################################

# W weight matrix of toxicity burden for each type of toxicity
# TOX list of matrices giving probability of observing given toxicity grade for each toxicity type
# ntox number of toxicity types
# dose dose to simulate nTTP for

nTTP.sim <- function(W, TOX, ntox, dose){
  Tox <- NA
  random.tox <- runif(ntox) #draw random number from uniform(0,1)
  comb.curr <- dose #Define current dose level
  
  W.1 <- W / sum(W[,4]) # normalize toxicity burden matrix (equivalent of tox.weights from Ezzalfani)
  
  for (k in 1:ntox) {
    Tox[k] <- ifelse(random.tox[k] < TOX[comb.curr, 1, k] || TOX[comb.curr,1,k] == 0, 0, 
                     ifelse(random.tox[k] < TOX[comb.curr, 2, k] + TOX[comb.curr, 1, k], 1, 
                            ifelse(random.tox[k] < TOX[comb.curr, 2, k] + TOX[comb.curr, 1, k] + TOX[comb.curr, 3, k], 2,
                                   ifelse(random.tox[k] < TOX[comb.curr, 1, k] + TOX[comb.curr, 2, k] + TOX[comb.curr, 3, k] + TOX[comb.curr, 4, k], 3, 4))))
  }
  #Tox #vector of observed toxicity's grades for each toxicity given the current dose
  toxscores <- NA
  for (k in 1:ntox) {
    toxscores[k] <- ifelse(max(Tox[k]) == 0, 0, W.1[k, max(Tox[k])]) 
  }
  
  nTTP <- sum(toxscores) #The observed toxicity score the patient 
  return(nTTP)
}






# dose  number of doses to be tested (scalar)
# dose.tox  vector of true toxicities for each dose. Values range from 0 - 1.
# p1  toxicity under null (unsafe DLT rate). Values range from 0 - 1.
# p2  toxicity under alternative (safe DLT rate). Values range from 0 - 1; p1 > p2
# K  threshold for LR. Takes integer values: 1,2,...(recommended K=2)
# coh.size  cohort size (number of patients) per dose (Stage 1) 

tox.profile.nTTP <- function(dose, p1, p2, K, coh.size, W, TOX, ntox){ 
  
  dose   <- c(1:dose)                                           # vector of counts up to number of doses given
  stop   <- 0
  cohort <- 0
  i      <- 1
  x      <- c()
  
  # bounds for nTTP (truncated normal distribution)
  a = 0
  b = 1
  
  while ((stop == 0) & (i <= length(dose))) {
    cohort <- cohort + 1                                        # current cohort corresponding to dose
    dltsi  <- replicate(coh.size, nTTP.sim(W = W, 
                                           TOX = TOX, 
                                           ntox = ntox, 
                                           dose = dose[i]))	# nTTPs for that dose based on tox prob
    
    l.p2   <- prod(sapply(dltsi, FUN = function(i){ dnorm((i - p2)/sigma) })) / (sigma*(pnorm((b - p2)/sigma) - pnorm((a - p2)/sigma)))^coh.size # likelihood of acceptable/alternative hypothesis 
    l.p1   <- prod(sapply(dltsi, FUN = function(i){ dnorm((i - p1)/sigma) })) / (sigma*(pnorm((b - p1)/sigma) - pnorm((a - p1)/sigma)))^coh.size # likelihood of unacceptable/null hypothesis
    LR     <- round(l.p2/l.p1, 2)    
    
    x <- c(x, dose[i], mean(dltsi), cohort, LR)                              			
    
    if (LR <= (1/K))                        # stop escalation
      stop <- 1
    
    if (LR > (1/K))                         # escalate to next dose
      i <- i + 1	 
    
  }       
  return(matrix(x, ncol = 4, byrow = TRUE))
} 
