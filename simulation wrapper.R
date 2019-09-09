source(file.path(getwd(), "source fns.R"))

ntox <- 3 #Number of unique toxicities
d <- 6 #Number of dose levels

#### Define the weight Matrix ####
W <- matrix(c(0.1, 0.35, 0.7, 1.00, #Burden weight for grades 1-4 for toxicity 1
              0.08, 0.23,  0.6, 0.80, #Burden weight for grades 1-4 for toxicity 2
              0.00, 0.15, 0.45, 0.80 ##Burden weight for grades 1-4 for toxicity 3
), nrow=ntox, byrow=T)

#################################################################
####  Adjust W to normalized range of [0,1]              ########
#### This ensures the highest toxicity score wont be > 1 ########
#################################################################
for(i in 1:ntox){
  W <- tox.weights(W) 
}

#Define array to hold toxicitiy probabilities 
TOX <- array(NA, c(d, 5, ntox)) 

#### DEFINE ALL THE PROBABILITY SCENARIOS FOR 3x2 drug combo ######
#Scenario 1 = from Ezzalfani et al. paper
TOX[,,1] <- matrix(c(0.791, 0.172, 0.032, 0.004, 0.001, #probability of grades 0-4 toxicities for dose 1 and toxicity 1 in scenario 1
                     0.738, 0.195, 0.043, 0.015, 0.009, #probability of grades 0-4 toxicities for dose 2 and toxicity 1 in scenario 1
                     0.685, 0.190, 0.068, 0.044, 0.013, #probability of grades 0-4 toxicities for dose 3 and toxicity 1 in scenario 1
                     0.662, 0.200, 0.078, 0.046, 0.014, #probability of grades 0-4 toxicities for dose 4 and toxicity 1 in scenario 1
                     0.605, 0.223, 0.081, 0.071, 0.020, #probability of grades 0-4 toxicities for dose 5 and toxicity 1 in scenario 1
                     0.390, 0.307, 0.201, 0.074, 0.028),  #probability of grades 0-4 toxicities for dose 6 and toxicity 1
                   nrow=6, byrow=T)

TOX[,,2] <- matrix(c(0.968, 0.029, 0.002, 0.001, 0.000, #probability of grades 0-4 toxicities for dose 1 and toxicity 2 in scenario 1
                     0.763, 0.192, 0.006, 0.039, 0.000, #probability of grades 0-4 toxicities for dose 2 and toxicity 2 in scenario 1
                     0.762, 0.183, 0.041, 0.010, 0.004, #probability of grades 0-4 toxicities for dose 3 and toxicity 2 in scenario 1
                     0.682, 0.195, 0.108, 0.010, 0.005, #probability of grades 0-4 toxicities for dose 4 and toxicity 2 in scenario 1
                     0.397, 0.258, 0.276, 0.061, 0.008, #probability of grades 0-4 toxicities for dose 5 and toxicity 2 in scenario 1
                     0.260, 0.377, 0.281, 0.073, 0.009),  #probability of grades 0-4 toxicities for dose 6 and toxicity 2 in scenario 1
                   nrow=6, byrow=T)
TOX[,,3] <- matrix(c(0.907, 0.080, 0.008, 0.000, 0.005, #probability of grades 0-4 toxicities for dose 1 and toxicity 3 in scenario 1
                     0.602, 0.281, 0.035, 0.040, 0.042, #probability of grades 0-4 toxicities for dose 2 and toxicity 3 in scenario 1
                     0.536, 0.208, 0.031, 0.091, 0.134, #probability of grades 0-4 toxicities for dose 3 and toxicity 3 in scenario 1
                     0.015, 0.134, 0.240, 0.335, 0.276, #probability of grades 0-4 toxicities for dose 4 and toxicity 3 in scenario 1
                     0.005, 0.052, 0.224, 0.372, 0.347, #probability of grades 0-4 toxicities for dose 5 and toxicity 3 in scenario 1
                     0.004, 0.022, 0.223, 0.345, 0.406),  #probability of grades 0-4 toxicities for dose 6 and toxicity 3 in scenario 1
                   nrow=6, byrow=T)

#Obtain the true mean toxicity score at each dose level
get.thresh(3, W, TOX)

#Define at which grade of each toxicity is a DLT based on the traditional definition
dlts <- c(3,3,4) #i.e. A DLT is defined as a grade 3 of toxicity 1, a grade 3 for toxicity 2, and a grade 4 for toxicity 3

#Get the true probability of a least 1 DLT at each dose level
dlt.prob(TOX[,,], 3, dlts) 



################################################
###### To simualte data in a trial #############
################################################

N = 18
n.dose = 6
coh.size = 3

d.alloc = sort(rep(1:n.dose, N/n.dose))

Tox <- NA #Initializes vector to hold toxicity outcomes for single patient
Y.tox <- list() #Initializes list object to hold toxicity outcomes for all patients
Y.grade <- list()

#Simulate the three toxicities and their respective grades for a patient at the current dose level
for (j in 1:N){ 
  random.tox <- runif(ntox) #draw random number from uniform(0,1)
  comb.curr <- d.alloc[j] #Define current dose level
  
  for(k in 1:ntox){
    Tox[k] <- ifelse(random.tox[k] < TOX[comb.curr, 1, k] || TOX[comb.curr,1,k]==0, 0, 
                     ifelse(random.tox[k] < TOX[comb.curr, 2,k]+TOX[comb.curr,1,k], 1, 
                            ifelse(random.tox[k] < TOX[comb.curr,2,k]+TOX[comb.curr,1,k]+TOX[comb.curr, 3,k], 2,
                                   ifelse(random.tox[k] < TOX[comb.curr,1,k]+TOX[comb.curr,2,k]+TOX[comb.curr, 3,k]+TOX[comb.curr,4,k], 3, 4))))
  }
  #Tox #vector of observed toxicity's grades for each toxicity given the current dose
  toxscores <- NA
  for(k in 1:ntox){
    toxscores[k] <- ifelse(max(Tox[k]) == 0, 0, W[k,max(Tox[k])]) 
  }
  
  toxscore <- sum(toxscores) #The observed toxicity score the patient 
  #toxscore
  
  Y.grade[[j]] <- Tox
  Y.tox[[j]] <- toxscore
}


# calculate mean nTTP for each dose level
nTTP.bar = NA
for (i in 1:(N/coh.size)) {
    index = (i - 1)*coh.size
    nTTP.bar[i] = mean(unlist(Y.tox[index+(1:coh.size)]))
}

nTTP.bar
