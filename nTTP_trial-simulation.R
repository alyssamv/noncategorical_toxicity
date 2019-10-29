source(file.path(getwd(), "nTTP_iAdapt-fns.R"))
library(tidyverse)

p1 = 0.35 # unsafe nTTP (null hypothesis)
p2 = 0.1 # acceptable nTTP (alternative hypothesis)
sigma = 0.15 # variance of nTTP on each dose, assumed constant and estimated form simulated nTTP data
tru_tox = c(0.015, 0.05, 0.1, 0.2, 0.27, 0.32) # true (approx) nTTP for each dose, taken from simulations
DLT_rate = c(0.01, 0.07, 0.15, 0.28, 0.40, 0.49) # scenario equivalent for binary tox

coh.size = 3 # number pts per dose
ntox <- 3 # Number of unique toxicities
d <- 6 # Number of dose levels

#### Define the weight Matrix
W <- matrix(c(0.10, 0.35, 0.70, 1.00, # Burden weight for grades 1-4 for toxicity 1
              0.08, 0.23, 0.60, 0.80, # Burden weight for grades 1-4 for toxicity 2
              0.00, 0.15, 0.45, 0.80), ## Burden weight for grades 1-4 for toxicity 3
            nrow = ntox, byrow = T)

#### Define array to hold toxicitiy probabilities 
TOX <- array(NA, c(d, 5, ntox)) 

# Scenario 1 = from Ezzalfani et al. paper
TOX[,,1] <- matrix(c(0.791, 0.172, 0.032, 0.004, 0.001, #probability of grades 0-4 toxicities for dose 1 and toxicity 1 in scenario 1
                     0.738, 0.195, 0.043, 0.015, 0.009, #probability of grades 0-4 toxicities for dose 2 and toxicity 1 in scenario 1
                     0.685, 0.190, 0.068, 0.044, 0.013, #probability of grades 0-4 toxicities for dose 3 and toxicity 1 in scenario 1
                     0.662, 0.200, 0.078, 0.046, 0.014, #probability of grades 0-4 toxicities for dose 4 and toxicity 1 in scenario 1
                     0.605, 0.223, 0.081, 0.071, 0.020, #probability of grades 0-4 toxicities for dose 5 and toxicity 1 in scenario 1
                     0.390, 0.307, 0.201, 0.074, 0.028),  #probability of grades 0-4 toxicities for dose 6 and toxicity 1
                   nrow = 6, byrow = T)

TOX[,,2] <- matrix(c(0.968, 0.029, 0.002, 0.001, 0.000, #probability of grades 0-4 toxicities for dose 1 and toxicity 2 in scenario 1
                     0.763, 0.192, 0.006, 0.039, 0.000, #probability of grades 0-4 toxicities for dose 2 and toxicity 2 in scenario 1
                     0.762, 0.183, 0.041, 0.010, 0.004, #probability of grades 0-4 toxicities for dose 3 and toxicity 2 in scenario 1
                     0.682, 0.195, 0.108, 0.010, 0.005, #probability of grades 0-4 toxicities for dose 4 and toxicity 2 in scenario 1
                     0.397, 0.258, 0.276, 0.061, 0.008, #probability of grades 0-4 toxicities for dose 5 and toxicity 2 in scenario 1
                     0.260, 0.377, 0.281, 0.073, 0.009),  #probability of grades 0-4 toxicities for dose 6 and toxicity 2 in scenario 1
                   nrow = 6, byrow = T)

TOX[,,3] <- matrix(c(0.907, 0.080, 0.008, 0.000, 0.005, #probability of grades 0-4 toxicities for dose 1 and toxicity 3 in scenario 1
                     0.602, 0.281, 0.035, 0.040, 0.042, #probability of grades 0-4 toxicities for dose 2 and toxicity 3 in scenario 1
                     0.536, 0.208, 0.031, 0.091, 0.134, #probability of grades 0-4 toxicities for dose 3 and toxicity 3 in scenario 1
                     0.015, 0.134, 0.240, 0.335, 0.276, #probability of grades 0-4 toxicities for dose 4 and toxicity 3 in scenario 1
                     0.005, 0.052, 0.224, 0.372, 0.347, #probability of grades 0-4 toxicities for dose 5 and toxicity 3 in scenario 1
                     0.004, 0.022, 0.223, 0.345, 0.406),  #probability of grades 0-4 toxicities for dose 6 and toxicity 3 in scenario 1
                   nrow = 6, byrow = T)


###########################################################
############# iAdapt::tox.profile() for nTTP ##############
###########################################################

tox.profile.nTTP(dose = d, # number of doses
                 p1 = p1, # unsafe nTTP
                 p2 = p2, # acceptable nTTP
                 K = 2, # LR threshold
                 coh.size = coh.size, # n per dose
                 W = W, # weight matrix for TTP
                 TOX = TOX, # tox probabilities for each toxicity type
                 ntox = ntox) # number of toxicity types


###########################################################
########## simulation: compare nTTP and binary ############
###########################################################
K = 2
N = 1e3

## Using nTTP
stg1.nTTP <- list()
for (i in 1:N) {
  sim = tox.profile.nTTP(dose = d,
                         p1 = p1,
                         p2 = p2,
                         K = K,
                         coh.size = coh.size,
                         W = W,
                         TOX = TOX,
                         ntox = ntox)
  
  # is the dose safe (1) or unsafe (0)? (not including weak evidence)
  stg1.nTTP[[i]] = ifelse(sim[,4] >= K, 1, 0) 
  if (length(stg1.nTTP[[i]]) < d) {
    # designate for all 6 doses (allocation stops if dose declared unsafe)
    stg1.nTTP[[i]] = append(stg1.nTTP[[i]], rep(0, d - length(stg1.nTTP[[i]])))
  }
  
  # stg1[[i]]$w.weak = ifelse(sim[,4] >= K, 1, 
  #                        ifelse(sim[,4] <= 1/K, 0, NA))
  # 
}

## Binary tox equivalent
stg1 <- list()
for (i in 1:N) {
  sim = iAdapt::tox.profile(dose = d,
                            dose.tox = DLT_rate,
                            p1 = p1,
                            p2 = p2,
                            K = K,
                            coh.size = coh.size)
  
  # is the dose safe (1) or unsafe (0)? (not including weak evidence)
  stg1[[i]] = ifelse(sim[,4] >= K, 1, 0) 
  if (length(stg1[[i]]) < d) {
    # designate for all 6 doses (allocation stops if dose declared unsafe)
    stg1[[i]] = append(stg1[[i]], rep(0, d - length(stg1[[i]])))
  }
  
  # stg1[[i]]$w.weak = ifelse(sim[,4] >= K, 1, 
  #                        ifelse(sim[,4] <= 1/K, 0, NA))
  # 
}

print("Cohort of 3")
# % of times dose is declared safe with nTTP
rbind(tru_tox, colMeans(do.call(rbind, stg1.nTTP))) %>%
  `rownames<-`(c("True nTTP", "Rate declared safe")) %>%
  knitr::kable(row.names = TRUE, col.names = 1:6)
# % of times dose is declared safe for binary toxicity
rbind(DLT_rate, colMeans(do.call(rbind, stg1))) %>%
  `rownames<-`(c("True DLT rate", "Rate declared safe")) %>%
  knitr::kable(row.names = TRUE, col.names = 1:6)



###########################################################
####################### Cohort of 1 #######################
###########################################################

coh.size = 1 # number pts per dose
## Using nTTP
stg1.nTTP <- list()
for (i in 1:N) {
  sim = tox.profile.nTTP(dose = d,
                         p1 = p1,
                         p2 = p2,
                         K = K,
                         coh.size = coh.size,
                         W = W,
                         TOX = TOX,
                         ntox = ntox)
  
  # is the dose safe (1) or unsafe (0)? (not including weak evidence)
  stg1.nTTP[[i]] = ifelse(sim[,4] >= K, 1, 0) 
  if (length(stg1.nTTP[[i]]) < d) {
    # designate for all 6 doses (allocation stops if dose declared unsafe)
    stg1.nTTP[[i]] = append(stg1.nTTP[[i]], rep(0, d - length(stg1.nTTP[[i]])))
  }
  
  # stg1[[i]]$w.weak = ifelse(sim[,4] >= K, 1, 
  #                        ifelse(sim[,4] <= 1/K, 0, NA))
  # 
}

## Binary tox equivalent
stg1 <- list()
for (i in 1:N) {
  K = 1 # with cohort of 1, need K = 1
  sim = iAdapt::tox.profile(dose = d,
                            dose.tox = DLT_rate,
                            p1 = p1,
                            p2 = p2,
                            K = K, 
                            coh.size = coh.size)
  
  # is the dose safe (1) or unsafe (0)? (not including weak evidence)
  stg1[[i]] = ifelse(sim[,4] >= K, 1, 0) 
  if (length(stg1[[i]]) < d) {
    # designate for all 6 doses (allocation stops if dose declared unsafe)
    stg1[[i]] = append(stg1[[i]], rep(0, d - length(stg1[[i]])))
  }
  
  # stg1[[i]]$w.weak = ifelse(sim[,4] >= K, 1, 
  #                        ifelse(sim[,4] <= 1/K, 0, NA))
  # 
}

print("Cohort of 1")
# % of times dose is declared safe with nTTP
rbind(tru_tox, colMeans(do.call(rbind, stg1.nTTP))) %>%
  `rownames<-`(c("True nTTP", "Rate declared safe")) %>%
  knitr::kable(row.names = TRUE, col.names = 1:6)
# % of times dose is declared safe for binary toxicity
rbind(DLT_rate, colMeans(do.call(rbind, stg1))) %>%
  `rownames<-`(c("True DLT rate", "Rate declared safe")) %>%
  knitr::kable(row.names = TRUE, col.names = 1:6)


###########################################################
###################### Low toxicity #######################
###########################################################

p1 = 0.15 # unsafe nTTP (null hypothesis)
p2 = 0.05 # acceptable nTTP (alternative hypothesis)
coh.size = 3

## Using nTTP
stg1.nTTP <- list()
for (i in 1:N) {
  K = 2
  
  sim = tox.profile.nTTP(dose = d,
                         p1 = p1,
                         p2 = p2,
                         K = K,
                         coh.size = coh.size,
                         W = W,
                         TOX = TOX,
                         ntox = ntox)
  
  # is the dose safe (1) or unsafe (0)? (not including weak evidence)
  stg1.nTTP[[i]] = ifelse(sim[,4] >= K, 1, 0) 
  if (length(stg1.nTTP[[i]]) < d) {
    # designate for all 6 doses (allocation stops if dose declared unsafe)
    stg1.nTTP[[i]] = append(stg1.nTTP[[i]], rep(0, d - length(stg1.nTTP[[i]])))
  }
  
  # stg1[[i]]$w.weak = ifelse(sim[,4] >= K, 1, 
  #                        ifelse(sim[,4] <= 1/K, 0, NA))
  # 
}

## Binary tox equivalent
stg1 <- list()
for (i in 1:N) {
  K = 1
  
  sim = iAdapt::tox.profile(dose = d,
                            dose.tox = DLT_rate,
                            p1 = p1,
                            p2 = p2,
                            K = K,
                            coh.size = coh.size)
  
  # is the dose safe (1) or unsafe (0)? (not including weak evidence)
  stg1[[i]] = ifelse(sim[,4] >= K, 1, 0) 
  if (length(stg1[[i]]) < d) {
    # designate for all 6 doses (allocation stops if dose declared unsafe)
    stg1[[i]] = append(stg1[[i]], rep(0, d - length(stg1[[i]])))
  }
  
  # stg1[[i]]$w.weak = ifelse(sim[,4] >= K, 1, 
  #                        ifelse(sim[,4] <= 1/K, 0, NA))
  # 
}

print("Low toxicity")
# % of times dose is declared safe with nTTP
rbind(tru_tox, colMeans(do.call(rbind, stg1.nTTP))) %>%
  `rownames<-`(c("True nTTP", "Rate declared safe")) %>%
  knitr::kable(row.names = TRUE, col.names = 1:6)
# % of times dose is declared safe for binary toxicity
rbind(DLT_rate, colMeans(do.call(rbind, stg1))) %>%
  `rownames<-`(c("True DLT rate", "Rate declared safe")) %>%
  knitr::kable(row.names = TRUE, col.names = 1:6)


###########################################################
###################### High toxicity ######################
###########################################################

p1 = 0.5 # unsafe nTTP (null hypothesis)
p2 = 0.3 # acceptable nTTP (alternative hypothesis)
coh.size = 3

## Using nTTP
stg1.nTTP <- list()
for (i in 1:N) {
  K = 2
  
  sim = tox.profile.nTTP(dose = d,
                         p1 = p1,
                         p2 = p2,
                         K = K,
                         coh.size = coh.size,
                         W = W,
                         TOX = TOX,
                         ntox = ntox)
  
  # is the dose safe (1) or unsafe (0)? (not including weak evidence)
  stg1.nTTP[[i]] = ifelse(sim[,4] >= K, 1, 0) 
  if (length(stg1.nTTP[[i]]) < d) {
    # designate for all 6 doses (allocation stops if dose declared unsafe)
    stg1.nTTP[[i]] = append(stg1.nTTP[[i]], rep(0, d - length(stg1.nTTP[[i]])))
  }
  
  # stg1[[i]]$w.weak = ifelse(sim[,4] >= K, 1, 
  #                        ifelse(sim[,4] <= 1/K, 0, NA))
  # 
}

## Binary tox equivalent
stg1 <- list()
for (i in 1:N) {
  K = 2
  
  sim = iAdapt::tox.profile(dose = d,
                            dose.tox = DLT_rate,
                            p1 = p1,
                            p2 = p2,
                            K = K,
                            coh.size = coh.size)
  
  # is the dose safe (1) or unsafe (0)? (not including weak evidence)
  stg1[[i]] = ifelse(sim[,4] >= K, 1, 0) 
  if (length(stg1[[i]]) < d) {
    # designate for all 6 doses (allocation stops if dose declared unsafe)
    stg1[[i]] = append(stg1[[i]], rep(0, d - length(stg1[[i]])))
  }
  
  # stg1[[i]]$w.weak = ifelse(sim[,4] >= K, 1, 
  #                        ifelse(sim[,4] <= 1/K, 0, NA))
  # 
}

print("High toxicity")
# % of times dose is declared safe with nTTP
rbind(tru_tox, colMeans(do.call(rbind, stg1.nTTP))) %>%
  `rownames<-`(c("True nTTP", "Rate declared safe")) %>%
  knitr::kable(row.names = TRUE, col.names = 1:6)
# % of times dose is declared safe for binary toxicity
rbind(DLT_rate, colMeans(do.call(rbind, stg1))) %>%
  `rownames<-`(c("True DLT rate", "Rate declared safe")) %>%
  knitr::kable(row.names = TRUE, col.names = 1:6)
