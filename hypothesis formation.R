
source(file.path(getwd(), "source fns.R"))


ntox <- 3 #Number of unique toxicities
#### Define the weight Matrix ####
W <- matrix(c(0, 0.10, 0.35, 0.70, 1.00, #Burden weight for grades 0-4 for toxicity 1
              0, 0.08, 0.23, 0.60, 0.80, #Burden weight for grades 0-4 for toxicity 2
              0, 0.00, 0.15, 0.45, 0.80), ##Burden weight for grades 0-4 for toxicity 3
            nrow = ntox, byrow = T)

thetamax = sum(W[,ncol(W)])
# for (i in 1:ntox) {
#   W <- tox.weights(W) 
# }


### Calculate all possible nTTP values ###
### DLT defined as grade >= 3 in any toxicity category ###
tt = list(safe = c(),
          tox = c())
for (i in 1:5) {
  for (j in 1:5) {
    for (k in 1:5) {
      if (i >= 4 || j >= 4 || k >= 4 ){
        tt[["tox"]] = append(tt[["tox"]], 
                             sqrt(sum(W[1, i]^2, W[2, j]^2, W[3, k]^2))/thetamax)
      } else {
        tt[["safe"]] = append(tt[["safe"]], 
                              sqrt(sum(W[1, i]^2, W[2, j]^2, W[3, k]^2))/thetamax)
      }
    }
  }
}


# nTTP for safe and toxic (>= 1 DLT) 
lapply(tt, sort)
(max(tt[[1]]) + min(tt[[2]])) / 2 # target nTTP
lapply(tt, mean) # hypothesis values


par(mfrow = c(2, 1))
hist(tt[[1]], xlim = c(0, 0.8))
abline(v = mean(tt[[1]]), lty = 3, lwd = 2)
hist(tt[[2]], xlim = c(0, 0.8))
abline(v = mean(tt[[2]]), lty = 3, lwd = 2)
