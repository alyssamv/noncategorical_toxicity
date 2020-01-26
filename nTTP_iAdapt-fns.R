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





#### NOTE: no 'dose.tox' object needed for nTTP case ####

# dose  number of doses to be tested (scalar)
# p1  toxicity under null (unsafe DLT rate). Values range from 0 - 1.
# p2  toxicity under alternative (safe DLT rate). Values range from 0 - 1; p1 > p2
# K  threshold for LR. Takes integer values: 1,2,...(recommended K=2)
# coh.size  cohort size (number of patients) per dose (Stage 1) 
# W  weight matrix for toxicities
# TOX  probability of toxicity grades for each toxicity type
# ntox  number of toxicity types

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
    
    if (LR <= (1/K)) {                       # stop escalation
      stop <- 1
    } else if (LR > (1/K)) {                        # escalate to next dose
      i <- i + 1	 
    }
    
  }       
  return(matrix(x, ncol = 4, byrow = TRUE))
} 


# dose  number of doses to be tested (scalar)
# p1  toxicity under null (unsafe DLT rate). Values range from 0 - 1.
# p2  toxicity under alternative (safe DLT rate). Values range from 0 - 1; p1 > p2
# K  threshold for LR. Takes integer values: 1,2,...(recommended K=2)
# coh.size  cohort size (number of patients) per dose (Stage 1) 
# W  weight matrix for toxicities
# TOX  probability of toxicity grades for each toxicity type
# ntox  number of toxicity types

safe.dose.nTTP <- function(dose, p1, p2, K, coh.size, W, TOX, ntox) {
  
  res <- tox.profile.nTTP(dose, p1, p2, K, coh.size, W, TOX, ntox)
  alloc.total <- sort(rep(res[, 1], coh.size))
  n1 <- nrow(res) * coh.size
  unsafe.dose <- which(res[, 4] <= (1/K))
  if (length(unsafe.dose) == 0) {
    alloc.safe <- res[, 1:2]
  } else {
    alloc.safe <- res[res[, 1] != unsafe.dose, 1:2]
  }
  return(list(alloc.safe = alloc.safe, alloc.total = alloc.total, 
              n1 = n1))
}



eff.stg1.nTTP <- function(dose, p1, p2, K, coh.size, m, v, nbb = 100, W, TOX, ntox) {
  
  res <- safe.dose.nTTP(dose, p1, p2, K, coh.size, W, TOX, ntox)
  d.alloc <- res$alloc.total
  val.safe <- res$alloc.safe
  Y.safe <- d.safe <- tox.safe <- Y.alloc <- NULL
  n1 <- res$n1
  for (i in 1:length(d.alloc)) {
    ab <- beta.ab(m[d.alloc[i]]/100, v[d.alloc[i]])
    p <- stats::rbeta(1, ab$a, ab$b)
    Y.alloc[i] <- 100 * stats::rbinom(1, nbb, p)/nbb
  }
  if (length(val.safe) > 2) {
    d.safe <- sort(rep(val.safe[, 1], coh.size))
    tox.safe <- res$alloc.safe[, 2]
    Y.safe <- Y.alloc[1:length(d.safe)]
  } else if (length(val.safe) == 2) {
    d.safe <- sort(rep(val.safe[1], coh.size))
    tox.safe <- res$alloc.safe[2]
    Y.safe <- Y.alloc[1:length(d.safe)]
  } else {
    Y.safe <- d.safe <- NULL
    tox.safe <- res$alloc.safe[, 2]
  }
  return(list(Y.safe = Y.safe, d.safe = d.safe, tox.safe = tox.safe, 
              n1 = n1, Y.alloc = Y.alloc, d.alloc = d.alloc))
}


beta.ab <- function(m, v) {
  
  a <- seq(0.5, 20, 0.01)                            # a is a seq of alpha in beta distr.
  b <- a * (1 - m) / m
  
  vfit  <- a * b / ((a + b + 1) * (a + b)^2)
  diff  <- abs(vfit - v)
  index <- (1:length(diff))[diff == min(diff)]       # return the index of the min var.
  
  return(list(a = a[index],
              b = b[index]))                         # return alpha and beta for the min.var.
}



rand.stg2.nTTP <- function(dose, p1, p2, K, coh.size, m, v, N, stop.rule = 9, 
          cohort = 1, samedose = TRUE, nbb = 100, W, TOX, ntox) {
  
  res <- eff.stg1.nTTP(dose, p1, p2, K, coh.size, m, v, nbb, W, TOX, ntox)
  dose   <- c(1:dose) 
  yk.safe <- res$Y.safe
  yk.final <- res$Y.alloc
  dk.safe <- res$d.safe
  dk.final <- dk1 <- dk2 <- res$d.alloc
  toxk <- res$tox.safe
  n1 <- res$n1
  nmore <- N - n1
  nd <- length(unique(dk.safe))
  rp <- NULL
  stop <- 0
  
  if (nd == 0) {
    yk.final <- yk.final
    dk.final <- dk.final
    stop <- 1
  }
  
  if (nd == 1) {
    extra <- stop.rule - length(dk.safe)
    ab <- beta.ab(m[1]/100, v[1])
    y.extra <- 100 * stats::rbinom(extra, nbb, stats::rbeta(1, 
                                                            ab$a, ab$b))/nbb
    yk.final <- c(yk.final, y.extra)
    dk.final <- c(dk.final, rep(1, extra))
    stop <- 1
  }
  
  if (nd > 1) {
    coh.toxk <- cbind(matrix(dk.safe, ncol = coh.size, byrow = TRUE)[, 
                                                                     1], toxk)
    for (k in 1:nmore) {
      if (stop == 0) {
        reg <- stats::lm(log(yk.safe + 1) ~ factor(dk.safe))
        fit <- as.vector(reg$fitted.values)
        dose.unique <- duplicated(dk.safe)
        fitp <- exp(fit)
        fitp <- fitp[dose.unique == FALSE]
        rp <- fitp/sum(fitp)
        rp <- ifelse(rp < 0.02, 0.02, rp)
        dj <- stats::rmultinom(1, 1, prob = rp)
        if (samedose) {
          dj <- rep((1:length(dj))[dj == 1], cohort)
        } else {
          dosemat <- as.vector(dj * matrix(1:nd, ncol = cohort, 
                                           nrow = nd))
          dj <- dosemat[dosemat > 0]
        }
        ab <- beta.ab(m[dj]/100, v[dj])
        p <- stats::rbeta(1, ab$a, ab$b)
        yj <- 100 * stats::rbinom(1, nbb, p)/nbb
        toxj <- nTTP.sim(W = W, 
                         TOX = TOX, 
                         ntox = ntox, 
                         dose = dose[dj]) # stats::rbinom(1, size = 1, dose.tox[dj])
        coh.toxj <- c(dj, toxj)
        yk.safe <- c(yk.safe, yj)
        yk.final <- c(yk.final, yj)
        dk.safe <- c(dk.safe, dj)
        dk.final <- c(dk.final, dj)
        coh.toxk <- rbind(coh.toxk, coh.toxj)
        toxk <- c(toxk, toxj)
        n.obsk <- table(dk.safe)
        if (toxj == 0) {
          dk.safe <- dk.safe
          yk.safe <- yk.safe
        } else {
          LR.table.temp <- table(coh.toxk[, 1], coh.toxk[, 2])
          if (ncol(LR.table.temp) == 2) {
            LR.table <- cbind(LR.table.temp[, 2], n.obsk)
          } else {
            LR.table <- cbind(LR.table.temp[, 1], n.obsk)
          }
          loglik.p2 <- NULL
          loglik.p1 <- NULL
          lik.diff <- NULL
          accept.dose <- NULL
          for (j in 1:nrow(LR.table)) {
            loglik.p2[j] <- LR.table[j, 1] * log(p2) + 
              (LR.table[j, 2] - LR.table[j, 1]) * log(1 - 
                                                        p2)
            loglik.p1[j] <- LR.table[j, 1] * log(p1) + 
              (LR.table[j, 2] - LR.table[j, 1]) * log(1 - 
                                                        p1)
            lik.diff[j] <- exp(loglik.p2[j] - loglik.p1[j])
            accept.dose[j] <- ifelse(lik.diff[j] > (1/K), 
                                     1, 0)
          }
          dk.safe[dk.safe >= which(accept.dose == 0)] <- NA
          new.model <- cbind(dk.safe, yk.safe)
          new.model <- stats::na.omit(new.model)
          dk.safe <- new.model[, 1]
          yk.safe <- new.model[, 2]
          yk.final <- yk.final
          dk.final <- dk.final
          coh.toxk <- coh.toxk[!apply(coh.toxk, 1, function(x) {
            any(x >= which(accept.dose == 0))
          }), ]
        }
        if (length(unique(dk.safe)) > 1) {
          dk.safe <- dk.safe
          yk.safe <- yk.safe
          dk.final <- dk.final
          yk.final <- yk.final
        } else if (length(unique(dk.safe)) == 1) {
          new.size <- nmore + length(dk2)
          length.dk1 <- length(dk.final)
          if ((length(dk.safe) < stop.rule) && (length.dk1 < 
                                                new.size)) {
            extra.one <- min(new.size - length.dk1, 
                             stop.rule - length(dk.safe))
            ab <- beta.ab(m[1]/100, v[1])
            yj.one <- 100 * stats::rbinom(extra.one, 
                                          nbb, stats::rbeta(1, ab$a, ab$b))/nbb
            yk.final <- c(yk.final, yj.one)
            dk.final <- c(dk.final, rep(1, extra.one))
            stop <- 1
          } else {
            dk.final <- dk.final
            yk.final <- yk.final
            stop <- 1
          }
        }
        if (length(unique(dk.safe)) < 1) {
          dk.final <- dk.final
          yk.final <- yk.final
          stop <- 1
        }
      } else {
        dk.final <- dk.final
        yk.final <- yk.final
      }
    }
  }
  return(list(Y.final = yk.final, d.final = dk.final, n1 = n1))
}


sim.trials.nTTP <- function (numsims, dose, p1, p2, K, coh.size, m, v, 
                             N, stop.rule = 9, cohort = 1, samedose = TRUE, nbb = 100, W, TOX, ntox) 
{
  sim.yk <- sim.dk <- matrix(NA, nrow = numsims, ncol = N)
  sim.doses <- matrix(NA, nrow = numsims, ncol = dose)
  for (i in 1:numsims) {
    fstudy.out <- rand.stg2.nTTP(dose, p1, p2, K, coh.size, 
                            m, v, N, stop.rule, cohort, samedose, nbb, W, TOX, ntox)
    n.safe <- max(fstudy.out$d.final[(fstudy.out$n1 + 1):length(fstudy.out$d.final)], 
                  na.rm = TRUE)
    sim.doses[i, ] <- c(rep(1, n.safe), rep(0, dose - n.safe))
    if (length(fstudy.out$Y.final) < N) {
      sim.yk[i, ] <- c(fstudy.out$Y.final, rep(NA, N - 
                                                 length(fstudy.out$Y.final)))
      sim.dk[i, ] <- c(fstudy.out$d.final, rep(NA, N - 
                                                 length(fstudy.out$d.final)))
    } else {
      sim.yk[i, ] <- fstudy.out$Y.final
      sim.dk[i, ] <- fstudy.out$d.final
    }
    cat(i, "\n")
  }
  return(list(sim.Y = sim.yk, sim.d = sim.dk, safe.d = sim.doses))
}




sim.summary.new <- function(sims, print = TRUE){
  sim.doses = sims$sim.d
  n.doses = max(sim.doses, na.rm = TRUE)
  sim.eff = sims$sim.Y
  dose.mat.a <- matrix(NA, nrow(sim.doses), n.doses)
  for (i in 1:nrow(sim.doses)) {
    dose.no.na <- na.omit(sim.doses[i, ])
    dose.mat.a[i, ] <- table(factor(dose.no.na, levels = 1:n.doses))/length(dose.no.na)
  }
  est.dose1 <- matrix(NA, n.doses, 5)
  
  for (j in 1:n.doses) {
    est.dose1[j, ] <- c(j/100, 
                        quantile(dose.mat.a[, j], 
                                 prob = c(0.25, 0.5, 0.75), na.rm = TRUE),
                        round(mean(dose.mat.a[, j]), 2))
  }
  dose.IQR = round(est.dose1 * 100, 1)

  pers.hat.a <- matrix(NA, nrow(sim.eff), n.doses + 1)
  for (i in 1:nrow(sim.eff)) {
    for (j in 1:(n.doses + 1)) {
      pers.hat.a[i, j] <- (median(sim.eff[i, sim.doses[i, 
                                                       ] == j - 1]))
    }
  }
  est.pers1 <- matrix(NA, (n.doses + 1), 5)
  for (j in 1:(n.doses + 1)) {
    est.pers1[j, ] <- c((j - 1), 
                        quantile(pers.hat.a[, j], 
                                 prob = c(0.25, 0.5, 0.75), na.rm = TRUE),
                        round(mean(pers.hat.a[, j]), 2))
  }
  Y = est.pers1[-1, ]
  
  if (print == TRUE) {
    print(knitr::kable(dose.IQR, caption = "Percent allocation per dose level", 
                       col.names = c("Dose", "25th percentile", "Median", "75th percentile", "Mean")))
    print(knitr::kable(Y, caption = "Estimated efficacy per dose level", 
                       col.names = c("Dose", "25th percentile", "Median", "75th percentile", "Mean")))
  }

  return(list(pct.treated = dose.IQR, efficacy = Y))
}
