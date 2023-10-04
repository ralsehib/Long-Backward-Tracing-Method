

# This file contains functions to estimate the infection rate and
# the recovery rate in a partially known network using the Long Back-Tracing fill-up method.
 
# Last modified: 26 September 2023

# Authors: Omar De La Cruz Cabrera and Razan Alsehibani.



## Simulation parameters
n = 200   # Number of individuals
m = 200   # Number of time steps
q = 100   # Number of repetitions



## SIS parameters
beta  =  0.04 # Rate of infection
gamma =  0.02 # Rate of recovery


## UK parameters
alpha = 0.007 # Rate of contact at random
delta = 0.006 # Rate of contact through network
rao = alpha + delta # Rate of contact at random or through network 


library(igraph)
start_time<- Sys.time()


AB2<- readMM(file="./matrixData/AdjBAComplex_Power_1_.mtx") 
AE <-AB2*1

listBetaHatE <- vector("integer", q)
listGammaHatE <- vector("integer", q)

for (b in 1:q){
  print(b)
  sE = c(1,rep(0,n-1) )  # initial status (1: Infected, 0: Susceptible)
  kE = c(1, rep(0,n-1) ) # initial status (1: Known, 0: Unknown)
  ninE = c(1,rep(0,n-1)) # initial status for the neighbors (1: Infected, 0: Susceptible)
  
  KE = matrix(0,nrow = m,ncol = n)
  KE[1,] = kE
  
  SE = matrix(0,nrow = m,ncol = n)
  SE[1,] = sE
  
  NINE = matrix(0,nrow = m,ncol = n)
  NINE[1,] = ninE
  
  
  
  for (t in 2:m){   ## loop over time
    
    
    for (i in 1:n){  ## loop over nodes
      
      if (sE[i] == 0 && kE[i] == 0){
        probOfNoInfection = exp(-beta*sum(AE[i,]*ifelse(sE==1,1,0)))  # "ifelse" is not needed for SIS; but it would be needed for SIR
        sE[i] = sample(x = 0:1, replace = FALSE, 
                       size=1,
                       prob = c(probOfNoInfection,1-probOfNoInfection))
        
        probOfNotContacted = exp(-rao*sum(AE[i,]))
        kE[i] = sample(x = 0:1,
                       size=1,replace = FALSE, 
                       prob = c(probOfNotContacted,1-probOfNotContacted))
        
      }
      
      else    if (sE[i] == 1 && kE[i] == 0 ){
        probOfNoRecovery = exp(-gamma)
        sE[i] = sample(x = 0:1, replace = FALSE, 
                       size = 1,
                       prob = c(1-probOfNoRecovery,probOfNoRecovery))
        
        probOfNotContacted = exp(-rao*sum(AE[i,]))
        kE[i] = sample(x = 0:1,
                       size=1,
                       prob = c(probOfNotContacted,1-probOfNotContacted))
        
      }
      
      else if (sE[i] == 0 && kE[i] == 1){
        probOfNoInfection = exp(-beta*sum(AE[i,]*ifelse(sE==1,1,0)))  # "ifelse" is not needed for SIS; but it would be needed for SIR
        sE[i] = sample(x = 0:1, replace = FALSE, 
                       size=1,
                       prob = c(probOfNoInfection,1-probOfNoInfection))
        
      }
      
      else    if (sE[i] == 1 && kE[i] == 1 ){
        probOfNoRecovery = exp(-gamma)
        sE[i] = sample(x = 0:1, replace = FALSE, 
                       size = 1,
                       prob = c(1-probOfNoRecovery,probOfNoRecovery))
        
      }
      
      
    }  ## i loop 
    
    
    SE[t,] = sE
    KE[t,] = kE
    
    
  }  ## t loop
  
  
  
  
  
  order(colSums(KE)) ## The order of the individuals based on the discovery date
  KsE=subset(KE,select = order(colSums(KE))) ## Sort K based on the order above
  SsE=subset(SE,select = order(colSums(KE))) ## Sort S based on the order above
  Ks2E <- replace(KsE, KsE == 0, -1) ## replace any unknown(0) with (-1)
  final2E <- replace(SsE, Ks2E == -1, -1) ## replace the infectious status for any unknown individual in Ss with (-1)
  
  
  
  niE = length(which(final2E==1)) ## number of infected individuals in final2E 
  nsE = length(which(final2E==0)) ## number of susceptible individuals in final2E 
  totE = niE+nsE ## total number of individuals with known status
  prevE = niE / totE ## the percentage of infected individuals 
  final2PE <- replace(final2E, final2E == -1, prevE)  ## replace the unknown values with prevelance
  
  
  ## calculate gamma
  counterE = 0 ## A counter of how many times an individual remains I given it was I.
  
  for(t in 3:m){
    
    for(i in 1:n){
      if (final2E[t-1,i]== 1 && final2E[t,i]== 1){
        counterE = counterE + 1
      }
      
      
    } # i loop
  } # t loop  
  
  
  gammaHatE = log(niE/counterE)
  listGammaHatE[[b]] <- gammaHatE
  
  
  
  ### start the backward tracing method   
  
  sampleWithoutSurprisesE <- function(listOfUISE) {
    if (length(listOfUISE) <= 1) {
      return(listOfUISE)
    } else {
      return(sample(listOfUISE,1))
    }
    
  } # function
  
  for(t in (m-1):1){
    prob0 = 0
    indexE = 0
    listOfUISE <- list()
    
    for(i in 1:n){
      if(final2PE[t,i] == prevE){
        listOfUISE <- c(listOfUISE,i)
        indexE <- indexE + 1
      }
    } # i loop 
    
    if(length(listOfUISE) > 0) {
      
      
      for (w in 1:indexE){
        
        # pick a random sample from the listOfUIS
        
        delE<- as.integer(sampleWithoutSurprisesE(listOfUISE)[1])
        
        
        # remove the individual from the list of listOfUIS
        listOfUISE <-listOfUISE[-c(match(delE,listOfUISE))]
        
        
        MME =  final2PE[t,] %*% AE[delE,] 
        
        # Update the probability weights
        
        
        if (final2PE [t+1, delE] == 1){
          probI = (exp(-gamma)* prevE)/((exp(-gamma) *prevE) + (1-exp(-beta*sum(MME)) * (1-prevE)))
          
          
          final2PE[t,delE] = sample(x = 0:1,
                                    size=1,replace = FALSE, 
                                    prob = c(probI,1-probI))
        }
        
        
        else if ( final2PE[t+1,delE]== 0){
          probS = ((1-exp(-gamma))* prevE)/(((1-exp(-gamma)) *prevE) + (exp(-beta*sum(MME)) * (1-prevE)))
          
          
          final2PE[t,delE] = sample(x = 0:1,
                                    size=1,replace = FALSE, 
                                    prob = c(probS,1-probS))
          
        }
        
      }  # w loop 
    } ## if > 0 
  } #t loop 
  
  
  loglikbetaE = function(beta){
    sumlog = 0
    prob = 0
    
    for (t in 2:m){   ## loop over time
      
      for (i in 1:n){  ## loop over nodes
        
        ninE =  final2PE[t,]%*%AE
        NINE[t,] = ninE[1,]
        
        if (final2PE[t-1,i]== 0 && final2PE[t,i]== 0){
          prob = exp(-beta*ninE[,i])
          sumlog = sumlog + log(prob)
        }
        
        else if (final2PE[t-1,i]== 0 && final2PE[t,i]== 1){
          prob = (1- exp(-beta*ninE[,i]))
          sumlog = sumlog + log(prob)
        }
      } ## i loop 
    } # t loop 
    
    return(sumlog)
  } ## The loglikbetaE function
  
  oBE = optimize(loglikbetaE, c(0,1), maximum = TRUE)[1]
  listBetaHatE[[b]] <- oBE
  
} ## v loop 


plot(unlist(listBetaHatE),  type="o", col="blue", pch="o", lty=1, xlab = "Iteration", ylab = "Beta_hat", ylim = c(0.009, 0.050) )
legend("bottomright", legend =c("Erdos-Renyi"), col=c("blue"), lty =1, bty="n")
abline(h= c(0,0.01))




end_time<- Sys.time()
Total_Time = end_time - start_time


outliersBetaHatE <- boxplot(unlist(listBetaHatE), plot=FALSE)$out
OBHE<-unlist(listBetaHatE)
OBHE<- OBHE[-which(OBHE %in% outliersBetaHatE)]
boxplot(unlist(OBHE), main = expression(paste(hat(beta)," for ER")))


outliersGammaHatE <- boxplot(unlist(listGammaHatE), plot=FALSE)$out
OGHE<-unlist(listGammaHatE)
OGHE<- OGHE[-which(OGHE %in% outliersGammaHatE)]
boxplot(unlist(OGHE), main = expression(paste(hat(gamma)," for ER")))
