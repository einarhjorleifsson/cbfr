
library(Rcpp)
source("R/BetaBin.R")
cppFunction('NumericVector GrowthTransition(NumericMatrix Growthprob,NumericVector Nold) {
  int Maxgr = Growthprob.ncol(), nlgr = Growthprob.nrow();
  int le,gr;
  NumericVector Nnew(nlgr);
  for(le=0; le < nlgr; le++) Nnew(le) = 0; // Not needed
  for(le=0; le < nlgr-Maxgr+1; le++) {
    for(gr = 0;gr < Maxgr ; gr++){
      Nnew(le+gr) += Nold(le)*Growthprob(le,gr);
    }
  }
  for(le=nlgr-Maxgr+1; le < nlgr ; le++) {
    for(gr = 0;gr < Maxgr ; gr++){
      if((le+gr) < nlgr)
        Nnew(le+gr) += Nold(le)*Growthprob(le,gr);
      else
        Nnew(nlgr-1)  += Nold(le)*Growthprob(le,gr); // Maximum growth reached.
    }
  }
  return Nnew;
}')



AddRecruitment <- function(MeanLrec,sdLrec,Nrec){
  maxl <- trunc(meanLrec+3*sdLrec+1) #could use if we do not want too large fish
  N1 <- diff(c(0,pnorm(Lengths,MeanLrec,sdLrec)))*Nrec
  return(N1)
}

#Natm will only be called once but at every timestep if it is related
#density
CalcNatM <- function(Length,Mr,Lr,Rho) #
  return(Mr*(Lr/Length)^(3*Rho)) # use Lr/L instead of L/Lr to get rid of power-

simulateGrowthOneTimestep <- function(N,i) { # i is index to time
# Weights Linf,K,d,dt,maxgr,beta, growthprob are global variables.
# growthprob is changed but the global variable is always 0.
  Biom <- sum(N*Weights) # Per hectar
  LinfB <- Linf - d*Biom # Density corrected Linf
  LinfDensDep[i] <<- LinfB  # Store densdep growth
  Growth <- K*(LinfB-Lengths)*dt/dl # Growth has to be in the right units dl used in case something else than cm was used.
  Growth[Growth < 0] <- 0.0001 # Growth can not be 0, lgamma fails
  for(i in 1:length(Growth)) {
    alpha <- getalpha(Growth[i],beta,maxgr)
    growthprob[,i] <<- betabin(0:maxgr,maxgr,alpha,beta)
  }

  N1 <- GrowthTransition(t(growthprob),N)
  return(N1)
}

CalcSel <- function(lengths,params) # Params[1] selslope params[2] l50
  return(1/(1+exp(-params[1]*(lengths-params[2]))))


Lengths <- seq(0,75,by=1) # CM
dl <- mean(diff(Lengths))
wtcoeffs <- c(2.5e-5,3) # Length weight relationship.
Weights <- wtcoeffs[1]*Lengths^wtcoeffs[2]
Time <- seq(0,6,by=1/12)
dt <- 1/12
nsubsteps <- 10  #addition for very high F and M
Year <- trunc(Time)



N <- Catch <- NatMnumbers <- matrix(0,length(Time),length(Lengths))
tmpCatch <- tmpNatMnumbers <- rep(0,length(Lengths)) # Temp storage for substeps

meanLrec <- 5 # Seed (rec) mean length
sdLrec <- 1 # Seed sd
Nrec <- 600 # 1000 recruits/Ha
Linf <- 80 # Base Linf
LinfDensDep <- rep(0,length(Time)) # to store
K <- 0.4 # Von bertalanfy's
# DensDep growth
d <- 0.5 #cm ha/kg


#NatmPar per year
Mr <- 0.4
Lr <- 20
Rho <- 0.5 # reduction with size

beta <- 10  # parameter of betabinomial distribution higher values less dispersion
maxgr <- 3  # Maximum growth per time interval in number of length groups

# Selection of fisheries

Sell50 <- 25 # L50
Selslope <- 3 # Steepness
Effort <- 1

# Growth update  probabilities.
growthprob <- matrix(0,(maxgr+1),length(Lengths))




Selpath <- CalcSel(Lengths,c(Selslope,Sell50))
# Do only calculate Natm once as it is not density dependent.
Natm <- CalcNatM(Lengths,Mr,Lr,Rho)
for(i  in 1:length(Time)){
#  if(((Time[i]) %% 1) == 0) # Once a year rec
#    N[i,] <- N[i,]+AddRecruitment(meanLrec,sdLrec,Nrec)
  N[i,] <- N[i,]+AddRecruitment(meanLrec,sdLrec,Nrec/12) # Every month
  for(j in 1:(nsubsteps)) {
    tmpCatch  <-  N[i,]*(1-exp(-dt/nsubsteps*Effort*Selpath))
    tmpNatMnumbers <-  N[i,]*(1-exp(-dt/nsubsteps*Natm))
    N[i,] <- N[i,] - tmpNatMnumbers - tmpCatch
    Catch[i,] <- Catch[i,] + tmpCatch
    NatMnumbers[i,] <- NatMnumbers[i,] + tmpNatMnumbers
  }
  if(i < length(Time)) N[i+1,] <- simulateGrowthOneTimestep(N[i,],i)

  # Latter part of substep  if we use growth in middle of timestep
  # in that cast nsubsteps would have to be even.
  #for(j in (nsubsteps/2+1):nsubsteps) {
  #  tmpCatch  <- N[i,]*(1-exp(-dt/nsubsteps*Effort*Selpath))
  #  tmpNatMnumbers <-  N[i,]*(1-exp(-dt/nsubsteps*Natm))
  #  N[i,] <- N[i,] - tmpNatMnumbers - tmpCatch
  #  Catch[i,] <- Catch[i,] + tmpCatch
  #  NatMnumbers[i,] <- NatMnumbers[i,] + tmpNatMnumbers
  #}
}


# Results (by time and length stored in the matrices N,NatMnumbers and catch.
# The vectors Time, Years and Weights are also usefule.

print(tapply(apply(Catch,1,sum),Year,sum))

i <- Year==5
cbylen <- apply(Catch[i,],2,sum)
plot(Lengths,cbylen,type="l")

