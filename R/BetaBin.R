betabin <- function(k,n,alpha,beta)  return(exp(lgamma(n+1)-lgamma(k+1)-lgamma(n-k+1) +lgamma(k+alpha)+lgamma(n-k+beta)-lgamma(n+alpha+beta)+lgamma(alpha+beta)-lgamma(alpha)-lgamma(beta)))

meanbeta <- function(n,alpha,beta) return(n*alpha/(alpha+beta))
varbeta  <- function(n,alpha,beta) return(n*alpha*(n*(1+alpha)+beta)/((alpha+beta)*(1+alpha+beta)) -  (n*alpha/(alpha+beta))^2)


vardivmu <- function(n,alpha,beta)  return(varbeta(n,alpha,beta)/meanbeta(n,alpha,beta))

beta <- 2
n <- 5
getalpha <- function(meanvalue,beta,n)
  return(beta*meanvalue/(n-meanvalue))


# Need average growth of 5.
alpha <- getalpha(1.5,beta,n)
# print(vardivmu(n,alpha,beta)) # 4.11

TEST <- FALSE
if(TEST) {
  probs <- betabin(0:n,n,alpha,beta)
  sum(probs)
  sum(0:5*probs)

  beta <- 99
  alpha <- getalpha(1.5,beta,n)
  print(varbeta(n,alpha,beta)) #1.7
  probs <- betabin(0:n,n,alpha,beta)
  print(sum((0:5)^2*probs)-sum(0:5*probs)^2) # should get the same value

#Dæmi um notkun (Hærra gilti á beta gefur lægri meiri varíanse fyrir gefið meðaltal.
# Vil fá meðalvöxt upp á 1.5 lengdarbil, max möstur á tímaeiningu 5 lengdarbil.
  beta <- 2
  alpha <- getalpha(1.5,beta,n)
  probs <- betabin(0:n,n,alpha,beta)

# probs eru þá vaxtarmatrixan fyrir 0, 1, 2, 3, 4 og 5 lengarbil.
# má auðvitað ekki reyna eitthvað sem gengur ekki þ.e meðalvöxtur < n
# dæmi þar sem meðalvöxtur fer mjög nærri n
# meðalvöxtur 2.8 lengdarbil
  beta <- 2
  n <- 3
  alpha <- getalpha(2.8,beta,n)
  probs <- betabin(0:n,n,alpha,beta)
}
