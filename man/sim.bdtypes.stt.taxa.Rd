\name{sim.bdtypes.stt.taxa}
\alias{sim.bdtypes.stt.taxa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{sim.bdtypes.stt.taxa: Simulating multitype birth-death trees with a fixed number of tips sampled through time.}
\description{
sim.bdtypes.stt.taxa simulates trees on n tips sampled through time under a multitype birth-death process.
}
\usage{
sim.bdtypes.stt.taxa(n,lambdavector,deathvector,
sampprobvector,init=-1,EI=FALSE,eliminate=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{Number of sampled tips.}
  \item{lambdavector}{Matrix of dimension kxk, where k is the number of different states. The entry (i,j) is the rate with which an individual in state i gives rise to a new lineage of state j.}
  \item{deathvector}{Vector of dimension k, the entry i is the death rate of an individual in state i.}
  \item{sampprobvector}{Vector of dimension k, the entry i is the probability of an individual in state i being sampled upon death, i.e. being included into the final tree.}
  \item{init}{Default is -1, meaning the initial individual is in a random state (which is chosen from the equilibrium distribution of states). If init>0, then the initial state is 'init'.}
  \item{EI}{If EI=TRUE a model with two types, namely exposed and infectious individuals, is assumed. Infectious individuals transmit and give rise to exposed individuals with rate lambdavector[2,1], and exposed individuals become infectious with rate lambdavector[1,2]. Exposed individuals have a 0 death rate and cannot be sampled. For an example simulation see below.}
  \item{eliminate}{Only relevant if EI=TRUE. Under EI=TRUE all sampled tips are in state 2. If eliminate>0, the first eliminate tips are marked with state 1. This facilitates further analysis, e.g. we now can easily prune these first eliminate tips to mimic no sampling close to the epidemic outbreak.}
}
%\details{
%todo}
\value{
\item{out}{Phylogenetic tree with 'n' sampled tips. In out$states, the states for the tips are stored.}
}
\references{ T. Stadler, S. Bonhoeffer. Uncovering epidemiological dynamics in heterogeneous host populations using phylogenetic methods. Phil. Trans. Roy. Soc. B, 368 (1614): 20120198, 2013.
}
\author{Tanja Stadler}
\note{ 
A large number of trees can be obtained using the R function lapply. The tree can be plotted using the R package ape function plot(out). sim.bdtypes.stt.taxa function extends the simulator in the R package diversitree to trees which contain tips being sampled sequentially.}
\seealso{ sim.bdsky.stt }
\examples{
# Simulate two trees with 10 tips
set.seed(1)
n<-10
lambda <- rbind(c(2,1),c(3,4))
mu <- c(1,1)
sampprob <-c(0.5,0.5)
numbsim<-2
trees<-lapply(rep(n,numbsim),sim.bdtypes.stt.taxa,
lambdavector=lambda,deathvector=mu,sampprobvector=sampprob)

# Testing the model with exposed class (EI = TRUE)
set.seed(2)
# simulate tree with expected incubation period of 14 days, 
# infectious period of 7 days, and R0 of 1.5:
mu <- c(0,1/7)
lambda <- rbind(c(0,1/14),c(1.5/7,0))
# sampling probability of infectious individuals is 0.35:
sampprob <-c(0,0.35)
# we stop once we have 20 samples:
n <- 20
# we simulate one tree:
numbsim<-1
# We mark first eliminate=10 tips such that we can easily drop them later
# (if deleting these 10 tips, we mimic no sampling close to the outbreak)
trees<-lapply(rep(n,numbsim),sim.bdtypes.stt.taxa,lambdavector=lambda,deathvector=mu,
sampprobvector=sampprob,EI=TRUE,eliminate=10)
}