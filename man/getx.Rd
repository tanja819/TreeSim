\name{getx}
\alias{getx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
getx: Calculating the vector of speciation / transmission times and sampling times for a phylogenetic tree.}
\description{
getx calculates the vector of branching (speciation / transmission) times and sampling times for a phylogenetic tree (which may have polytomies). This vector is the input for the TreePar methods.
}
\usage{
getx(datatree,sersampling)
}
\arguments{
  \item{datatree}{
Phylogenetic tree.
}
\item{sersampling}{Set sersampling = 0 if all tips are from one timepoint; 1 otherwise.}
}
\value{
\item{x}{Vector of branching times where 0 is the present and time increasing into the past; If sersampling = 1: Vector of branching and tip sampling times. Second column indicates for each time if branching event (1) or tip (0).}
}
\author{
Tanja Stadler 
}
%\keyword{diversification rates}
%\keyword{rate shifts}
%\keyword{mass extinction}
%\keyword{maximum likelihood rates}
%\keyword{phylogeny}
\examples{
### tree with tips sampled at one timepoint
n<-10
lambda <- 2.0
mu <- 0.5
frac <-0.6
numbsim<-1
trees<-sim.bd.taxa(n, numbsim, lambda, mu, frac,complete=FALSE,stochsampling=TRUE)
branching<-getx(trees[[1]])

### tree with tips sampled sequentially through time
set.seed(1)
n<-10
lambda <- c(2,1,2)
mu <- c(1,0.5,1.5)
sampprob <-c(0.5,0.5,0.5)
times<-c(0,1,2)
numbsim<-2
trees<-lapply(rep(n,numbsim),sim.bdsky.stt,lambdasky=lambda,deathsky=mu,
timesky=times,sampprobsky=sampprob,rho=0,timestop=0)
branchingserial<-getx(trees[[1]][[1]],sersampling=1)
}