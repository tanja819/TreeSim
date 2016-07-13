\name{corsim}
\alias{corsim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{corsim: Simulating the missing speciation events in an incomplete phylogenies.}
\description{
corsim simulates the missing speciation event in an incomplete phylogeny assuming a constant speciation and extinction rate. These rates can be estimated with the functions bd.shifts.optim (if random speciation events are missing) and bd.groups.optim (if only young speciation events are missing) provided in TreePar. corsim allows to specify an upper and lower bound for the times of the missing speciation events.
}
\usage{
corsim(x,lambda,mu,missing,told=0,tyoung=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Vector of the speciation times in the incomplete phylogeny (where time is measured such that 0 is the present and increasing going into the past).}
  \item{lambda}{Speciation rate.}
  \item{mu}{Extincion rate.}
  \item{missing}{Number of missing species (i.e. missing speciation events).}
  \item{told}{Upper bound for the time of missing speciation events. Default is 0 which means no upper bound.}
  \item{tyoung}{Lower bound for the time of missing speciation events. Default is 0 which means no lower bound. tyoung<told unless tyoung=told=0 is required.}
}
\value{
\item{x}{Vector of speciation times: input and simulated speciation times.}
}
\references{N. Cusimano, T. Stadler, S. Renner. A new method for handling missing species in diversification analysis applicable to randomly or non-randomly sampled phylogenies. Syst. Biol., 61(5): 785-792, 2012.
}
\author{Tanja Stadler }
\seealso{ sim.bd.age, sim.rateshift.taxa, sim.gsa.taxa, birthdeath.tree }
\examples{
# Speciation times of a tree with five species:
x<-c(1,1.5,3,5)
# We simulate using the following parameters:
lambda<-2
mu<-1
tyoung<-0.5
told<-4.5
# We simulate 5 additional speciation times (i.e. five additional species):
missing<-5

# xcompleted is x plus 5 additional speciation events between 0.5 and 4.5 timesteps 
# in the past. xcompleted corresponds to a 10-species tree:
xcompleted<-corsim(x,lambda,mu,missing,told,tyoung)
}