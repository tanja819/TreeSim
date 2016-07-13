\name{sim.genespeciestree}
\alias{sim.genespeciestree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{sim.genespeciestree: Simulating birth-death species trees with nested coalescent gene trees.}
\description{
sim.genespeciestree simulates birth-death species trees (using sim.bd.taxa or sim.bd.taxa.age). Within each species tree, a gene tree is simulated, assuming a coalescent with coalescent rate being 1. The method returns summary statistics for the gene tree.
}
\usage{
sim.genespeciestree(n, numbsim, lambda, mu, frac = 1, age=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{Number of extant sampled tips.}
  \item{numbsim}{Number of trees to simulate.}
  \item{lambda}{Speciation rate.}
  \item{mu}{Extinction rate.}
  \item{frac}{ Each tip is included into the final species tree with probability frac.}
  \item{age}{The time since origin / most recent common ancestor. If age = 0 (default) a uniform prior for the time since origin is assumed.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\item{statistics }{For each simulated gene tree the following statistics are returned (with "gammaspecies" being the gamma statistic for the corresponding species tree): "Colless","s","Sackin","cherries","matching of species tree and gene tree","gammaspecies","gamma".}
}

\references{T. Stadler, J. Degnan, N. Rosenberg. Manuscript.}
\author{Tanja Stadler }

\examples{
#Simulate two gene trees within two species trees:
n<-10
lambda <- 2.0
mu <- 0.5
frac <-0.6
numbsim<-2
age<-2

# Simulation is conditioned on 10 final tips
sim.genespeciestree(n, numbsim, lambda, mu, frac, 0)

# Simulation is conditioned on 10 final tips and tree age 2
sim.genespeciestree(n, numbsim, lambda, mu, frac, age)
}