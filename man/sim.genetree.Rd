\name{sim.genetree}
\alias{sim.genetree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{sim.genetree: Simulating coalescent gene trees.}
\description{
sim.genetree simulates a gene tree assuming the coalescent with coalescent rate being 1. The method returns summary statistics of the gene tree.
}
\usage{
sim.genetree(n, numbsim)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{Number of extant sampled tips.}
  \item{numbsim}{Number of trees to simulate.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\item{statistics }{For each simulated gene tree the following statistics are returned: "Colless","s","Sackin","cherries".}
}

\references{T. Stadler, J. Degnan, N. Rosenberg. Manuscript.}
\author{Tanja Stadler }

\examples{
n<-10
numbsim<-2

sim.genetree(n, numbsim)

}