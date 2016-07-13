\name{cuttree}
\alias{cuttree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
cuttree: Cutting off the tree to prune recent branches.
}
\description{
cuttree takes as input a tree and a cuttime, and then prunes all lineages more recent than cuttime.
}
\usage{
cuttree(tree,cuttime)
}
\arguments{
  \item{tree}{
Phylogenetic tree.
}
\item{cuttime}{Time before present at which all descendent lineages are deleted. Value between 0 (the present i.e. nothing is deleted) and the age of the tree (i.e. the whole tree is deleted).}
}
\value{
\item{tree}{Tree where all branches more recent than cuttime are pruned from the input tree.}
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
n<-10
lambda <- 2.0
mu <- 0.5
frac <-0.6
numbsim<-2
age<-2

##
# Simulating numbsim reconstructed trees with n sampled species under a
# birth-death process with speciation rate lambda, extinction rate mu,
# sampling probability frac, and time age since origin:

trees<-sim.bd.taxa.age(n, numbsim, lambda, mu, frac, age, mrca = FALSE)

treec<-cuttree(trees[[1]],0.3)

plot(trees[[1]])
plot(treec)
}