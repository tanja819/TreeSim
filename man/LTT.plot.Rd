\name{LTT.plot}
\alias{LTT.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{LTT.plot: Plots the lineages through time of a set of phylogenetic trees.}
\description{
LTT.plot plots the lineages through time (LTT) for a set of phylogenetic trees in black (complete or reconstructed; with or without polytomies) together with the average LTT plot in red. The trees may be simulated using any function in TreeSim, or may be empirical trees. The method works for ultrametric and non-ultrametric trees which are binary or have polytomies. NOTE: you probably need to adapt the code such that the plot is pretty for your particular data (range of axes etc).}
\usage{
LTT.plot(trees,width,precalc,bound=10^(-12),timemax,nmax,avg)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trees}{List with one or two entries. First component: list of phylogenetic trees; second component: vector with time of origins (can be empty).}
  \item{width}{Width of lines in plot.}
  \item{precalc}{Default = 0. If = 1, then parse 'LTT.plot.gen(trees)' instead of 'trees' for the input variable 'trees'.}
  \item{bound}{Determines the value by which leaf times may differ in an ultrametric tree. If two tips are further apart than 'bound', they are considered as sequentially sampled tips.}
  \item{timemax}{Time axis is drawn from present=0 to timemax years in the past.}
  \item{nmax}{Axis with number of species is drawn from 1 to nmax.}
  \item{avg}{Default=FALSE. If true then the average LTT plot of all individual LTT plots is drawn in red.}
}
%\details{
%todo}
\references{ T. Stadler. Simulating trees on a fixed number of extant species. Syst. Biol., 60: 676-684, 2011.
}
\author{Tanja Stadler }
\seealso{LTT.plot.gen, sim.bd.taxa, sim.bd.age, sim.rateshift.taxa, sim.gsa.taxa, birthdeath.tree }
\examples{
# Simulation of a tree of age 10 under the density-dependent model
numbsim<-3
age<-10
lambda<-0.3
mu<-0
K<-40
tree<- sim.bd.age(age,numbsim,lambda,mu,mrca=FALSE,complete=FALSE,K=K)
# Plot of tree
LTT.plot(c(list(tree),list(c(age,age,age))))
#
# Simulation of a tree with 10 tips under the constant rate birth-death model
numbsim<-3
n<-10
lambda<-0.3
mu<-0
tree<- sim.bd.taxa(10,numbsim,lambda,mu,complete=FALSE,stochsampling=TRUE)
# Plot of tree
ages<-c()
for (i in 1:length(tree)){
	ages<-c(ages,tree[[i]]$root.edge+max(getx(tree[[i]])))
}
LTT.plot(c(list(tree),list(ages)),avg=TRUE)
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{birth-death tree}
%\keyword{ phylogenetic tree}