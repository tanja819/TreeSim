\name{LTT.plot.gen}
\alias{LTT.plot.gen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{LTT.plot.gen: Calculates the number of lineages through time for each input tree, as well as the average number of lineages over all trees.}
\description{
LTT.plot.gen calculates the number of lineages through time for each input tree, as well as the average number of lineages over all trees. The trees may be simulated using any function in TreeSim, or may be empirical trees. The method works for ultrametric and non-ultrametric trees which are binary or have polytomies.
}
\usage{
LTT.plot.gen(trees,bound=10^(-12))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trees}{List with one or two entries. First component: list of phylogenetic trees; second component: vector with time of origins (can be empty).}
  \item{bound}{Determines the value by which leaf times may differ in an ultrametric tree. If two tips are further apart, they are considered as sequentially sampled tips.}
}
%\details{
%todo}
\value{
	\item{out}{out[[1]]: First column are the branching times of ALL input trees. Second column is the number of lineages after the branching time. out[[i]]: Equivalent vector as out [[1]], but for tree i-1.}
}
\references{ T. Stadler. Simulating trees on a fixed number of extant species. Syst. Biol., 60: 676-684, 2011.
}
\author{Tanja Stadler }
\seealso{ LTT.plot,sim.bd.taxa, sim.bd.age, sim.rateshift.taxa, sim.gsa.taxa, birthdeath.tree }
\examples{
# Simulation of a tree with a mrca at time 10 in the past, 
# under the density-dependent model
numbsim<-10
age<-10
lambda<-0.3
mu<-0.2
K<-40
# You can produce LTT plots as follows. 
# (for now this is un-commented, as some combinations of geiger / TreeSim on certain platforms 
# produce problems. If this is the case for you, please report to tanjs.stadler@bsse.ethz.ch).
% tree<- sim.bd.age(age,numbsim,lambda,mu,mrca=TRUE,complete=FALSE,K=K)
% # Plot of tree
% ltttrees<-LTT.plot.gen(list(tree))
% plot(ltttrees[[1]],type="l")
% #alternatively use LTT.plot() for plotting:
% LTT.plot(list(tree),avg=TRUE)

% # Simulation of a tree with the first lineage at time 10 in the past, 
% # under the density-dependent model
% numbsim<-10
% age<-10
% lambda<-0.3
% mu<-0
% K<-40
% tree<- sim.bd.age(age,numbsim,lambda,mu,mrca=FALSE,complete=FALSE,K=K)
% # Plot of tree
% ltttrees<-LTT.plot.gen(c(list(tree),list((1:10)*0+age)))
% plot(ltttrees[[5]],type="l")
% #alternatively use LTT.plot() for plotting!
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{birth-death tree}
%\keyword{ phylogenetic tree}
