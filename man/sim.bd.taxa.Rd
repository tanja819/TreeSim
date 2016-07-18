\name{sim.bd.taxa}
\alias{sim.bd.taxa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{sim.bd.taxa: Simulating birth-death trees on a fixed number of extant taxa.}
\description{
sim.bd.taxa simulates trees on n species under the constant rate birth-death process. The method allows for incomplete sampling, i.e. (i) only a fixed fraction of all extant tips is included in the sampled tree or (ii) each extant tip from a complete tree is included with a fixed probability. In both cases, the tree is conditioned to have n tips after sampling. If you want to relax the assumption of constant rates, this function will not work. If you want to change rates through time use sim.rateshift.taxa. If you want to have species-age dependent rates, use sim.taxa in R package TreeSimGM.
}
\usage{
sim.bd.taxa(n, numbsim, lambda, mu, frac = 1, complete = TRUE,
 stochsampling = FALSE, fast = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{Number of extant sampled tips.}
  \item{numbsim}{Number of trees to simulate.}
  \item{lambda}{Speciation rate.}
  \item{mu}{Extinction rate.}
  \item{frac}{When complete = FALSE and stochsampling=FALSE: The actual number of tips is n/frac, but only n tips are included (incomplete sampling). When complete = FALSE and stochsampling=TRUE: Each tip is included into the final tree with probability frac. When complete = TRUE: all extinct and non-sampled lineages are included, i.e. the tree has n/frac extant tips.}
  \item{complete}{If TRUE, the tree with the extinct and non-sampled lineages is returned. If FALSE, the extinct lineages are suppressed.}
  \item{stochsampling}{See frac.}
  \item{fast}{Use a faster version of the simulation that takes advantage of compiled Fortran code. See fast.tree for additional functionality using this method.}
}
%\details{
%todo}
\value{
\item{out}{List of numbsim simulated trees with n extant sampled tips.}
}
\references{ T. Stadler. Simulating trees on a fixed number of extant species. Syst. Biol., 60: 676-684, 2011.

T. Stadler. On incomplete sampling under birth-death models and connections to the sampling-based coalescent. Jour. Theo. Biol. 261: 58-66, 2009.}
\author{Tanja Stadler, Nick Beeton }
\note{ 
For stochsampling = TRUE: The algorithm is fast for the critical process, lambda=mu.
}
\seealso{ sim.bd.age, sim.rateshift.taxa, sim.gsa.taxa, birthdeath.tree, sim.taxa, fast.tree }
\examples{
n<-10
lambda <- 2.0
mu <- 0.5
frac <-0.6
numbsim<-2

##
# Simulating numbsim trees with n species under a birth-death process with 
# speciation rate lambda an extinction rate mu:

sim.bd.taxa(n,numbsim,lambda,mu)

# Each extant species is included in final tree with probability frac 
# (the tree has n species AFTER sampling):

sim.bd.taxa(n,numbsim,lambda,mu,frac,complete=FALSE,stochsampling=TRUE)

# A fraction frac of the extant species is included into the final tree 
# (the tree has n species AFTER sampling):

sim.bd.taxa(n,numbsim,lambda,mu,frac,complete=FALSE,stochsampling=FALSE)

## 
# Demonstration of the two tree techniques
#

# use microbenchmark library to compare run times for a single tree
# (note that improvement decreases for multiple trees using sim.bd.taxa due to overheads)
library(microbenchmark)

# generate a single 1000 extant leaf random tree using raw fast.tree function
microbenchmark(fast.tree(n = 1000, lambda = 1.5, mu = 1, frac = 1), times = 1) 
# generate the same tree using the sim.bd.taxa framework (which uses fast.tree)
microbenchmark(sim.bd.taxa(n = 1000, numbsim = 1, lambda = 1.5, mu = 1, frac = 1, fast = TRUE), times = 1)
# now generate a tree using the original technique
microbenchmark(sim.bd.taxa(n = 1000, numbsim = 1, lambda = 1.5, mu = 1, frac = 1, fast = FALSE), times = 1)

##
# suite of tests to compare both tree generating techniques

# generate 10,000 random 100-leaf trees using fast and original techniques
# This will generally take a few minutes to run
N = 10000
test1 = sim.bd.taxa(100, N, 1.5, 1, fast = FALSE)
test2 = sim.bd.taxa(100, N, 1.5, 1, fast = TRUE)

# Examine number of nodes - histograms should look similar
nnode1 = unlist(lapply(test1, Nnode))
nnode2 = unlist(lapply(test2, Nnode))
hist(nnode1, breaks = 10*(10:80), col='#0000FF80', ylim = c(0,N/10), main = 'Number of nodes')
hist(nnode2, breaks = 10*(10:80), col='#FF000080', add = TRUE)

# Examine overall tree length
tl1 = unlist(lapply(test1, function(x) max(node.depth.edgelength(x))))
tl2 = unlist(lapply(test2, function(x) max(node.depth.edgelength(x))))
hist(tl1, breaks = 0.5*(0:60), col='#0000FF80', ylim = c(0,N/10), main = 'Tree length (time)')
hist(tl2, breaks = 0.5*(0:60), col='#FF000080', add = TRUE)

# Examine mean node depth
nd1 = unlist(lapply(test1, function(x) mean(node.depth(x))))
nd2 = unlist(lapply(test2, function(x) mean(node.depth(x))))
hist(nd1, breaks = 0.5*(0:60), col='#0000FF80', ylim = c(0,N/10), main = 'Mean node depth (time)')
hist(nd2, breaks = 0.5*(0:60), col='#FF000080', add = TRUE)

}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{birth-death tree}
%\keyword{ phylogenetic tree}