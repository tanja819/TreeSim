\name{sim.bd.taxa}
\alias{sim.bd.taxa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{sim.bd.taxa: Simulating birth-death trees on a fixed number of extant taxa.}
\description{
sim.bd.taxa simulates trees on n species under the constant rate birth-death process. The method allows for incomplete sampling, i.e. (i) only a fixed fraction of all extant tips is included in the sampled tree or (ii) each extant tip from a complete tree is included with a fixed probability. In both cases, the tree is conditioned to have n tips after sampling. If you want to relax the assumption of constant rates, this function will not work. If you want to change rates through time use sim.rateshift.taxa. If you want to have species-age dependent rates, use sim.taxa in R package TreeSimGM.
}
\usage{
sim.bd.taxa(n, numbsim, lambda, mu, frac = 1, complete = TRUE,
 stochsampling = FALSE)
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
}
%\details{
%todo}
\value{
\item{out}{List of numbsim simulated trees with n extant sampled tips.}
}
\references{ T. Stadler. Simulating trees on a fixed number of extant species. Syst. Biol., 60: 676-684, 2011.

T. Stadler. On incomplete sampling under birth-death models and connections to the sampling-based coalescent. Jour. Theo. Biol. 261: 58-66, 2009.}
\author{Tanja Stadler }
\note{ 
For stochsampling = TRUE: The algorithm is fast for the critical process, lambda=mu.
}
\seealso{ sim.bd.age, sim.rateshift.taxa, sim.gsa.taxa, birthdeath.tree, sim.taxa }
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

}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{birth-death tree}
%\keyword{ phylogenetic tree}