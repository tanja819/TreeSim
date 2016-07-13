\name{sim.rateshift.taxa}
\alias{sim.rateshift.taxa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{sim.rateshift.taxa: Simulating trees incorporating mass extinction events and rate shifts.}
\description{
sim.rateshift.taxa simulates trees on n species under the constant rate birth-death process. At user-specified points in the past, the rates can shift. Further, mass extinction events can be incorporated. The method further allows for incomplete sampling, i.e. only a fixed fraction of all tips is included in the sampled tree. The tree is conditioned to have n tips after sampling.
}
\usage{
sim.rateshift.taxa(n, numbsim, lambda, mu, frac, times, complete = TRUE, K=0, norm = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{Number of extant sampled tips.}
  \item{numbsim}{Number of trees to simulate.}
  \item{lambda}{Vector of speciation rates, the rate in entry i is the speciation rate prior (ancestral) to time times[i].}
  \item{mu}{Vector of extinction rates, the rate in entry i is the extinction rate prior (ancestral) to time times[i].}
  \item{frac}{Vector of proportion of species surviving mass extinction event. Entry i corresponds to the mass extinction at time times[i]. If frac[i]=1, only rate shift but no mass extinction at time times[i].
}
  \item{times}{Vector of mass extinction and rate shift times. Time is 0 today and increasing going backwards in time. Specify the vector as times[i]<times[i+1]. times[1]=0 (today).}
  \item{complete}{If TRUE, the tree including the extinct lineages and non-sampled lineages is returned (so the tree has round(n/frac[1]) extant tips). If FALSE, the extinct lineages and non-sampled lineages are suppressed.}
\item{K}{If K>0, then a density-dependent speciation rate = lambda*(1-numberspecies/K) is used. Only works currently for mu=0.}
\item{norm}{If norm = TRUE the simulations are exact. If norm = FALSE tree is always returned once N=0 in Stadler 2011, p.678, point (7).}
}
%\details{
%todo
%}
\value{
  \item{out}{List of numbsim simulated trees with n extant sampled tips.}
}
\references{
T. Stadler. Simulating trees on a fixed number of extant species. Syst. Biol., 60: 676-684, 2011.
}
\author{Tanja Stadler }
%\note{ %~~further notes~~ 
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{
sim.bd.age, sim.bd.taxa, sim.gsa.taxa, birthdeath.tree
}
\examples{
n<-10
numbsim<-1

##
# Simulating trees with a fixed number of species having shifts in rate
# and mass extinction events.
# Between today and time 0.3 in the past, we have speciation rate 2,
# extinction rate 0. At time 0.3, we have a mass extinction event which
# 10% of the species survive. Prior to 0.3, we have a speciation rate 
# of 1 and an extinction rate of 0.3:

sim.rateshift.taxa(n,numbsim,c(2,1),c(0,0.3),
c(1,0.1),c(0,0.3),complete=TRUE)

# The fraction 0.6 of the extant species is included into the final tree
# (the tree has n species AFTER sampling, extinct and
# non-sampled lineages are not included):

sim.rateshift.taxa(n,numbsim,c(2,1),c(0,0.3),
c(0.6,0.1),c(0,0.3),complete=FALSE)
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ mass extinction }
%\keyword{ simulating trees }
%\keyword{ phylogenetic tree }
%\keyword{rateshift}