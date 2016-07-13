\name{sim.bd.taxa.age}
\alias{sim.bd.taxa.age}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{sim.bd.taxa.age: Simulating birth-death trees with a given age on a fixed number of extant taxa.}
\description{
sim.bd.taxa.age simulates trees on n species with a (i) fixed time since origin or (ii) fixed time since the most recent common ancestor of the sampled species under the constant rate birth-death process. The method allows for incomplete sampling, i.e. each extant tip from a complete tree is included with a fixed probability. The tree is conditioned to have n tips after sampling and a fixed time since origin or since the most recent common ancestor of the sampled species.
}
\usage{
sim.bd.taxa.age(n, numbsim, lambda, mu, frac = 1, age, mrca = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{Number of extant sampled tips.}
  \item{numbsim}{Number of trees to simulate.}
  \item{lambda}{Speciation rate.}
  \item{mu}{Extinction rate.}
  \item{frac}{ Each tip is included into the final tree with probability frac.}
  \item{age}{The time since origin / most recent common ancestor.}
  \item{mrca}{If mrca = FALSE: The time since the origin of the process. If mrca = TRUE: The time since the most recent common ancestor of the sampled species.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\item{treearray }{Array of numbsim trees with n>1 tips with a given age. The extinct lineages are not included.}
}

\references{T. Stadler: On incomplete sampling under birth-death models and connections to the sampling-based coalescent.
J. Theo. Biol. (2009) 261: 58-66.}
\author{Tanja Stadler }
\note{
The algorithm is fast for the critical process, lambda=mu.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
sim.bd.age, sim.rateshift.taxa, sim.gsa.taxa, birthdeath.tree 
}
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

sim.bd.taxa.age(n, numbsim, lambda, mu, frac, age, mrca = FALSE)

# Simulating numbsim reconstructed trees with n sampled species under a
# birth-death process with speciation rate lambda, extinction rate mu,
# sampling probability frac, and time age since the most recent 
# common ancestor of the extant sampled species:

sim.bd.taxa.age(n, numbsim, lambda, mu, frac, age, mrca = TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{birth-death tree}
%\keyword{ phylogenetic tree}
