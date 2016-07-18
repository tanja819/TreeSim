\name{fast.tree}
\alias{fast.tree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{fast.tree: Fast version of sim.bd.taxa for a single tree, with Brownian character tracking}
\description{
fast.tree replicates the functionality of sim.bd.taxa for a single tree (more specifically, the internal function sim2.bd.reverse.single). When traits is set to TRUE it will also simulate the evolution of a character over the tree with Brownian motion, returning its result for every node and tip. In addition sample points across the length of the tree can be given and fast.tree will return a time slice of the character values for the extant branches given at these time points.
}
\usage{
fast.tree = function(n, lambda, mu, frac = 1, traits = FALSE, sigma = 1, sampleTimes = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{Number of extant sampled tips.}
  \item{lambda}{Speciation rate.}
  \item{mu}{Extinction rate.}
  \item{frac}{For use with incomplete sampling (see sim.bd.taxa). The actual number of tips is n/frac.}
  \item{traits}{When set to TRUE, enables character tracking and sampling.}
  \item{sigma}{The standard deviation of the Brownian motion of the character over time.}
  \item{sampleTimes}{A vector of times at which to sample, expressed as a proportion of the overall tree length (in time).}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  If traits = FALSE, a list with two values, the first containing the simulated tree with n/frac extant tips, and the second containing the tree length. This format is used for consistency with sim2.bd.reverse.single.
  
  If traits = TRUE, a list with items:
  \item{tree}{The simulated tree with n/frac extant tips.}
  \item{samples}{A max.leaves x n.samples matrix of characters taken at sampleTimes (padded with zeros), where max.leaves is the maximum number of leaves at a sampled time and n.samples is the length of vector sampleTimes.}
  \item{n.leaves}{The number of leaves extant at each of the sampled times.}
  \item{sampleTimes}{The vector of sampled times (no longer expressed as a proportion).}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
  Nick Beeton
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
sim.bd.taxa
}
\examples{
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
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
