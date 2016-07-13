\name{sim.gsa.taxa}
\alias{sim.gsa.taxa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
sim.gsa.taxa: Sampling trees on n tips from bigger trees.}
\description{
sim.gsa.taxa samples trees on n tips (using the GSA approach, see references) from trees with m tips where m>n, given the m-tip trees are simulated under the simple sampling approach (i.e. simulating until first m>>n tips are reached or the tree is extinct).
The TreeSim methods to simulate n-tip trees, sim.bd.taxa and sim.rateshift.taxa, are implemented such that sim.gsa.taxa is not necessary. sim.gsa.taxa is needed for post processing of trees generated NOT in TreeSim: if the aim is to simulate trees with n co-existing tips, then typically simulators stop once the first time n co-existing lineages are reached. However, due to death, we can observe n tips later (e.g. n+1 lineages followed by death leads n lineages). sim.gsa.taxa produces an appropriate set of n-tip trees where the input are m-tip trees with m>>n and the m-tip trees are simulated under these typical simulators.

 sim.gsa.taxa works for m-tip trees generated under a model where: (i) the number of tips eventually tends to zero or stays bigger than n and (ii) birth / death rate changes do not depend on the time between the change and the present - e.g. one cannot model a mass extinction event 1 million years BEFORE the present. But one can model a mass extinction event 1 million years AFTER the origin of the tree. The package TreeSimGM uses sim.gsa.taxa to obtain n-tip trees.
}
\usage{
sim.gsa.taxa(treearray, n, frac = 1, sampling = 1, complete = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{treearray}{Array of trees with a fixed number of tips.}
  \item{n}{Number of tips in sampled trees.}
  \item{frac}{Relevant when complete = FALSE: The actual number of tips is n/frac, but only n tips are included (incomplete sampling).
When complete = TRUE: We set frac = 1.}
  \item{sampling}{Parameter determining how close the returned trees in treearray are to the "true" distribution. The higher 'sampling', the closer the output trees to the 'true' distribution. Default is 40. Higher values of sampling return fewer output trees meaning a larger input treearray is needed.}
  \item{complete}{If TRUE, the tree with the extinct lineages is returned. If FALSE, the extinct lineages are suppressed.}
}
%\details{
  %~~ If necessary, more details than the description above ~~
%}
\value{
\item{treearray}{Array of sampled trees with n  extant sampled tips. Note that the number of trees in the output is significantly smaller than the number of trees in the input (in order to ensure correct tree sampling).}
}
\references{K. Hartmann, D. Wong, T. Stadler. Sampling trees from evolutionary models. Syst. Biol., 59(4): 465-476, 2010.

T. Stadler. Simulating trees on a fixed number of extant species. Syst. Biol., 60: 676-684, 2011.}
\author{Tanja Stadler}
%\note{ %~~further notes~~ 
%~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{sim.bd.age, sim.bd.taxa, sim.rateshift.taxa, birthdeath.tree
} %~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##
# First 100 trees on 9 tips under a birth-death process are generated.
# Then trees on 5 species are sampled from these 100 trees using the GSA
# (see references). 
# You can easily simulate trees on m species with the simple sampling
# approach (see references) under a variety of models. Then use the
# provided GSA algorithm to get a correct sample of trees on n<<m species:

m<-9
n<-5
numbsim<-100
lambda <- 2.0
mu <- 0.5

t<-sim.bd.taxa(m,numbsim,lambda,mu)
t2<-sim.gsa.taxa(t,n)
}

%\keyword{tree}
%\keyword{simulation}
