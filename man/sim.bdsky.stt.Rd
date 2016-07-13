\name{sim.bdsky.stt}
\alias{sim.bdsky.stt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{sim.bdsky.stt: Simulating sequentially sampled birth-death, SIS, SIR or SIRS trees where birth and death rates are changing through time.}
\description{
sim.bdsky.stt simulates birth-death trees with tips being sampled sequentially. The birth and death rates may change in a piecewise fashion. The birth rates may additionally depend on the number of susceptible individuals in an epidemic, corresponding to epidemiological SIS, SIR or SIRS dynamics. The trees are conditioned on a fixed number of tips or a fixed age.}
\usage{
sim.bdsky.stt(n,lambdasky,deathsky,timesky,sampprobsky,omegasky=rep(0,length(timesky)),
rho=0,timestop=0,model="BD",N=0,trackinfecteds=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{If n>0, then the simulation is stopped once 'n' tips are sampled sequentially through time. If n=0, the simulation is stopped after time 'timestop'.}
  \item{lambdasky}{Vector of dimension k, where k is the number of different birth rates. An individual between time (timesky[i],timesky[i+1]) has birth rate lambdasky[i].}
  \item{deathsky}{Vector of dimension k, where k is the number of different death rates. An individual between time (timesky[i],timesky[i+1]) has death rate deathsky[i].}
  \item{timesky}{Vector of dimension k, containing the times of rate shifts. Time is measured forward in time (unlike the function sim.rateshift.taxa where shifts are measured backward in time), with the origin of the tree being at time 0, i.e. timesky[1]=0.}
  \item{sampprobsky}{Vector of dimension k, an individual dying during time (timesky[i],timesky[i+1]) is sampled with probability sampprobsky[i], i.e. is being included into the final tree.}
  \item{omegasky}{Leave to default unless SIRS model simulation is being performed. omegasky is a vector of dimension k, where k is the number of different loose immunity rates (i.e. the rates of R->S transition in SIRS model). An individual between time (timesky[i],timesky[i+1]) has loose immunity rate omegasky[i].}
  \item{rho}{Default is rho=0. If rho>0 and timestop>0, then the process is stopped after timestop and each individual alive at time timestop is included into the final tree with probability rho.}
  \item{timestop}{Default is timestop=0, meaning the simulation is stopped once n tips are sampled. If timestop>0, then the simulation is stopped after time timestop.}
  \item{N}{Total population size is N. Set N>0 when simulating under either of SIS/SIR/SIRS models.}
  \item{model}{Should be set to desired model. Default="BD" (birth-death skyline model). The paramter accepts values "BD","SIS","SIR" or "SIRS". For all the models but "BD" N>0 should be set.}
  \item{trackinfecteds}{Set to TRUE records prevalence and incidence data, i.e. number of overall infecteds and times of infections in epidemiological terms, or overall number of species in that clade since time of origin, and times of speciations in macroevolutionary terms.}	
}
%\details{
%todo}
\value{
\item{out}{List containing the phylogenetic tree with n sampled tips or a fixed age timestop for trackinfecteds=FALSE. If trackinfecteds=TRUE, the list contains also a second item, a list that tracks numbers of susceptible/infected/recovered/sampled individuals over the course of the tree growth. This list consist of: $timesky - times of rate changes, $eventtimes - times dating events happening in the tree, i.e. bifurcation, death, sampling, or loose immunity, $infecteds - the number of infected infividuals at $eventtimes, $cumulativeinfecteds - the cumulative number of infected individuals at $eventtimes, $cumulativesampleds - the cumulative number of sampled individuals at $eventtimes, and in case of SIR/SIRS model, the list also contains: $susceptibles - the number of susceptible individuals at $eventtimes, and $recovereds - the number of recovered individuals at $eventtimes. Times $timesky and $eventtimes are stated backward-in-time, such that time=0 is the time of the most recent sample, and time is increasing into the past. This allows for determining precisely the times of rate changes for skyline tree analyses. The counts in $infecteds,$cumulativeinfecteds, $cumulativesampleds, $susceptibles and $recovereds represent the number of individuals in each category prior to (more ancestral than) the $eventtimes.}
}
\references{ T. Stadler, D. Kuehnert, S. Bonhoeffer, A. Drummond. Birth-death skyline plot reveals temporal changes of epidemic spread in HIV and hepatitis C virus (HCV). Proc. Nat. Acad. Sci., 110(1): 228-233, 2013.
		
 V. Boskova, S. Bonhoeffer, T. Stadler. Inference of epidemiological dynamics based on simulated phylogenies using birth-death and coalescent models. Manuscript.
}
\author{Tanja Stadler, Veronika Boskova}
\note{ 
A large number of trees can be obtained using the R function lapply. The tree can be plotted using the R package ape function plot(tree).  sim.bdsky.stt function extends the function sim.rateshift.taxa to trees which contain tips being sampled sequentially.}
\seealso{sim.bdtypes.stt.taxa}
\examples{
### Set the values for birth rates (lambda), deathrates (mu),
# sampling proportion (sampprob) and times of rate shifts (times).
# Also set the number of sampled tips in the final tree (n) and
# the number of simulations (numbsim).
set.seed(1)
n<-10
lambda <- c(2,1,2)
mu <- c(1,0.5,1.5)
sampprob <-c(0.5,0.5,0.5)
times<-c(0,1,2)
numbsim<-2
# Simulate trees under the birth-death skyline model
trees<-lapply(rep(n,numbsim),sim.bdsky.stt,lambdasky=lambda,deathsky=mu,
timesky=times,sampprobsky=sampprob,rho=0,timestop=0)

### Simulate 10 trees with 100 tips under the SIRS model with 
# total population size N=500
trees<-lapply(rep(100,10),sim.bdsky.stt,lambdasky=c(3,0.5,3,0.5,3),
deathsky=c(0.5,0.5,0.5,0.5,0.5),sampprobsky=c(0.5,0.5,0.5,0.5,0.5),
timesky=c(0,1,2,3,4),trackinfecteds=TRUE,model="SIRS",N=500,
omegasky=c(0,0.5,0.5,0.5,0))

### Simulate 1 tree with 100 tips under the SIRS model with 
# total population size N=500 and plot the S,I,R classes
trees<-sim.bdsky.stt(100,lambdasky=c(3,0.5,3,0.5,3),deathsky=c(0.5,0.5,0.5,0.5,0.5),
sampprobsky=c(0.5,0.5,0.5,0.5,0.5),timesky=c(0,2,2.5,3,3.2),trackinfecteds=TRUE,
model="SIRS",N=500,omegasky=c(0,0.5,0.5,0.5,0.5))
plot(trees[[2]]$eventtimes,trees[[2]]$infecteds,xlim=rev(range(trees[[2]]$eventtimes)),
type="l",col="red",ylim=c(min(trees[[2]]$recovereds,trees[[2]]$infecteds,trees[[2]]$susceptibles),
max(trees[[2]]$recovereds,trees[[2]]$infecteds,trees[[2]]$susceptibles)),
xlab="time",ylab="Number of individuals")
abline(v=trees[[2]]$timesky,lty=2)
points(trees[[2]]$eventtimes,trees[[2]]$recovereds,type="l",col="green")
points(trees[[2]]$eventtimes,trees[[2]]$susceptibles,type="l",col="blue")
points(trees[[2]]$eventtimes,trees[[2]]$cumulativesampleds,type="l",col="grey")
legend("topleft",c("S","I","R","samples","rate changes"),
col=c("blue","red","green","grey","black"),lty=c(1,1,1,1,2))
}