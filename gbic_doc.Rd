\name{gbic}
\alias{gbic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Structureal Bayesian information criterion for mutinormal data with given graph structure
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gbic(data, theta, prob, P)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{x} here~~
n by p dataframe representing observations 
}
\item{theta}{p by p matrix representing the given graph structure}
\item{prob}{expected error rate}
\item{P}{prior adjacency matrix}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{sbic}{the value of sbic with given temperature parameter and prior adjaceency matrix}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
set.seed(1)
d=simulate(n=100, p=100, m1 = 100, m2 = 30)
data=d$data
P=d$priornetwork
theta=d$realnetwork
prob=0.15
index=gbic(data=data, theta=theta, prob=prob, P=P)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
