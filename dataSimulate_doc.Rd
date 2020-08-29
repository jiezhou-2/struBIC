\name{simulate}
\alias{simulate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Randomly generate  a adjacency matrix based on which simulated data are generated
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
According to a given edge density, first generate the adjacency matrix P of a graph. Based on P, the simulated multinorml data are generated with  mean zero and given precision matrix}
\usage{
simulate(p,n,prob1,prob2,ka)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
the number of vertices in graph
}
\item{n}{sample size}
\item{m1}{the number of edges in the true graph}
\item{m2}{the number of elements in adjacency matrix that stay in different states, i.e., 0 or 1, in true and prior graph.}
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
a list including the simulated data, real adjacency matrix and a prior adjacency matrix
\item{data}{simulated data}
\item{realnetwork}{real adjacency matrix}
\item{priornetowrk}{prior adjacency matrix}
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
d=simulate(n=100,p=200, m1=100, m2=30)
d$data
d$realnetwork
d$priornetwork
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
