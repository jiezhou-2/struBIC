\name{modelset}
\alias{modelset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Construct the model pool based on two-step algorithm
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
For a given prior graph, the two-step algorithm, including edge enrichment and pruning,  is used to construct the model pool. 
}
\usage{
modelset(data, lambda, P)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{x} here~~
n by p data frame 
}
\item{lambda}{tuning parmeter vector}
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
\item{candidates}{a list including all the candidate models in the model pool. Each model is represented by a p by p adjacency matrix} }
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
library("glmnet")
set.seed(1)
d=simulate(n=100, p=100, m1 = 100, m2 = 30)
data=d$data
P=d$priornetwork
lambda=exp(seq(-5,5,length=100))
candidates=modelset(data=data,lambda=lambda, P=P)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
