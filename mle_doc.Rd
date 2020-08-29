\name{mle}
\alias{mle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Estimate the precision matrix for multinormal distribution with given adjacency matrix
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function find the maximum likelihood estimate of  the precision matrix with given adjacency matrix  for multinormal distribution.  
}
\usage{
mle(
data,
priori
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a p by p datataframe representing the observations 
%%     ~~Describe \code{x} here~~
}
\item{priori}{
a p by p matrix representing the given adjacency matrix
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
The methods  are based on the relationship between precision matrix of multinormal distribution and regression ecoefficients.  
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
A p by p matrix representing the estimate of the precision matrix. 
%% ...
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
d=simulate(p=10,n=100, prob1 = 0.1, prob2 = 0.2,ka=4)
data=d$data
priori=d$realnetwork
precision=mle(data=data,priori=priori)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
