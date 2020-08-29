\name{sggm}
\alias{sggm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Model selection of Gaussian graphical model based on SBIC
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
select the model  based on the SBIC and the model pool construted by two-step algorithm
}
\usage{
sggm(data=data, lambda, M, prob)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{x} here~~
n by p dataframe representing the observations
}
\item{lambda}{tuning parameter vector used to build the model pool}
\item{M}{the prior adjacency matrix}
\item{prob}{the mean error rate}
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
\item{networkhat}{the final selected adjacency matrix}
\item{candidates}{the model pool}
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
d=simulate(n=100,p=100 m1 = 100, m2 = 30)
lambda=exp(seq(-5,5,length=100))
data=d$data
M=d$priornetwork
r1=m2/m1
r2=m2/(p*(p-1)/2-m1)
r=(r1+r2)/2
result=sggm(data=data, lambda=lambda, M=M,prob=r)
result$networkhat
d$realnetwork
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
