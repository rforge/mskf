\name{mskf-utilities}
\alias{paW}
\alias{pac}
\alias{paH}
\alias{paG}
\alias{paK}
\alias{pap}
\alias{maW}
\alias{mac}
\alias{maH}
\alias{maG}
\alias{maK}
\alias{map}
\alias{paW<-}
\alias{pac<-}
\alias{paH<-}
\alias{paG<-}
\alias{paK<-}
\alias{pap<-}
\alias{maW<-}
\alias{mac<-}
\alias{maH<-}
\alias{maG<-}
\alias{maK<-}
\alias{map<-}
\alias{startValues.mskfModel}
\alias{startValues}
\alias{startValues<-}
\alias{lower}
\alias{lower<-}
\alias{upper}
\alias{upper<-}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Utility Function For Handling And Manipulating \code{mskfModel} Objects.
}
\description{
Extract model matrices from a \code{mskfModel} object. Assign values to model matrices in a \code{mskfModel} object.
}
\usage{
paW(x)
paK(x)
paG(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Object of class \code{mskfModel}}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Arrays with the model matrices.
}
\references{
Hamaker & Grasman (2012)
}
\author{
Raoul P. P. P. Grasman
}
\note{
No further notes.
}


\seealso{
\code{\link{mskf}}, \code{\link{mskfModel}}, \code{\link{mskfModelSkeleton}}
}
\examples{
x = list(paW=NA)
paW(x)
}



