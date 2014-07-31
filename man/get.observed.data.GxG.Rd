\name{get.observed.data.GxG}
\alias{get.observed.data.GxG}
\title{Generates exposure data with some error}
\usage{
get.observed.data.GxG(data = NULL, g1.error = NULL, g1.model = NULL,
  freq1 = NULL, g2.error = NULL, g2.model = NULL, freq2 = NULL)
}
\arguments{
  \item{data}{Input table of simulated data considered as
  true data}

  \item{g1.error}{Misclassification rates in the assessment
  of the first genetic variant: 1-sensitivity and
  1-specificity}

  \item{g1.model}{Genetic model of the first genetic
  variant; 0 for binary and 1 for additive}

  \item{freq1}{Minor allele frequency of the first genetic
  variant}

  \item{g2.error}{Misclassification rates in the assessment
  of the second genetic variant: 1-sensitivity and
  1-specificity}

  \item{g2.model}{Genetic model of the second genetic
  variant; 0 for binary and 1 for additive}

  \item{freq2}{Minor allele frequency of the second genetic
  variant}
}
\value{
A matrix
}
\description{
Uses THE function make.obs.geno to generate effect data
with a set level of error
}
\details{

}
\author{
Gaye A
}
\keyword{internal}

