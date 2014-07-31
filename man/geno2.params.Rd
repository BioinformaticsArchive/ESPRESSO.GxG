\name{geno2.params}
\alias{geno2.params}
\docType{data}
\title{Parameters to generate genetic exposures data}
\description{A table that contains parameters to set the genetic exposure}
\usage{data(geno2.params)}
\format{
  A data frame with 4 observations on the following 7 variables.
  \describe{
    \item{\code{scenario.id}}{scenario ID}
    \item{\code{geno2.model}}{genetic model: binary (0) or additive (1)}
    \item{\code{MAF2}}{minor allele frequency of the SNP}
    \item{\code{geno2.OR}}{odds ratio if outcome is binary}
    \item{\code{geno2.efkt}}{effect size if outcome is quantitative}
    \item{\code{geno2.sensitivity}}{sensitivity of the assessment of the SNP alleles}
    \item{\code{geno2.specificity}}{specificity of the assessment of the SNP alleles}
  }
}
\examples{data(geno2.params)}

