\name{get.critical.results.GxG}
\alias{get.critical.results.GxG}
\title{Summarizes the main results}
\usage{
get.critical.results.GxG(scenario = NULL, pheno.model = NULL,
  geno1.model = NULL, geno2.model = NULL, sample.sizes.required = NULL,
  empirical.power = NULL, modelled.power = NULL, mean.beta = NULL)
}
\arguments{
  \item{scenario}{Scenario number}

  \item{pheno.model}{Type of the outcome; 0 for binary and
  2 for continuous}

  \item{geno1.model}{Genetic model of the 1st SNP; 0 for
  binary and 1 for additive}

  \item{geno2.model}{Genetic model of the 1st SNP; 0 for
  binary and 1 for additive}

  \item{sample.sizes.required}{Number of cases and controls
  or number of subjects required to achieve the desired
  power}

  \item{empirical.power}{Estimated empirical power}

  \item{modelled.power}{Calculated theoretical power}

  \item{mean.beta}{Mean beta value of each of the
  determinants}
}
\value{
A table containing the following items: \code{gene1.model}
Model of the first genetic determinant \code{gene2.model}
Model of the second genetic determinant
\code{number.of.cases.required} Number of cases required to
achieve the desired power under binary outcome model
\code{number.of.controls.required} Number of controls
required to achieve the desired power under binary outcome
model \code{number.of.subjects.required} Number of subjects
required to achieve the desired power under a quantatative
outcome model \code{empirical.power} Estimated empirical
power under each model \code{modelled.power} Power achieved
under each model with specified sample size
\code{estimated.OR} Estimated odds-ratios due to shrinkage
toward the null resulting from misclassification
\code{estimated.effect} Estitmated effect size if the
outocme is continuous
}
\description{
Gets the number of cases and controls or subjects and the
empirical and theoretical power under each model and prints
a summary on the screen
}
\details{

}
\author{
Gaye A.
}
\keyword{internal}

