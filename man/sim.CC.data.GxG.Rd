\name{sim.CC.data.GxG}
\alias{sim.CC.data.GxG}
\title{Simulates case and controls}
\usage{
sim.CC.data.GxG(n = NULL, ncases = NULL, ncontrols = NULL,
  max.sample.size = NULL, pheno.prev = NULL, freq1 = NULL,
  g1.model = NULL, g1.OR = NULL, freq2 = NULL, g2.model = NULL,
  g2.OR = NULL, i.OR = NULL, b.OR = NULL, ph.error = NULL)
}
\arguments{
  \item{n}{Number of observations to generate per
  iteration}

  \item{ncases}{Number of cases to simulate}

  \item{ncontrols}{Number of controls to simulate}

  \item{max.sample.size}{Maximum number of observations
  allowed}

  \item{pheno.prev}{Prevalence of the binary outcome}

  \item{freq1}{Minor allele frequency of the 1st genetic
  determinant}

  \item{g1.model}{Genetic model of the 1st genetic
  determinant; 0 for binary and 1 for additive}

  \item{g1.OR}{Odds ratios of the 1st genetic determinant}

  \item{freq2}{Minor allele frequency of the 2st genetic
  determinant}

  \item{g2.model}{Genetic model of the 2st genetic
  determinant; 0 for binary and 1 for additive}

  \item{g2.OR}{Odds ratios of the 2st genetic determinant}

  \item{i.OR}{Odds ration of the interaction}

  \item{b.OR}{Baseline odds ratio for subject on 95 percent
  population centile versus 5 percentile. This parameter
  reflects the heterogeneity in disease risk arising from
  determinates that have not been measured or have not been
  included in the model.}

  \item{ph.error}{misclassification rates: 1-sensitivity
  and 1-specificity}
}
\value{
A matrix
}
\description{
Generates affected and non-affected subjects until the set
sample size is achieved.
}
\details{

}
\author{
Gaye A.
}
\keyword{internal}

