\name{sim.pheno.bin.GxG}
\alias{sim.pheno.bin.GxG}
\title{Generates phenotype status}
\usage{
sim.pheno.bin.GxG(num.obs = NULL, disease.prev = NULL, genotype1 = NULL,
  genotype2 = NULL, interaction = NULL, subject.effect.data = NULL,
  geno1.OR = NULL, geno2.OR = NULL, int.OR = NULL)
}
\arguments{
  \item{num.obs}{Number of observations to generate per
  iteration}

  \item{disease.prev}{Prevalence of the binary outcome}

  \item{genotype1}{Exposure data for 1st genetic
  determinant}

  \item{genotype2}{Exposure data for 2st genetic
  determinant}

  \item{interaction}{data}

  \item{subject.effect.data}{Subject effect data, reflects
  the heterogenity in baseline disease risk}

  \item{geno1.OR}{Odds ratios of the 1st genetic
  determinant}

  \item{geno2.OR}{Odds ratios of the 2st genetic
  determinant}

  \item{int.OR}{Odds ration of the interaction}
}
\value{
A dataframe of phenotype
}
\description{
Generates affected and non-affected subjects
}
\details{

}
\author{
Gaye A.
}
\keyword{internal}

