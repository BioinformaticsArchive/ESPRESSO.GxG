\name{sim.pheno.bin.GxE}
\alias{sim.pheno.bin.GxE}
\title{Generates phenotype status}
\usage{
sim.pheno.bin.GxE(num.obs = NULL, disease.prev = NULL, genotype = NULL,
  environment = NULL, interaction = NULL, subject.effect.data = NULL,
  geno.OR = NULL, env.OR = NULL, int.OR = NULL)
}
\arguments{
  \item{num.obs}{Number of observations to generate per
  iteration}

  \item{disease.prev}{Prevalence of the binary outcome}

  \item{genotype}{Exposure data for genetic determinate}

  \item{environment}{Exposure data for environment}

  \item{interaction}{Effect model: main effects=0,
  Gene-Environment interaction=1, Gene-Gene interaction=2
  and Environment-Enviroment interaction=3}

  \item{subject.effect.data}{Subject effect data, reflects
  the heterogenity in baseline disease risk}

  \item{geno.OR}{Odds ratios of the two genetic
  determinates}

  \item{env.OR}{Odds ratios of the two environments}

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

