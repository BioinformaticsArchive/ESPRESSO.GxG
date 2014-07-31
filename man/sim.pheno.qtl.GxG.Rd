\name{sim.pheno.qtl.GxG}
\alias{sim.pheno.qtl.GxG}
\title{Simulates continuous outcome data}
\usage{
sim.pheno.qtl.GxG(numsubjects = NULL, pheno.mean = NULL, pheno.sd = NULL,
  genotype1 = NULL, geno1.efkt = NULL, genotype2 = NULL,
  geno2.efkt = NULL, interaction = NULL, int.efkt = NULL)
}
\arguments{
  \item{num.subjects}{Number of subjects to simulate}

  \item{pheno.mean}{statistical mean}

  \item{pheno.sd}{standard deviation}

  \item{genotype1}{Genotypes of the 1st variant}

  \item{genotype2}{Genotypes of the 2nd variant}

  \item{geno1.efkt}{Effects of the 1st genetic variant}

  \item{geno2.efkt}{Effects of the 2nd genetic variant}

  \item{interaction}{Effect data}

  \item{int.efkt}{Interaction effect}
}
\value{
A dataframe of phenotype
}
\description{
Uses the effects data of the determinants to construct a
linear predictor(LP). The outcome is normally distributed
variable generated with a mean equal to LP and a standard
deviation of 1. Some error is then added to the simulated
outcome to obtained the observed outcome.
}
\details{

}
\author{
Gaye A.
}
\keyword{internal}

