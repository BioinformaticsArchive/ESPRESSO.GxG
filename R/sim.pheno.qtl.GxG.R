#' 
#' @title Simulates continuous outcome data
#' @description Uses the effects data of the determinants to construct a linear predictor(LP). The outcome is normally distributed variable generated with a mean equal to LP and a standard deviation of 1. Some error is then added to the simulated outcome to obtained the observed outcome.
#' @param num.subjects Number of subjects to simulate
#' @param pheno.mean statistical mean
#' @param pheno.sd standard deviation
#' @param genotype1 Genotypes of the 1st variant
#' @param genotype2 Genotypes of the 2nd variant
#' @param geno1.efkt Effects of the 1st genetic variant
#' @param geno2.efkt Effects of the 2nd genetic variant
#' @param interaction Effect data
#' @param int.efkt Interaction effect
#' @return A dataframe of phenotype
#' @keywords internal
#' @author Gaye A.
#'
sim.pheno.qtl.GxG <-
function(numsubjects=NULL,pheno.mean=NULL,pheno.sd=NULL,genotype1=NULL,geno1.efkt=NULL,
         genotype2=NULL,geno2.efkt=NULL,interaction=NULL,int.efkt=NULL)
{  

  # ALPHA IS EQUAL TO THE MEAN OF THE TRAIT, WHICH IS 0
   num.obs <- numsubjects
   alpha <- pheno.mean 

   # GENERATE THE LINEAR PREDICTOR
   lp <- alpha + (geno1.efkt*genotype1) + (geno2.efkt*genotype2) + (int.efkt*interaction)

   # GENERATE THE TRUE PHENOTYPE DATA
   phenotype <- rnorm(num.obs,lp,pheno.sd)
   
   return(phenotype)
}

