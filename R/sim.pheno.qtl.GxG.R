#' 
#' @title Simulates continuous outcome data
#' @description Uses the effects data of the determinants to construct a linear predictor(LP). The outcome is normally distributed variable generated with a mean equal to LP and a standard deviation of 1. Some error is then added to the simulated outcome to obtained the observed outcome.
#' @param num.subjects Number of subjects to simulate
#' @param pheno.mean statistical mean
#' @param pheno.sd standard deviation
#' @param genotype Genotype
#' @param geno.efkt Effects if the genetic variants
#' @param environment Exposure data for environment
#' @param env.efkt Effects of the environmental determiants
#' @param interaction Effect model: main effects=0, Gene-Environment interaction=1, Gene-Gene interaction=2 and Environment-Enviroment interaction=3
#' @param int.efkt Interaction effect
#' @return A dataframe of phenotype
#' @keywords internal
#' @author Gaye A.
#'
sim.pheno.qtl.GxE <-
function(numsubjects=NULL,pheno.mean=NULL,pheno.sd=NULL,genotype=NULL,geno.efkt=NULL,
         environment=NULL,env.efkt=NULL,interaction=NULL,int.efkt=NULL)
{  

  # ALPHA IS EQUAL TO THE MEAN OF THE TRAIT, WHICH IS 0
   num.obs <- numsubjects
   alpha <- pheno.mean 

   # GENERATE THE LINEAR PREDICTOR
   lp <- alpha + (geno.efkt*genotype) + (env.efkt*environment) + (int.efkt*interaction)

   # GENERATE THE TRUE PHENOTYPE DATA
   phenotype <- rnorm(num.obs,lp,pheno.sd)
   
   return(phenotype)
}

