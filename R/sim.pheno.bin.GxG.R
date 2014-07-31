#'
#' @title Generates phenotype status
#' @description Generates affected and non-affected subjects
#' @param num.obs Number of observations to generate per iteration
#' @param disease.prev Prevalence of the binary outcome
#' @param genotype1 Exposure data for 1st genetic determinant
#' @param genotype2 Exposure data for 2st genetic determinant
#' @param interaction data
#' @param subject.effect.data Subject effect data, reflects the heterogenity 
#' in baseline disease risk
#' @param geno1.OR Odds ratios of the 1st genetic determinant
#' @param geno2.OR Odds ratios of the 2st genetic determinant
#' @param int.OR Odds ration of the interaction
#' @return A dataframe of phenotype
#' @keywords internal
#' @author Gaye A.
#'
sim.pheno.bin.GxG <-
function(num.obs=NULL, disease.prev=NULL, genotype1=NULL, genotype2=NULL, interaction=NULL, 
         subject.effect.data=NULL, geno1.OR=NULL, geno2.OR=NULL, int.OR=NULL)
{ 
   # GET THE ALPHA AND BETA VALUES
   alpha <- log(disease.prev/(1-disease.prev))
   geno1.beta <-	log(geno1.OR)
   geno2.beta <-  log(geno2.OR)
   int.beta <- log(int.OR)

   # GENERATE THE LINEAR PREDICTOR
   lp <- alpha + (geno1.beta*genotype1) + (geno2.beta*genotype2) + (int.beta*interaction) + subject.effect.data
   # GET THE 'mu' THE PROBABILITY OF DISEASE THROUGH LOGISTIC TRANSFORMATION
   mu <- exp(lp)/(1 + exp(lp))
   
   # GENERATE THE PHENOTYPE DATA AND RETURN IT AS A DATAFRAME
   phenotype <- rbinom(num.obs,1,mu)
   
   return(phenotype)
}

