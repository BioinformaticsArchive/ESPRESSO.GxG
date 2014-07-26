#' 
#' @title Carries out regression analysis
#' @description Fits a conventional unconditional logistic regression model with a binary or continuous phentype as outcome and the genetic, environmental, interaction determinants as covariates.
#' @param pheno.model Type of outcome; 0=binary and 1=continuous
#' @param observed.data A dataframe that contains covariates and outcome data
#' @return A vector containing the beta, standard-error and z-statistic of each of the covariates
#' @keywords internal
#' @author Gaye A.
#'
glm.analysis.GxE <-
function(pheno.model=NULL, observed.data=NULL)
{

  # BINARY OUTCOME
  if(pheno.model == 0){
	   # FIT CONVENTIONAL UNCONDITIONAL LOGISTIC REGRESSION MODEL
	   mod.glm <- glm(phenotype ~ 1+genotype+environment+interaction,
	                  family=binomial,data=observed.data)
	   mod.sum <- summary(mod.glm)
  }
  
  # QUANTITATIVE OUTCOME     
  if(pheno.model == 1){
	    # FIT A GLM FOR A GAUSSIAN OUTCOME
	    mod.glm <- glm(phenotype ~ 1+genotype+environment+interaction,
	    family=gaussian,data=observed.data)
	    mod.sum <- summary(mod.glm)     
  }
  
	 beta.value <- mod.sum$coefficients[4,1]
	 se.value <- mod.sum$coefficients[4,2]
	 z.value <- mod.sum$coefficients[4,3]
	 
  # RETURN A VECTOR
  return(list(beta=beta.value, se=se.value, z=z.value))
}

