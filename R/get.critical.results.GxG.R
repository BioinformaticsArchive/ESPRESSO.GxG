#'
#' @title Summarizes the main results
#' @description Gets the number of cases and controls or subjects and the empirical and theoretical power under each model and prints a summary on the screen
#' @param scenario Scenario number
#' @param pheno.model Type of the outcome; 0 for binary and 2 for continuous
#' @param geno1.model Genetic model of the 1st SNP; 0 for binary and 1 for additive
#' @param geno2.model Genetic model of the 1st SNP; 0 for binary and 1 for additive
#' @param sample.sizes.required Number of cases and controls or number of subjects required to achieve the desired power
#' @param empirical.power Estimated empirical power
#' @param modelled.power Calculated theoretical power
#' @param mean.beta Mean beta value of each of the determinants
#' @return A table containing the following items:
#' \code{gene1.model} Model of the first genetic determinant
#' \code{gene2.model} Model of the second genetic determinant
#' \code{number.of.cases.required} Number of cases required to achieve the desired power under binary outcome model
#' \code{number.of.controls.required} Number of controls required to achieve the desired power under binary outcome model
#' \code{number.of.subjects.required} Number of subjects required to achieve the desired power under a quantatative outcome model
#' \code{empirical.power} Estimated empirical power under each model
#' \code{modelled.power} Power achieved under each model with specified sample size
#' \code{estimated.OR} Estimated odds-ratios due to shrinkage toward the null resulting from misclassification
#' \code{estimated.effect} Estitmated effect size if the outocme is continuous
#' @keywords internal
#' @author Gaye A.
#'
get.critical.results.GxG <-
function(scenario=NULL, pheno.model=NULL,geno1.model=NULL,geno2.model=NULL,sample.sizes.required=NULL,
         empirical.power=NULL,modelled.power=NULL,mean.beta=NULL)
{
		 
 	 if(geno1.model==0){
 	   g1.model <- "binary"
 	 }else{
 	   g1.model <- "additive"
 	 }
    
 	 if(geno2.model==0){
 	   g2.model <- "binary"
 	 }else{
 	   g2.model <- "additive"
 	 }

   if(pheno.model == 0){
			numcases <- sample.sizes.required[[1]]
			numcontrols <- sample.sizes.required[[2]]
   }else{
			numsubjects <- sample.sizes.required[[1]]
   }

  if(pheno.model==0){
     # estimated ORs
     estimated.OR <- exp(mean.beta)
     estimated.effect <- 'NA'

					cat("\n---- SUMMARY OF SCENARIO",scenario,"----\n")
					cat("\nModels\n")
					cat("------\n")
					cat("  Outcome: binary \n")
					cat("  First genetic determinant:",g1.model,"\n")
					cat("  Second genetic determinant:",g2.model)

					cat("\n\nNumber of cases required\n")
					cat("------------------------\n")
					cat(" ", numcases)

					cat("\n\nNumber of controls required\n")
					cat("---------------------------\n")
					cat(" ", numcontrols)

					cat("\n\nEmpirical power\n")
					cat("---------------\n")
					cat(" ",empirical.power)

					cat("\n\nModel power\n")
					cat("-----------\n")
					cat(" ",round(modelled.power,2))

					cat("\n\nEstimated interaction OR\n")
					cat("-----------\n")
					cat(" ",round(estimated.OR,2))


					cat("\n\n---- END OF SUMMARY ----\n")

		   crit.res <- c(g1.model,g2.model,numcases,numcontrols,round(empirical.power,2),round(modelled.power,2),
		                 round(estimated.OR,2), estimated.effect)
		   return(list(gene1.model=crit.res[1], gene2.model=crit.res[2],number.of.cases.required=crit.res[3],
		               number.of.controls.required=crit.res[4],empirical.power=crit.res[5], 
		               modelled.power=crit.res[6],estimated.OR=crit.res[7], estimated.effect=crit.res[8]))

  }else{
     # estimated ORs
     estimated.effect <- mean.beta
     estimated.OR <- 'NA'

					cat("\n---- SUMMARY OF SCENARIO",scenario,"----\n")
					cat("\nModels\n")
					cat("------\n")
					cat(" Outcome: quantitative; ")
					cat(" First genetic determinant: ",g1.model,"\n")
					cat(" Second genetic determinant: ",g2.model)

					cat("\n\nNumber of subjects required\n")
					cat("------------------------\n")
					cat(" ",numsubjects)

					cat("\n\nEmpirical power\n")
					cat("---------------\n")
					cat(" ",empirical.power)

					cat("\n\nModel power\n")
					cat("-----------\n")
					cat(" ",round(modelled.power,2))

					cat("\n\nEstimated interaction effect\n")
					cat("-----------\n")
					cat(" ",estimated.effect)

					cat("\n\n---- END OF SUMMARY ----\n")

		   crit.res <- c(g1.model,g2.model,numsubjects,round(empirical.power,2),round(modelled.power,2),estimated.OR,round(estimated.effect,2))
		   return(list(gene1.model=crit.res[1], gene2.model=crit.res[2],number.of.subjects.required=crit.res[3],
		              empirical.power=crit.res[4], modelled.power=crit.res[5], estimated.OR=crit.res[6],estimated.effect=crit.res[7]))
  }
}

