#' 
#' @title Calculates the empirical and theoretical power
#' @description The function determines the empirical and theoretical power. 
#' The empirical power is the proportion of simulations in which 
#' the z-statistic for the parameter of interest exceeds the z-statistic 
#' for the desured level if statistical significance. 
#' The theoretical power is the power of the study.
#' @param pval cut-off p-value defining statistical significance.
#' @param z.values z-statistic of the determinant.
#' @param mean.model.z mean z-statistic of the environmental determinant.
#' @return a list that contains the computed empirical power and theoretical power.
#' @keywords internal
#' @author Gaye A.
#'
power.calc <- function(pval=NULL, z.values=NULL, mean.model.z=NULL){
  
  if(is.null(z.values)){
    cat("\n\n ALERT!\n")
    cat(" No z-statistics found\n")
    cat(" Check the argument 'z.values'.\n")
    stop(" End of process!\n\n", call.=FALSE)
  }
  
  if(is.null(mean.model.z)){
    cat("\n\n ALERT!\n")
    cat(" The argument 'mean.model.z' is set to NULL.\n")
    cat(" This argument should be the ratio 'mean.beta/mean.se'.\n")
    stop(" End of process!\n\n", call.=FALSE)
  }
  
  # CALCULATE Z STATISTIC THRESHOLD FOR DESIRED P-VALUE 
  z.pval <- qnorm(1-pval/2)
  
  # GET EMPIRICAL POWER: THE PROPORTION OF SIMULATIONS IN WHICH THE 
  # Z STATISTIC FOR THE PARAMETER OF INTEREST EXCEEDS THE Z STATISTIC 
  # FOR THE DESIRED LEVEL OF STATISTICAL SIGNIFICANCE
  empirical.power <- round(mean((z.values > z.pval), na.rm=TRUE),3)
  
  # GET THE MODELLED POWER
  modelled.power <- pnorm(mean.model.z-z.pval)
  
  return(list(empirical=empirical.power, modelled=modelled.power))
}
