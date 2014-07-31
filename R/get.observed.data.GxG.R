#'
#' @title Generates exposure data with some error
#' @description Uses THE function make.obs.geno to generate effect data with a set level of error
#' @param data Input table of simulated data considered as true data
#' @param g1.error Misclassification rates in the assessment of the first genetic variant: 1-sensitivity and 1-specificity
#' @param g1.model Genetic model of the first genetic variant; 0 for binary and 1 for additive
#' @param freq1 Minor allele frequency of the first genetic variant
#' @param g2.error Misclassification rates in the assessment of the second genetic variant: 1-sensitivity and 1-specificity
#' @param g2.model Genetic model of the second genetic variant; 0 for binary and 1 for additive
#' @param freq2 Minor allele frequency of the second genetic variant
#' @return A matrix
#' @keywords internal
#' @author Gaye A
#'
get.observed.data.GxG <-
function(data=NULL,g1.error=NULL,g1.model=NULL,freq1=NULL,
         g2.error=NULL,g2.model=NULL,freq2=NULL)
{
		sim.df <- data      

    # GET THE OBSERVED GENOTYPES FOR THE FIRST SNP
    true.genotype1 <- sim.df$genotype1
    obs.genotype1 <- get.obs.geno(allele.A=sim.df$allele.A1,allele.B=sim.df$allele.B1,
                                 geno.model=g1.model,MAF=freq1,geno.error=g1.error)
    
		# GET THE OBSERVED GENOTYPES FOR THE SECOND SNP
		true.genotype2 <- sim.df$genotype2
		obs.genotype2 <- get.obs.geno(allele.A=sim.df$allele.A2,allele.B=sim.df$allele.B2,
		                              geno.model=g2.model,MAF=freq2,geno.error=g2.error)
		
    # GET THE OBSERVED INTERACTION DATA
    obs.interaction <- obs.genotype1$observed.genotype * obs.genotype2$observed.genotype

    # REPLACE THE TRUE DATA BY THE NOW GENERATED OBSERVED GENOTYPES
    # IN THE INITIAL MATRIX THAT HELD THE TRUE DATA
    sim.df$genotype1 <- obs.genotype1$observed.genotype
		sim.df$genotype2 <- obs.genotype2$observed.genotype
    sim.df$allele.A1 <- obs.genotype1$observed.allele.A
    sim.df$allele.B1 <- obs.genotype1$observed.allele.B
		sim.df$allele.A2 <- obs.genotype2$observed.allele.A
		sim.df$allele.B2 <- obs.genotype2$observed.allele.B
		sim.df$interaction <- obs.interaction
    
		# RETURN THE MATRIX WHICH NOW CONTAINS ONLY THE OBSERVED DATA TO ANALYSE BY GLM
    colnames(sim.df) <- c("id", "phenotype", "genotype1", "allele.A1", "allele.B1", "genotype2", "allele.A2", "allele.B2", "interaction")
    return(sim.df)
}

