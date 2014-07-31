#'
#' @title Simulates subjects for continuous outcome
#' @description Generates the specified number of subjects
#' @param n Number of subjects to simulate
#' @param ph.mean statistical mean
#' @param ph.sd standard deviation
#' @param freq1 Minor allele frequency if the 1st of the two genetic variant
#' @param freq2 Minor allele frequency if the 2nd of the two genetic variant
#' @param g1.model Genetic model; 0 for binary and 1 for additive
#' @param g1.efkt Effects of the 1st genetic variant
#' @param g2.model Genetic model; 0 for binary and 1 for additive
#' @param g2.efkt Effects of the 2nd genetic variant
#' @param i.efkt Interaction effect
#' @param pheno.rel reliability of the assessment for a quantitative outcome.
#' @return A matrix
#' @keywords internal
#' @author Gaye A.

sim.QTL.data.GxG <-
function(n=NULL,ph.mean=NULL,ph.sd=NULL,freq1=NULL,freq2=NULL,g1.model=NULL,g1.efkt=NULL,
         g2.model=NULL,g2.efkt=NULL,i.efkt=NULL,pheno.rel=NULL)
{
	   
   # GENERATE THE TRUE GENOTYPE DATA FOR THE 1st VARIANT
   geno1.data <- sim.geno.data(num.obs=n, geno.model=g1.model, MAF=freq1)
   allele.A1 <- geno1.data$allele.A
   allele.B1 <- geno1.data$allele.B
   geno1 <- geno1.data$genotype
			
   # GENERATE THE TRUE GENOTYPE DATA FOR THE 2nd VARIANT
   geno2.data <- sim.geno.data(num.obs=n, geno.model=g2.model, MAF=freq2)
   allele.A2 <- geno2.data$allele.A
   allele.B2 <- geno2.data$allele.B
   geno2 <- geno2.data$genotype
          
   # GENERATE THE TRUE INTERACTION DATA           
   int <- geno1 * geno2
           
   # GENERATE THE TRUE OUTCOME DATA
   pheno.data <- sim.pheno.qtl.GxG(numsubjects=n,pheno.mean=ph.mean,pheno.sd=ph.sd,
                                   genotype1=geno1,genotype2=geno2,geno1.efkt=g1.efkt,geno2.efkt=g2.efkt,
                                   interaction=int,int.efkt=i.efkt)
   true.phenotype <- pheno.data
   
   # GENERATE THE OBSERVED OUTCOME DATA 
   obs.phenotype <- get.obs.pheno(phenotype=true.phenotype, pheno.model=1, 
                                  pheno.sd=ph.sd, pheno.reliability=pheno.rel)
   pheno <- obs.phenotype
   

   # STORE THE GENERATED TRUE DATA INTO AN OUTPUT MATRIX 
   sim.matrix <- cbind(pheno,geno1,allele.A1,allele.B1,geno2,allele.A2,allele.B2,int)

   # ADD IDs (JUST A ROW COUNT)
   totalnumrows <- dim(sim.matrix)[1]
   sim.matrix <- cbind(1:totalnumrows, sim.matrix)

   # ADD COLUMN NAMES AND RETURN A DATAFRAME
   colnames(sim.matrix) <- c("id","phenotype","genotype1","allele.A1","allele.B1","genotype2","allele.A2","allele.B2","interaction")
   mm <- data.frame(sim.matrix)
}

