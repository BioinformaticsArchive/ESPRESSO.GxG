#'
#' @title Simulates subjects for continuous outcome
#' @description Generates the specified number of subjects
#' @param n Number of subjects to simulate
#' @param ph.mean statistical mean
#' @param ph.sd standard deviation
#' @param freq Minor allele frequencies of the two genetic variants
#' @param g.model Genetic model; 0 for binary and 1 for continuous
#' @param g.efkt Effects of the genetic variants
#' @param e.model Model of the environmental exposure
#' @param e.efkt Effects of the environment determinats
#' @param e.prev Prevalences of the environmental determinants
#' @param e.mean Mean under quantitative-normal model
#' @param e.sd Standard deviation under quantitative-normal model
#' @param e.low.lim lower limit under quantitative-uniform model
#' @param e.up.lim upper limit under quantitative-uniform model
#' @param i.efkt Interaction effect
#' @param pheno.rel reliability of the assessment for a quantitative outcome.
#' @return A matrix
#' @keywords internal
#' @author Gaye A.

sim.QTL.data.GxE <-
function(n=NULL,ph.mean=NULL,ph.sd=NULL,freq=NULL,g.model=NULL,g.efkt=NULL,e.model=NULL,
         e.efkt=NULL,e.prev=NULL,e.mean=NULL,e.sd=NULL,e.low.lim=NULL,e.up.lim=NULL,i.efkt=NULL,
         pheno.rel=NULL)
{
	   
   # GENERATE THE TRUE GENOTYPE DATA
   geno.data <- sim.geno.data(num.obs=n, geno.model=g.model, MAF=freq)
   allele.A <- geno.data$allele.A
   allele.B <- geno.data$allele.B
   geno <- geno.data$genotype
			
   # GENERATE THE TRUE ENVIRONMENTAL EXPOSURE DATA			
   env <- sim.env.data(num.obs=n,env.model=e.model,env.prev=e.prev,env.mean=e.mean,
                       env.sd=e.sd,env.low.lim=e.low.lim,env.up.lim=e.up.lim)
          
   # GENERATE THE TRUE INTERACTION DATA           
   int <- geno*env
           
   # GENERATE THE TRUE OUTCOME DATA
   pheno.data <- sim.pheno.qtl.GxE(numsubjects=n,pheno.mean=ph.mean,pheno.sd=ph.sd,
                                   genotype=geno,geno.efkt=g.efkt,environment=env,env.efkt=e.efkt,
                                   interaction=int,int.efkt=i.efkt)
   true.phenotype <- pheno.data
   
   # GENERATE THE OBSERVED OUTCOME DATA 
   obs.phenotype <- get.obs.pheno(phenotype=true.phenotype, pheno.model=1, 
                                  pheno.sd=ph.sd, pheno.reliability=pheno.rel)
   pheno <- obs.phenotype
   

   # STORE THE GENERATED TRUE DATA INTO AN OUTPUT MATRIX 
   sim.matrix <- cbind(pheno,geno,allele.A,allele.B,env,int)

   # ADD IDs (JUST A ROW COUNT)
   totalnumrows <- dim(sim.matrix)[1]
   sim.matrix <- cbind(1:totalnumrows, sim.matrix)

   # ADD COLUMN NAMES AND RETURN A DATAFRAME
   colnames(sim.matrix) <- c("id","phenotype","genotype","allele.A","allele.B","environment","interaction")
   mm <- data.frame(sim.matrix)
}

