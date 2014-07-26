#' 
#' @title Simulates case and controls
#' @description Generates affected and non-affected subjects until the set sample 
#' size is achieved.
#' @param n Number of observations to generate per iteration
#' @param ncases Number of cases to simulate
#' @param ncontrols Number of controls to simulate
#' @param max.sample.size Maximum number of observations allowed
#' @param pheno.prev Prevalence of the binary outcome
#' @param freq Minor allele frequency
#' @param g.model Genetic model; 0 for binary and 1 for continuous
#' @param g.OR Odds ratios of the genetic determinants
#' @param e.model Model of the environmental exposure
#' @param e.prev Prevelance of the environmental determinates
#' @param e.mean Mean under quantitative-normal model
#' @param e.sd Standard deviation under quantitative-normal model
#' @param e.low.lim Lower limit under quantitative-uniform model
#' @param e.up.lim Upper limit under quantitative-uniform model
#' @param e.OR Odds ratios of the environmental determinants
#' @param i.OR Odds ration of the interaction
#' @param b.OR Baseline odds ratio for subject on 95 percent population 
#' centile versus 5 percentile. This parameter reflects the heterogeneity in disease 
#' risk arising from determinates that have not been measured or have not been 
#' included in the model.
#' @param ph.error misclassification rates: 1-sensitivity and 1-specificity
#' @return A matrix
#' @keywords internal
#' @author Gaye A.
#'
sim.CC.data.GxE <-
function(n=NULL, ncases=NULL, ncontrols=NULL, max.sample.size=NULL, pheno.prev=NULL,
freq=NULL, g.model=NULL, g.OR=NULL, e.model=NULL, e.prev=NULL, e.mean=NULL, e.sd=NULL, 
e.low.lim=NULL, e.up.lim=NULL, e.OR=NULL, i.OR=NULL, b.OR=NULL, ph.error=NULL)
{
   # SET UP ZEROED COUNT VECTORS TO DETERMINE WHEN ENOUGH CASES AND CONTROLS HAVE BEEN GENERATED
   complete <- 0
   complete.absolute <- 0
   cases.complete <- 0
   controls.complete <- 0
   block <- 0

   # SET UP A MATRIX TO STORE THE GENERATED DATA
   sim.matrix <- matrix(numeric(0), ncol=6)

   # SET LOOP COUNTER
   numloops <- 0

   # LOOP UNTIL THE SET NUMBER OF CASES AND OR CONTROLS IS ACHIEVED OR THE 
   # THE SET POPULATION SIZE TO SAMPLE FROM IS REACHED
   while(complete==0 && complete.absolute==0)
     {

       # GENERATE THE TRUE GENOTYPE DATA
       geno.data <- sim.geno.data(num.obs=n, geno.model=g.model, MAF=freq)
       allele.A <- geno.data$allele.A
       allele.B <- geno.data$allele.B
       geno <- geno.data$genotype
     
       # GENERATE THE TRUE ENVIRONMEANTAL EXPOSURE DATA
       env <- sim.env.data(num.obs=n, env.model=e.model, env.prev=e.prev, env.mean=e.mean, 
                                   env.sd=e.sd, env.low.lim=e.low.lim,env.up.lim=e.up.lim)
       
       # GENERATE THE TRUE INTERACTION DATA
       int <- geno*env       

       # GENERATE SUBJECT EFFECT DATA THAT REFLECTS BASELINE RISK: 
       # NORMALLY DISTRIBUTED RANDOM EFFECT VECTOR WITH APPROPRIATE 
       # VARIANCE ON SCALE OF LOG-ODDS
       s.effect.data <- sim.subject.data(n, b.OR)
					  
       # GENERATE THE TRUE OUTCOME DATA
       pheno.data <- sim.pheno.bin.GxE(num.obs=n, disease.prev=pheno.prev, genotype=geno, environment=env, 
                                       interaction=int, subject.effect.data=s.effect.data, geno.OR=g.OR, 
                                       env.OR=e.OR, int.OR=i.OR)
       true.phenotype <- pheno.data
       
       # GENERATE THE OBSERVED OUTCOME DATA FROM WHICH WE SELECT CASES AND CONTROLS
       obs.phenotype <- get.obs.pheno(phenotype=true.phenotype, pheno.model=0, pheno.error=ph.error)
       pheno <- obs.phenotype
       
       # STORE THE TRUE OUTCOME, GENETIC AND ENVIRONMENT AND ALLELE DATA IN AN OUTPUT MATRIX 
       # WHERE EACH ROW HOLDS THE RECORDS OF ONE INDIVUDAL
       sim.matrix.temp <- cbind(pheno,geno,allele.A,allele.B,env,int)

       # UPDATE THE MATRIX THAT HOLDS ALL THE DATA GENERATED SO FAR, AFTER EACH LOOP
       sim.matrix <- rbind(sim.matrix, sim.matrix.temp)

       # SELECT OUT CASES
       sim.matrix.cases <- sim.matrix[pheno==1,]

       # SELECT OUT CONTROLS
       sim.matrix.controls <- sim.matrix[pheno==0,]

       # COUNT THE NUMBER OF CASES AND CONTROLS THAT HAS BEEN GENERATED
       cases.simulated <- dim(sim.matrix.cases)[1]
       controls.simulated <- dim(sim.matrix.controls)[1]

       # TEST IF THERE ARE AT LEAST ENOUGH CASES ALREADY SIMULATED
       # IF THERE ARE, DEFINE THE CASE ELEMENT OF THE DATA MATRIX
       if(cases.simulated >= ncases)
       {
         sim.matrix.cases <- sim.matrix.cases[1:ncases,]
         cases.complete <- 1
       }

       # TEST IF THERE ARE AT LEAST ENOUGH CONTROLS ALREADY SIMULATED
       # IF THERE ARE, DEFINE THE CONTROL ELEMENT OF THE DATA MATRIX
       if(controls.simulated>=ncontrols)
       {
         sim.matrix.controls <- sim.matrix.controls[1:ncontrols,]
         controls.complete <- 1
       }

       # HAVE WE NOW GENERATED THE SET NUMBER OF CASES AND CONTROLS?
       complete <- cases.complete*controls.complete		

       # HAVE WE EXCEEDED THE TOTAL SAMPLE SIZE ALLOWED?
       complete.absolute <- (((block+1)*n)>=max.sample.size)
       if(complete.absolute==1) {sample.size.excess <- 1}else{sample.size.excess <- 0}

        # INCREMENT LOOP COUNTER
        numloops <- numloops + 1
       
   }

   # STACK FINAL DATA MATRIX WITH CASES FIRST
   sim.matrix <- rbind(sim.matrix.cases,sim.matrix.controls)
   totalnumrows <- dim(sim.matrix)[1]
   sim.matrix <- cbind(1:totalnumrows, sim.matrix)

   # NAME THE COLUMNS OF THE MATRIX AND RETURN IT AS A DATAFRAMEDATAFRAME
   colnames(sim.matrix) <- c("id", "phenotype", "genotype", "allele.A", "allele.B", "environment", "interaction")
   mm <- list(data=data.frame(sim.matrix), allowed.sample.size.exceeded=sample.size.excess)
}

