DataPath <- "./data_"
rep  <- 1
cat("Replicate number=",rep,"\n\n")

library(AlphaSimR)

# ensure reproducible using seed
set.seed(202600 + rep)

## Initial parameters
nt     <-  2  # number of traits
ngener <-  20 # number of simulation generations/years

h2   <- c(0.5, 0.3) # Trait heritabilities
varP <- c(100,100)  # Phenotypic variances
varE <- (1-h2)*varP # Residual variances
varG <- h2*varP     # Genetic variances
rg   <-  0.3        # genetic correlation
re   <-  0.1        # residual correlation

cat("N traits=",nt,"\n")
cat("h2      =",h2,"\n")
cat("Var(P)  =",varP,"\n")
cat("Var(E)  =",varE,"\n")
cat("Var(G)  =",varG,"\n")
cat("gen cor =",rg,"\n")
cat("res cor =",re,"\n")

genCor <- matrix(c( 1, rg,
                   rg,  1), ncol=nt) # Genetic correlation matrix

resVar <- matrix(c(varE[1], re*sqrt(varE[1]*varE[2]), # Residual covariance matrix
                   re*sqrt(varE[1]*varE[2]), varE[2]), ncol=nt)

w  = c(1, 1)     # Trait economic weights
cat("Eco weights=",w,"\n")

nsel_female <- 2000       # number of selected females
nsel_male   <- 200        # number of selected males
nbase       <- 2*nsel_female # base population size, has about equal number of both females and males
nFounders   <- 280        # founder population size

cat("\nBase size       :", nbase,"\n")
cat("Selected males  :", nsel_male,"\n")
cat("Selected females:", nsel_female,"\n")

Nchr   <- 30              # number of chromosomes
nsites <- 3600            # number of segregation sites
jmp    <- 2               # SNP jump to segregation sites
Qjmp   <- 120             # QTL jump to segregation sites
Nsnp   <- Nchr*nsites/jmp # number of SNP markers
NsnpPerChr <- Nsnp/Nchr   # number of SNP markers per chromosome
NQTL       <- Nchr*nsites/Qjmp # number of QTLs
NQTLperChr <- max(NQTL/Nchr,1) # number of QTLs per chromosome

cat("\nN chromosomes   :", Nchr,"\n")
cat("N segreg. sites   :", nsites,"\n")
cat("N SNP             :", Nsnp,"\n")
cat("N SNP per chr     :", NsnpPerChr,"\n")
cat("N QTL             :", NQTL,"\n")
cat("N QTL per chr     :", NQTLperChr,"\n")

# Create founder haplotypes
tmp = paste0(DataPath, rep, "/founderPopRep.RData")       
if (!file.exists(tmp)) {                                     
  cat("\nrunMacs...\n")
  founderPop = runMacs(nInd = nFounders,        # Number of individuals
                     nChr = Nchr,        # Number of chromosomes
                     segSites = nsites,  # Number of segregating sites to keep per chromosome
                     species = "CATTLE") # The CATTLE history
  save(founderPop, file=tmp)             # save the founder population
} else {
  cat("\nload founder population...\n")
  load(tmp)
  print(tmp)
}


cat("SimParam...\n")
SP = SimParam$
  new(founderPop)$
  addTraitA(NQTLperChr, mean=c(0,0), var=varG, corA=genCor,  # number of QTL per chr; genetic var&correlation
	    gamma=TRUE, shape=0.4)$ # gamma density with shape leads to large QTLs
  setVarE(varE=resVar)$    # residual covariance matrix  
  setSexes("yes_sys")$    # assign systematically balanced sexes
  addSnpChip(NsnpPerChr, minSnpFreq=0.05)$  # Randomly assigns eligible SNPs to a SNP chip
  setTrackPed(TRUE)        # pedigree tracking ON


snp_map <- getSnpMap(simParam = SP)

qtl_map <- getQtlMap(trait = 1, simParam = SP)


# take a look on chromosome 1
snp_chr1 <- subset(snp_map, chr == 1)
qtl_chr1 <- subset(qtl_map, chr == 1)
head(snp_chr1)
head(qtl_chr1)

# Counts
nrow(snp_chr1)  # 1800 SNP markers
nrow(qtl_chr1)  # 30 QTL

# Check whether SNP and QTL positions overlap
intersect(snp_map$id, qtl_map$id)   # character(0)

# write out
tmp = paste0(DataPath, rep, "/SNP_map.txt")
write.table(snp_map, file = tmp, quote = FALSE, row.names = FALSE, col.names = FALSE)

tmp = paste0(DataPath, rep, "/QTL_map.txt")
write.table(qtl_map, file = tmp, quote = FALSE, row.names = FALSE, col.names = FALSE)


cat("newPop...\n")
pop = newPop(founderPop, simParam=SP) # Create new population

# Check simulated phenotypes and genetic values
cat("\n=== Initial founder population ===\n")

cat("\nFirst six phenotype records:\n")
print(head(pheno(pop)))

cat("\nPhenotypic variance-covariance matrix:\n")
print(var(pheno(pop)))

cat("\nGenetic variance-covariance matrix:\n")
print(var(gv(pop)))




# Generate base population
pop_base = randCross(pop = pop, nCrosses = nbase)  # Generate 4000 offspring from 4000 random crosses


cat("\n=== Base offspring population ===\n")

print(var(pheno(pop_base)))
print(var(gv(pop_base)))





pop_all <- c(pop, pop_base) 
pop     <- pop_base       
cat("Pop so far, N indiv(pop_all)=", nInd(pop_all),"\n")
cat("Breeding pop, N indiv(pop)  =", nInd(pop),"\n")

# Calculation average true economic value in the base population
# mean(gv(pop)%*%w)

gmean   <- meanG(pop) # Returns the mean genetic values for all traits
gMean_1 <- gmean[1]
gMean_2 <- gmean[2]

id_now <- as.numeric(slot(pop_all,"id"))
nn     <- length(id_now)
gener  <- matrix(1,nn,1)
id_num <- matrix(id_now, nn,1)

for(year in 2:(ngener+1)){
  cat("Year...", year,"N ind=",nInd(pop_all),"\n")
  pop_cross   <- randCross(pop, nCrosses=nsel_female) # Make random crosses
  id_now      <- as.numeric(slot(pop_cross,"id"))
  nn          <- length(id_now)
  taulu       <- matrix(year,nn,1)
  gener       <- c(gener, taulu)
  id_num      <- c(id_num, id_now)
  gmean       <- meanG(pop_cross)                     # Mean genetics of the crosses
  # sukup       <- slot(pop_cross, "sex")
  gMean_1     = c(gMean_1, gmean[1])                  # Add the mean by trait
  gMean_2     = c(gMean_2, gmean[2])
  
  pop_tot     <- c(pop, pop_cross)                    # current population and new crosses
  pop_all     <- c(pop_all, pop_cross)                # complete population
  pop_males   <- selectInd(pop_tot, nInd=nsel_male,   sex="M", trait=selIndex, b=w) # male select
  pop_females <- selectInd(pop_tot, nInd=nsel_female, sex="F", trait=selIndex, b=w) # female select
  pop         <- c(pop_males, pop_females)            # new current population
}

cat("Final     , N indiv(pop_all)=", nInd(pop_all),"\n")
cat("Writing files...\n")
sp             <- slot(pop_all, "sex") # get sex of individuals in the pedigree
n              <- length(sp)
sukup          <- matrix(1,n,1)
sukup[sp=="F"] <- 2

pedi    <- as.data.frame(cbind(getPed(pop_all), gener, sukup))# pedigree of population individuals
id      <- as.matrix(pedi$id, n,1)
bvs     <- as.data.frame(cbind(id, bv(pop_all)))     # breeding values of individuals
geno    <- as.data.frame(cbind(id,pullSnpGeno(pop_all, simParam = SP))) # SNP genotype data
pheno   <- as.data.frame(cbind(id,gener, sukup, pheno(pop_all)))   # phenotypes
QTLgeno <- as.data.frame(cbind(id,pullQtlGeno(pop_all, simParam = SP))) # QTLs

tmp = paste0(DataPath, rep, "/sim.ped")
write.table(pedi,    file=tmp,   quote=FALSE, col.names = FALSE, row.names = FALSE)
tmp = paste0(DataPath, rep, "/TBV.txt")
write.table(bvs,     file=tmp,   quote=FALSE, col.names = FALSE, row.names = FALSE)
tmp = paste0(DataPath, rep, "/pheno.txt")
write.table(pheno,   file=tmp, quote=FALSE, col.names = FALSE, row.names = FALSE)
tmp = paste0(DataPath, rep, "/QTL.txt")
write.table(QTLgeno, file=tmp,   quote=FALSE, col.names = FALSE, row.names = FALSE)

# get the QTL values
h      <- SP$traits
ade1   <- slot(h[[1]],"addEff")
ade2   <- slot(h[[2]],"addEff")
QTLeff <- cbind(ade1, ade2)
tmp = paste0(DataPath, rep, "/QTLeff.txt")
write.table(QTLeff, file=tmp, quote=FALSE, col.names = FALSE, row.names = FALSE)

tmp = paste0(DataPath, rep, "/SNP.txt")
write.table(geno,    file=tmp,   quote=FALSE, col.names = FALSE, row.names = FALSE)

cat("Done.\n")

