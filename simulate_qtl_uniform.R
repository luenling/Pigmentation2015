library(msm)
library(foreach)
library(parallel)
library(doMC)
library(reshape2)

registerDoMC(18)
#registerDoMC(4)
getDoParWorkers()
setwd("/Volumes/Temp/Lukas/Data/A7_pigmentation/Simulation")
gt.vals = c(2,1,0) ### value for each gt: aa:2, aA: 1, AA: 0

### get gt for locus
get.gt.locus <- function(gt.freq,gt.vals,total.ind){
  ### gets genotype frequencies and gt values for a locus and returns random genotypes for all ind. for this locus
  return(colSums(gt.vals*rmultinom(total.ind,1,gt.freq)))
}

### locus effects
loc.eff <- function(h.vals,g.type) { # takes the dominance and genotype vectors and returns a vector of gt effects (aa:1,aA:h,AA:0)
  return(g.type/2 + (g.type==1)*(h.vals-0.5))
}

### environmental variance of phenotype
v.E = 0.1
afs_cosm=c(0.2,0.83) # cosmopolitan allele frequencies
h2=0.3
get.V.E <- function(afs_loc,eff_loc,h_loc,h2=0.5,F.in=0,n=1000,eps.int=0.0){
  gt.freqs = sapply(afs_loc,function (p) c(p^2 + (p-p^2)*F.in,2*(p-p^2)*(1-F.in),(1-p)^2 + (p-p^2)*F.in))
  pos <- c()
  pos <- apply(gt.freqs,2,function (x) get.gt.locus(x,c(2,1,0),n))
  pos <- matrix(pos,nrow=length(afs_loc),ncol=n,byrow=T)
  ## pigmentation phenotypes resulting from genotypes
  ## additive effects (eff centered around 0)
  ind.eff <- apply(pos,2,function(x) loc.eff(h_loc,x)) 
  ## interaction component
  intxn <- ( eps.int*ind.eff[1,]*ind.eff[2,] - eps.int/2.0  )*mean(abs(eff_loc))
  ## additive effects (eff centered around 0)
  pheno <- eff_loc %*% ( ind.eff - 0.5 ) 
  pheno <- pheno + intxn
  return((1-h2)/h2*var(as.vector(pheno)))
}

simulate.eps <- function(v.E=0.0,total.ind =2100,n.replicates=3,n.pools=100,eps.int=0.0,
                       n.freqs=c(0.5,0.5),F.inbreed=0.0,eff=c(0.5,0.5),h.loci = c(0.5,0.5) ) {  
  #######  Initialize
  n.loci = length(n.freqs)
  ### genotype frequencies for each locus, aa, aA, AA
  gt.freqs = sapply(n.freqs,function (p) c(p^2 + (p-p^2)*F.inbreed,2*(p-p^2)*(1-F.inbreed),(1-p)^2 + (p-p^2)*F.inbreed)) 
  mean_eff=sum(gt.vals * gt.freqs %*% eff)
  d.replicates=total.ind/n.replicates # ind. per replicate
  ########results = foreach (i=1:n.runs,.combine="rbind") %dopar% {
  ## generate genotypes randomly (independent loci)
  # alleles in sample by binomial sampling from population. alleles
  # pos <- matrix(rbinom(n.loci*total.ind,2,p),nrow=n.loci,ncol=total.ind)    ### alleles A a, a:light allele, A: dark, gt.value: 2:aa 1:aA 0:AA
  pos <- c()
  #pos <- replicate(total.ind,get.gt(gt.freqs,gt.vals))
  #pos <- matrix(pos,nrow=n.loci,ncol=total.ind,byrow=T)
  pos <- apply(gt.freqs,2,function (x) get.gt.locus(x,gt.vals,total.ind))
  pos <- matrix(pos,nrow=n.loci,ncol=total.ind,byrow=T)
  # ind.eff only holds the dominance corrected effect of each locus, goes from 0 to 1 (or higher/lower in case of under/overdominance)
  ind.eff <- apply(pos,2,function(x) loc.eff(h.loci,x)) 
  ## pigmentation phenotypes resulting from genotypes
  ## interaction component, mean(abs(eff)) should always be total effect of both loci
  ## centered around zero roughly too
  #int <- (eps.int*ind.eff[1,]*ind.eff[2,] - eps.int/2 )*mean(abs(eff))
  ## as in the paper, not the symmetriced version 
  int <- (eps.int*(ind.eff[1,]-0.5)*(ind.eff[2,]-0.5))*mean(abs(eff))
  ## additive effects (ind.eff centered around 0)
  pheno <- eff %*% (ind.eff - 0.5)
  ## environmental effect
  pheno <- pheno + int + rnorm(total.ind,0,sqrt(v.E))
  ## get the randomized array of indices 
  indices=sample(1:total.ind)
  p.val <-  rep(0,n.loci)
  odds.ratio_ld <- rep(0,n.loci)
  odds.ratio_l <- rep(0,n.loci)
  odds.ratio_d <- rep(0,n.loci)
  # log.odds.ratio_ajb <- matrix(0,nrow=n.runs,ncol=n.loci)
  #     cors <- matrix(0,nrow=n.runs,ncol=4)
  af_avg <- rep(0,n.loci)
  af_light <- rep(0,n.loci)
  af_dark <- rep(0,n.loci)      
  dat = array(dim=c(n.loci,4,n.replicates)) # table for light dark cmh
  dat_l =  array(dim=c(n.loci,4,n.replicates)) # light base cmh
  dat_d =  array(dim=c(n.loci,4,n.replicates)) # dark base cmh
  alldark = c()
  alllight = c()
  ## get replicates
  for (j in 1:n.replicates) {
    ## indviduals in replicates
    sample.indices = indices[(1+(j-1)*d.replicates):(j*d.replicates)]
    ##  phenotypes in repl.  	
    sampled.pheno = pheno[sample.indices]
    ## order indices by pheno and take extreme ones
    ## light allele high gt value: reverse order as light have hihger phenotypic value!
    o.pheno <- sample.indices[order(sampled.pheno,decreasing=T)]
    light   <- o.pheno[1:n.pools]   # take "n.pools" lightest colored individuals for light pool
    dark    <- o.pheno[((d.replicates-n.pools)+1):d.replicates] # take "n.pools" darkest individuals for dark pool
    alllight=c(alllight, light)
    alldark=c(alldark, dark)
    # par(mfrow=c(1,1)); b.points=seq(min(pheno[sample.indices]),max(pheno[sample.indices]),by=0.5); hist(pheno[sample.indices],breaks=b.points); hist(pheno[light], col="yellow",add=T,breaks=b.points); hist(pheno[dark], col="blue",add=T,breaks=b.points)
    
    ## genotypes for the pools
    geno.light  <-  pos[,light]
    geno.dark   <-  pos[,dark]
    ## allele counts       
    ac.light=rowSums(geno.light)
    ac.dark =rowSums(geno.dark)
    ac.base =rowSums(pos[,sample.indices])
    ac.rep=rowSums(pos[,sample.indices])
    dat[,,j]=array(c(ac.light,ac.dark, n.pools*2-ac.light, n.pools*2-ac.dark), dim=c( n.loci,4))
    dat_l[,,j]=array(c(ac.light,ac.base, n.pools*2-ac.light,d.replicates*2-ac.base), dim=c( n.loci,4))
    dat_d[,,j]=array(c(ac.dark,ac.base, n.pools*2-ac.dark,d.replicates*2-ac.base), dim=c( n.loci,4))
  }
  af_avg=rowSums(pos)/(2*total.ind)
  ## calculate OR and AFs
  for (x in 1:n.loci){
    mhtable = array(dat[x,,],dim=c(2,2,n.replicates))
    cmh.res= mantelhaen.test(mhtable)
    odds.ratio_ld[x] = cmh.res$estimate
    p.val[x]=cmh.res$p.value
#     mhtable = array(dat_l[x,,],dim=c(2,2,n.replicates))
#     odds.ratio_l[x] = mantelhaen.test(mhtable)$estimate
#     mhtable = array(dat_d[x,,],dim=c(2,2,n.replicates))
#     odds.ratio_d[x] = mantelhaen.test(mhtable)$estimate
    af_light[x]=sum(dat[x,1,])/(2*n.replicates*n.pools)
    af_dark[x]=sum(dat[x,2,])/(2*n.replicates*n.pools)
  }
  results=c(p.val,odds.ratio_ld)
  names(results)=c("p1","p2","or1","or2")
#results=c(p.val,odds.ratio_ld,af_light,af_dark)
#names(results)=c("p1","p2","or1","or2","af1","af2","afl1","afl2","afd1","afd2")

  return(results)
}

draw.params <- function(n,t1,t2) {
  return((runif(n,min=-0.5,max=0.5)*t1)+t2)
}


### Runs for assessing parameter combinations
### Fixed afs
afs.eu=c(0.2,0.83)
afs.sa=c(0.17,0.47)
### Fis
###F.sa=0.25
F.sa=0.35
F.eu=0.0
### base effect
### light allele effect positive
eff1=1.0
eps.eff=c(1)
### loci
n.loci=2
### heritability
h2=0.30
afs.cosm=c(0.2,0.83)
### Experimental setup
n.pools.sa = 65 # size of light and dark pools
n.rep.sa = 3 #number of replicates
total.ind.sa =700*n.rep.sa # total individuals
n.pools.eu = 100 # size of light and dark pools
n.rep.eu = 6 #number of replicates
total.ind.eu =1500*n.rep.eu # total individuals

# set random seed
seed=as.integer(runif(1)*1000)
set.seed(seed)
print(paste("random seed set to",seed,"at",date()))
# number of runs
n.runs=500000
#n.runs=5000000
# limits:
h.lims=c(-0.125,1.125)
eff.lims=c(0.25,4.5)
#eps.lims=c(-3,3)
## for no eps.int
eps.lims=c(0,0)
# transformations
h.t1=h.lims[2]-h.lims[1]
eff.t1=eff.lims[2]-eff.lims[1]
eps.t1=eps.lims[2]-eps.lims[1]
h.t2=mean(h.lims)
eff.t2=mean(eff.lims)
eps.t2=mean(eps.lims)
t1=c(h.t1,h.t1,eff.t1,eps.t1)
t2=c(h.t2,h.t2,eff.t2,eps.t2)
params=draw.params(4,t1,t2)
names(params)=c("h1","h2","eff","eps.int")

# #a=t(replicate(1000,draw.params(4,t1,t2)))
# #summary(a)
# print(paste("starting noepes",n.runs,"sims at",date()))
# simres.noeps <- foreach (i=1:n.runs,.combine="rbind") %dopar% {
#   # draw random parameters
#   if (i%%100000 == 0) print(paste(i,"iterations performed at",date()))
#   params=draw.params(4,t1,t2)
#   names(params)=c("h1","h2","eff","eps.int")
#   # assign parameters
#   eff=c(eff1/(1+abs(params["eff"])),eff1*params["eff"]/(1+abs(params["eff"])))
#   h.loci=params[c("h1","h2")]
#   # get variance
#   v.E=get.V.E(afs.cosm,eff,h.loci,h2,eps.int=params["eps.int"])
#   # do europe
#   res.eu <- simulate.eps(v.E=v.E,total.ind=total.ind.eu,n.replicates=n.rep.eu,n.pools=n.pools.eu,n.freqs=afs.eu,
#                       F.inbreed=F.eu,eff=eff,h.loci=h.loci,eps.int=params["eps.int"] )  
#   # do SA
#   res.sa <- simulate.eps(v.E=v.E,total.ind=total.ind.sa,n.replicates=n.rep.sa,n.pools=n.pools.sa,n.freqs=afs.sa,
#                       F.inbreed=F.sa,eff=eff,h.loci=h.loci,eps.int=params["eps.int"] )
#   return(c(params,v.E,res.eu,res.sa))
# }  
# #colnames(simulation_results)=c(names(params),"v.E","p1Eu","p2Eu","or1Eu","or2Eu","af1Eu","af2Eu","afl1Eu","afl2Eu","afd1Eu","afd2Eu","p1Sa","p2Sa","or1Sa","or2Sa","af1Sa","af2Sa","afl1Sa","afl2Sa","afd1Sa","afd2Sa")
# colnames(simres.noeps)=c(names(params),"v.E","P.eu.tan","P.eu.bab","or.eu.tan","or.eu.bab","P.sa.tan","P.sa.tan","or.sa.tan","or.sa.bab")
# 
# #save(simulation_results,file="simres.unif.eps.RData.gz",compress=T)
# saveRDS(simres.noeps,file="simres.unif.noeps.2mio.rds",compress=TRUE)
# #load("simres.unif.eps.RData.gz")

#rm(simres.noeps)

## for eps.int
eps.lims=c(-1.5,3.5)
eps.t1=eps.lims[2]-eps.lims[1]
eps.t2=mean(eps.lims)
t1=c(h.t1,h.t1,eff.t1,eps.t1)
t2=c(h.t2,h.t2,eff.t2,eps.t2)
params=draw.params(4,t1,t2)
names(params)=c("h1","h2","eff","eps.int")

print(paste("starting eps",n.runs,"sims at",date()))
simres.eps <- foreach (i=1:n.runs,.combine="rbind") %dopar% {
  # draw random parameters
  if (i%%100000 == 0) print(paste(i,"iterations performed at",date()))  
  params=draw.params(4,t1,t2)
  names(params)=c("h1","h2","eff","eps.int")
  # assign parameters
  eff=c(eff1/(1+abs(params["eff"])),eff1*params["eff"]/(1+abs(params["eff"])))
  h.loci=params[c("h1","h2")]
  # get variance
  v.E=get.V.E(afs.cosm,eff,h.loci,h2,eps.int=params["eps.int"])
  # do europe
  res.eu <- simulate.eps(v.E=v.E,total.ind=total.ind.eu,n.replicates=n.rep.eu,n.pools=n.pools.eu,n.freqs=afs.eu,
                         F.inbreed=F.eu,eff=eff,h.loci=h.loci,eps.int=params["eps.int"] )  
  # do SA
  res.sa <- simulate.eps(v.E=v.E,total.ind=total.ind.sa,n.replicates=n.rep.sa,n.pools=n.pools.sa,n.freqs=afs.sa,
                         F.inbreed=F.sa,eff=eff,h.loci=h.loci,eps.int=params["eps.int"] )
  return(c(params,v.E,res.eu,res.sa))
}  
#colnames(simulation_results)=c(names(params),"v.E","p1Eu","p2Eu","or1Eu","or2Eu","af1Eu","af2Eu","afl1Eu","afl2Eu","afd1Eu","afd2Eu","p1Sa","p2Sa","or1Sa","or2Sa","af1Sa","af2Sa","afl1Sa","afl2Sa","afd1Sa","afd2Sa")
colnames(simres.eps)=c(names(params),"v.E","P.eu.tan","P.eu.bab","or.eu.tan","or.eu.bab","P.sa.tan","P.sa.tan","or.sa.tan","or.sa.bab")
print(paste("writing output file at",date()))
save(simres.eps,file="simres.unif.500K.special.eps.2.RData")
#save(simres.eps,file="simres.unif.eps.2.25.500K.RData",compress=T)















