##load project related functions
source('function.r')

##setup function
parseCommandArgs <- function (evaluate = TRUE) 
{
  arglist <- list()
  args <- commandArgs()
  i <- which(args == "--args")
  if (length(i) == 0 || length(args) < 1) 
    return(invisible())
  args <- args[(i + 1):length(args)]
  for (i in seq(1, length(args), by = 2)) {
    value <- NA
    tryCatch(value <- as.double(args[i + 1]), warning = function(e) {
    })
    if (is.na(value)) {
      value <- args[i + 1]
      if (substr(value, 1, 2) == "c(") 
        value <- eval(parse(text = args[i + 1]))
    }
    if (evaluate) 
      assign(args[i], value, inherits = TRUE)
    arglist[[length(arglist) + 1]] <- value
    names(arglist)[length(arglist)] <- args[i]
  }
  return(arglist)
}

#need to take in parameter
## passed in from the command prompt.
parseCommandArgs()
## Will disply the default values on the first run,
print(seed)
print(n_rep)
print(beta0)
print(beta1)
print(var_error)
print(n_family) #number of family
print(risk.variant)

#####################################################
## ... your simualtion codek
##the founder's haplotypes
#import the haplotypes generated by cosi
haplotype <- read.table("out_100k_10k_1kb.hap-1", header=F)
colnames(haplotype) <- c("HAP", "CHROM", paste("SNP", 1:(ncol(haplotype)-2), sep=""))
snp <-read.table("out_100k_10k_1kb.pos-1", header=T)
#make allele 1 is the minor allele, 2 is the common allele
temp.idx <- snp$FREQ1 > snp$FREQ2
temp.freq <- snp$FREQ2
snp$FREQ2[temp.idx] <- snp$FREQ1[temp.idx]
snp$FREQ1[temp.idx] <- temp.freq[temp.idx]
#also change the genotype file
haplotype[,which(temp.idx==T)+2] <- 3 - haplotype[,which(temp.idx==T)+2]

#allele frequency
nrow(snp) #total number of snp
sum(snp$FREQ1 < 0.05) # number of snp with f < 0.05
sum(snp$FREQ1 < 0.01) # number of snp with f < 0.01
sum(snp$FREQ1 == 0.0001) # number of singletons

##assign risk variants and the corresponding effect size (proportional to allele frequency)
# null <- FALSE
n_haplo <- 10000
n_snp <- ncol(haplotype)-2

#set up causual SNPs
#generate risk haplotypes
# risk.variant.id <- c(3, 8,19,21,23,27,44,47,49,50)
risk.variant.id <- c(risk.variant) #2 for common 7 for rare 39 for super rare
risk.haplo.id <- which(apply(2-as.matrix(haplotype[, risk.variant.id+2]), 1, sum)>0)
(risk.haplo.f <- mean(apply(2-as.matrix(haplotype[, risk.variant.id+2]), 1, sum)>0)) #carrier haplotype frequency
haplotype.risk <- rep(0, length=nrow(haplotype))
#assign mean relative risk calculate the haplotype variants p(A|h)
haplotype.risk[risk.haplo.id] <- beta1
mean(haplotype.risk[risk.haplo.id]) #mean relative risk

##gene drop simulation for two generations
family_strct.2g2c <- data.frame(family=c(1,1,1,1), person=c(1,2,3,4), father=c(0,0,1,1), 
                                mother=c(0,0,2,2), sex=c(1,2,1,1), affect=c(1,2,2,2)) #1=male, 2=female, 1=unaffected, 2=affected


rep.idx <<- 1
sim_result <- replicate(n_rep, {
  print(rep.idx)
  rep.idx <<- rep.idx + 1

  #trapq
  family_generated_2g2c <<- gene_family(family_strct=family_strct.2g2c, beta0 = beta0, var=var_error, n=n_family, haplotype.risk=haplotype.risk) 
    ##test based on new imputation framework
  result.trapq <- trapq(data=family_generated_2g2c, f=risk.variant.id)$p.value
#   result.trapq.weight <- trapq(data=family_generated_2g2c, f=risk.variant.id, weight=T)$p.value
  
  #population deseign
  generated_pop <- gene.data.pop(beta0=beta0, beta1=beta1, var=var_error, f=snp$FREQ1[risk.variant.id], n=n_family*4)
  result.pop <- lm(pheno~I(H1+H2), data=generated_pop)
  coef(summary(result.pop))
  (result.pop <- coef(summary(result.pop))[2,4])
  
  #only report p.value
  c(result.trapq, result.pop)
})
## Write out your results to a csv file
result.df <- as.data.frame(t(sim_result))
colnames(result.df) <- c("result.trapq", "result.pop")
result.df <- cbind(seed,beta0,beta1,var_error,risk.variant,risk.haplo.f,n_family,result.df)
write.csv(result.df, paste("res_",beta0,"_",beta1,"_",var_error,"_",n_family,"_",seed,".csv",sep=""), row.names=FALSE)
## R.miSniam
## End(Not run)
## Not run:

####################commands to run jobs in parallel
##passed in from the command prompt.
parseCommandArgs()
## Will disply the default values on the first run,
##initialize
seed=1000
#number of replications
n_rep=500
#number of family
n_family=1000
beta0=110
var_error=5
# print(null) #null=FALSE
#family structure
# print(family_strct) #family_strct='family_strct.2g3c'

print(getwd())

##command line
parallel <- function(...) {
  names.args <- names(list(...))
  value.args <- as.vector(list(...))
  ##commands to run
  for(i in 1:4) {
    rfile="mainSim.R"
    cmd <- paste("R --vanilla --args seed ", seed, " ", paste(names.args, value.args, collapse=" "),
                 " < ", rfile, " > ", "mainSim_", paste(value.args, collapse="_"), ".Rout", seed, " 2>&1", sep="")
    print(cmd)
    writeLines(cmd, fileConn) #write jobs to text
    
    #add seed number
    seed<<-seed+1
  }
}
##clean-up & initialization
system("rm *.csv")
system("rm *.Rout*")
system("rm *.out*")
fileConn<-file("jobs.txt", "w")

##create your jobs here
for(i in seq(0,1, length.out = 11)) {
  parallel(n_rep=n_rep, beta0=beta0, beta1=i, var_error=var_error, n_family=n_family, risk.variant=2) #common
}
for(i in seq(0,2, length.out = 11)) {
  parallel(n_rep=n_rep, beta0=beta0, beta1=i, var_error=var_error, n_family=n_family, risk.variant=7) #rare
}
for(i in seq(0,3, length.out = 11)) {
  parallel(n_rep=n_rep, beta0=beta0, beta1=i, var_error=var_error, n_family=n_family, risk.variant=39) #super rare
}

##write jobs.txt
close(fileConn)
##submit the jobs.txt using runslurm.pl
system("runslurm.pl -replace -logfile jobs.log -copts \"--time=12:00:00 --mem=256\" jobs.txt")
##submit jobs.txt via runbatch.pl i.e. mosic mode
# dir=paste("/net/frodo", substr(getwd(), 0, 500), sep="")
# cmd <- paste("runbatch.pl -replace -logfile jobs.log -concurrent 40 -copts \"-E", dir," -j2,4-8\" jobs.txt &", sep="")
# print(cmd)
# system(cmd)
###################################################
# The color-blind friednly palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#check power using true, imputed founder carrier, minimum offspring carrier
#three versions -- command and rare and super rare #2 for common, 7 for rare, 39 for super rare
result <- read.csv("2015-07-09 TRAP_Q.csv", header=T)
result <- result %>% melt(id.vars=c("seed", "beta0", "beta1", "var_error","risk.variant","risk.haplo.f","n_family"), variable.name="method", value.name="p.value")
result.plot <- result %>% group_by(risk.variant, risk.haplo.f, beta1, method) %>% 
  summarise(n=n(), power=mean(p.value<0.05, na.rm=T))
result.plot$risk.haplo.f <- factor(result.plot$risk.haplo.f, labels=c("f=0.0039","f=0.0178","f=0.202"))

#TRAP imputation
pd <- position_dodge(0.0)
filter(result.plot) %>% ggplot(aes(x=beta1, y=power, ymax=max(power), group=method, col=method)) +
  #   geom_point(size=3, alpha=1) +
  facet_wrap(~risk.haplo.f, ncol=3, scale="free_x") +
  geom_line(size=1.2, alpha=0.7, position=pd) +
  geom_point(size=1.2, position=pd) +
#   ggtitle("f=0.202, consider effect size of risk haplotypes, TRAP") +
#   ggtitle("f=0.0178, consider effect size of risk haplotypes, TRAP") +
  ggtitle("var=5") +
  labs(x="beta_1") +
  scale_y_continuous(limits=c(0,1)) +
  theme_gray(base_size = 20) +
  theme(legend.position="right", panel.background = element_rect(fill = 'grey85'))# +
#   scale_color_manual(values=cbbPalette)