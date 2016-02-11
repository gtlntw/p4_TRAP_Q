source("function.r")

#need to take in parameter
## passed in from the command prompt.
parseCommandArgs()
## Will disply the default values on the first run,
print(beta1)
print(f)
print(n_rep)

## ... your simualtion code
#joint test
dist <- function(sib_pheno_sum, sib_pheno_diff, n_dep_carrier, n_carrier_diff, beta0, beta1, var) {
  log(dnorm(sib_pheno_sum, mean=3*beta0+beta1*n_dep_carrier, sd=sqrt(3*var)))+log(dnorm(sib_pheno_diff, mean=beta1*n_carrier_diff,sd=sqrt(6*var)))
} 
fn_decomp <- function(beta, data=generated_family.3sib) {
  beta0 <- beta[1]
  beta1 <- beta[2]
  var <- beta[3]
  with(data, -sum(dist(sib_pheno_sum, sib_pheno_diff, n_dep_carrier-H1_f-H2_f-H1_m-H2_m, n_carrier_diff, beta0=beta0, beta1=beta1, var=var)))
}
fn_alt <- function(beta, data=generated_family.3sib) {
  beta0 <- beta[1]
  beta1 <- beta[2]
  var <- beta[3]
  with(data, -sum(dist(sib_pheno_sum, sib_pheno_diff, n_dep_carrier, n_carrier_diff, beta0=beta0, beta1=beta1, var=var)))
}
fn_null <- function(beta, data=generated_family.3sib) {
  beta0 <- beta[1]
  beta1 <- 0
  var <- beta[2]
  with(data, -sum(dist(sib_pheno_sum, sib_pheno_diff, n_dep_carrier, n_carrier_diff, beta0=beta0, beta1=beta1, var=var)))
}

sim <- function(beta1=1, f=f, n_rep=1000){
  #parameter
  beta0=110; beta1=beta1; var=5; f=f; n_family=666; n_pop=3*n_family
  n_rep=n_rep
  #simulation to see the distribution of beta1 and beta3(interaction) and chisquare test
  n_replicate <- 1
  sim.result <- replicate(n_rep, {
    if(n_replicate%%10==0) print(n_replicate)
    n_replicate <<- n_replicate + 1
    generated_family.3sib <- gene.data.family.3sib(beta0=beta0, beta1=beta1, var=var, f=f, n_family=n_family)
    #association test
    #sib_sum test
    M.sum <- summary(lm(sib_pheno_sum~n_dep_carrier, data=generated_family.3sib))
    (M.sum.coef <- coef(M.sum))
    (M.sum.p.value <- M.sum.coef[2,4])
    #sib_diff test
    M.diff <- summary(lm(sib_pheno_diff~n_carrier_diff, data=generated_family.3sib))
    (M.diff.coef <- coef(M.diff))
    (M.diff.p.value <- M.diff.coef[2,4])
    
    #joint tset
    fit_alt <- optim(beta<-c(110,1,5), fn_alt, data=generated_family.3sib, hessian=T, method="L-BFGS-B",lower=c(-Inf,-Inf, 0))
    wald.test <- fit_alt$par[2]/sqrt(1/fit_alt$hessian[2,2])
    (M.joint.wald.p.value <- pchisq(wald.test^2, 1, lower=F))
    #   print(M.joint.wald.p.value)
    fit_null <- optim(beta<-c(110,5), fn_null, data=generated_family.3sib, method="L-BFGS-B",lower=c(-Inf, 0))
    (M.joint.p.value <- pchisq(2*(fit_null$value - fit_alt$value), 1, lower=F))
    #   print(M.joint.p.value)
    
    #QTDT type decomposition
    fit_decomp <- optim(beta<-c(110,1,5), fn_decomp, data=generated_family.3sib, hessian=T, method="L-BFGS-B",,lower=c(-Inf,-Inf, 0))
    decomp.wald.test <- fit_decomp$par[2]/sqrt(1/fit_decomp$hessian[2,2])
    (M.joint.decom.wald.p.value <- pchisq(decomp.wald.test^2, 1, lower=F))
    
    #simulation in population to see power
    generated_pop <- gene.data.pop(beta0=beta0, beta1=beta1, var=var, f=f, n=n_pop)
    result.pop <- lm(pheno~I(H1+H2), data=generated_pop)
    coef(summary(result.pop))
    (p.value.pop <- coef(summary(result.pop))[2,4])
    
    #   c(M.full.beta=M.full.beta, n_info=nrow(generated_family_informative), M.full.beta.p.value=M.full.beta.p.value, M.14.p.value.family=M.14.p.value.family, M.15.p.value.family=M.15.p.value.family, M.45.p.value.family=M.45.p.value.family, M.145.p.value.family=M.145.p.value.family, p.value.pop=p.value.pop)
    c(M.sum.p.value=M.sum.p.value, M.diff.p.value=M.diff.p.value, M.joint.p.value=M.joint.p.value,
      M.joint.wald.p.value=M.joint.wald.p.value, p.value.pop=p.value.pop,
      M.joint.decom.wald.p.value=M.joint.decom.wald.p.value)
  })
  sim.result
}

sim_result <- sim(beta1=beta1, f=f, n_rep=n_rep)

## Write out your results to a csv file
write.csv(data.frame(seed=seed, beta1=beta1, f=f, M.sum.p.value=sim_result[1,], M.diff.p.value=sim_result[2,], M.joint.p.value=sim_result[3,], 
                     M.joint.wald.p.value=sim_result[4,], p.value.pop=sim_result[5,], M.joint.decom.wald.p.value=sim_result[6,]),
          paste("res_",beta1,"_",f,"_",seed,".csv",sep=""), row.names=FALSE)


# sim.result
sim.result <- read.csv("2014-12-02 idea expansion - vary f.csv")
alpha=0.05
# apply(sim.result[c(1:4),], 1, function(x) mean(x)) #mean of betas
# apply(sim.result[c(1:4),], 1, function(x) sd(x)) #se of betas
sim.result.01 <- subset(sim.result, f==0.01) 
sim.result.05 <- subset(sim.result, f==0.05) 
sim.result.20 <- subset(sim.result, f==0.20) 
apply(sim.result.01[,-c(1:3)], 2, function(x) tapply(x, list(sim.result.01$beta1), function(y) mean(y<alpha)))
apply(sim.result.05[,-c(1:3)], 2, function(x) tapply(x, list(sim.result.05$beta1), function(y) mean(y<alpha)))
apply(sim.result.20[,-c(1:3)], 2, function(x) tapply(x, list(sim.result.20$beta1), function(y) mean(y<alpha)))
# power.omnibus.test <- mean(sim.result[1,] < alpha/2 | sim.result[2,] < alpha/2) 
# power.result <- c(power.result[1:2], power.omnibus.test, power.result[-c(1:2,6)])
barplot(apply(sim.result.01[,-c(1:3)], 2, function(x) tapply(x, list(sim.result.01$beta1),
        function(y) mean(y<alpha))), beside=T, names=c("sib_sum", "sib_diff","joint-LRT","joint-Wald","Population", "joint-decomp-Wald"), 
        legend=c(expression(beta[1]=="0"), expression(beta[1]=="0.4"), expression(beta[1]=="0.8")), 
        ylim=c(0,1), ylab="Power", main="f = 0.01")

barplot(apply(sim.result.05[,-c(1:3)], 2, function(x) tapply(x, list(sim.result.05$beta1),
                                                             function(y) mean(y<alpha))), beside=T, names=c("sib_sum", "sib_diff","joint-LRT","joint-Wald","Population", "joint-decomp-Wald"), 
        legend=c(expression(beta[1]=="0"), expression(beta[1]=="0.2"), expression(beta[1]=="0.4")), 
        ylim=c(0,1), ylab="Power", main="f = 0.05")

barplot(apply(sim.result.20[,-c(1:3)], 2, function(x) tapply(x, list(sim.result.20$beta1),
                                                             function(y) mean(y<alpha))), beside=T, names=c("sib_sum", "sib_diff","joint-LRT","joint-Wald","Population", "joint-decomp-Wald"), 
        legend=c(expression(beta[1]=="0"), expression(beta[1]=="0.1"), expression(beta[1]=="0.2")), 
        ylim=c(0,1), ylab="Power", main="f = 0.20")


# xx <- barplot(power.result, names=c("sib_sum", "sib_diff","joint-LRT","joint-Wald","Population", "joint-decomp-Wald"), ylim=c(0,1), ylab="power")
# text(x = xx, y = power.result, label=round(power.result,3), pos = 3, col = "red")
