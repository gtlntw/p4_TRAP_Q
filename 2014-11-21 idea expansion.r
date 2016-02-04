source("function.r")
#joint test
dist <- function(sib_pheno_sum, sib_pheno_diff, n_dep_carrier, n_carrier_diff, beta0, beta1, var) {
  log(dnorm(sib_pheno_sum, mean=2*beta0+beta1*n_dep_carrier, sd=sqrt(2*var)))+log(dnorm(sib_pheno_diff, mean=beta1*n_carrier_diff,sd=sqrt(2*var)))
} 
fn_alt <- function(beta, data=generated_family) {
  beta0 <- beta[1]
  beta1 <- beta[2]
  var <- beta[3]
  with(data, -sum(dist(sib_pheno_sum, sib_pheno_diff, n_dep_carrier, n_carrier_diff, beta0=beta0, beta1=beta1, var=var)))
}
fn_null <- function(beta, data=generated_family) {
  beta0 <- beta[1]
  beta1 <- 0
  var <- beta[2]
  with(data, -sum(dist(sib_pheno_sum, sib_pheno_diff, n_dep_carrier, n_carrier_diff, beta0=beta0, beta1=beta1, var=var)))
}

#parameter
beta0=110; beta1=0; var=5; f=0.01; n_family=1000; n_pop=2*n_family
n_rep=2500

#generate familyi
generated_family <- gene.data.family(beta0=beta0, beta1=beta1, var=var, f=f, n_family=n_family, structure="2g2c")

#number of S=0,1,2 siblings
S0 <- subset(generated_family, S==0)
S1 <- subset(generated_family, S==1)
S2 <- subset(generated_family, S==2)
table(generated_family$S)
prop.table(table(generated_family$S))
with(generated_family, tapply(sib_pheno_sum, list(S), mean))
with(generated_family, tapply(sib_pheno_sum, list(S,n_dep_carrier), mean))
with(generated_family, tapply(sib_pheno_sum, list(S,n_IBD_carrier), mean))
with(generated_family, tapply(sib_pheno_sum, list(S,n_nonIBD_carrier), mean))
with(generated_family, tapply(sib_pheno_sum, list(S,n_carrier_diff), mean))
with(generated_family, tapply(sib_pheno_diff, list(S), mean))
with(generated_family, tapply(sib_pheno_diff, list(S,n_indep_carrier), mean))
with(generated_family, tapply(sib_pheno_diff, list(S,n_IBD_carrier), mean))
with(generated_family, tapply(sib_pheno_diff, list(S,n_nonIBD_carrier), mean))
with(generated_family, tapply(sib_pheno_diff, list(S,n_carrier_diff), mean))
with(generated_family, table(n_indep_carrier))
with(generated_family, table(n_carrier_diff))
with(generated_family, plot(sib_pheno_sum,sib_pheno_diff))
with(generated_family, cor(sib_pheno_sum,sib_pheno_diff))

#simulation to see the distribution of beta1 and beta3(interaction) and chisquare test
n_replicate <- 1
sim.result <- replicate(n_rep, {
  if(n_replicate%%10==0) print(n_replicate)
  n_replicate <<- n_replicate + 1
  generated_family <- gene.data.family(beta0=beta0, beta1=beta1, var=var, f=f, n_family=n_family, structure="2g2c")
  #association test
  #sib_sum test
  M.sum <- summary(lm(sib_pheno_sum~n_dep_carrier, data=generated_family))
  (M.sum.coef <- coef(M.sum))
  (M.sum.p.value <- M.sum.coef[2,4])
  #sib_diff test
  M.diff <- summary(lm(sib_pheno_diff~n_carrier_diff, data=generated_family))
  (M.diff.coef <- coef(M.diff))
  (M.diff.p.value <- M.diff.coef[2,4])
  
  #joint tset
  fit_alt <- optim(beta<-c(110,1,5), fn_alt, data=generated_family, hessian=T, method="L-BFGS-B",lower=c(-Inf,-Inf, 0))
  wald.test <- fit_alt$par[2]/sqrt(1/fit_alt$hessian[2,2])
  (M.joint.wald.p.value <- pchisq(wald.test^2, 1, lower=F))
#   print(M.joint.wald.p.value)
  fit_null <- optim(beta<-c(110,5), fn_null, data=generated_family, method="L-BFGS-B",lower=c(-Inf,0))
  (M.joint.p.value <- pchisq(2*(fit_null$value - fit_alt$value), 1, lower=F))
#   print(M.joint.p.value)
  
  
  #simulation in population to see power
  generated_pop <- gene.data.pop(beta0=beta0, beta1=beta1, var=var, f=f, n=n_pop)
  result.pop <- lm(pheno~I(H1+H2), data=generated_pop)
  coef(summary(result.pop))
  (p.value.pop <- coef(summary(result.pop))[2,4])
  
#   c(M.full.beta=M.full.beta, n_info=nrow(generated_family_informative), M.full.beta.p.value=M.full.beta.p.value, M.14.p.value.family=M.14.p.value.family, M.15.p.value.family=M.15.p.value.family, M.45.p.value.family=M.45.p.value.family, M.145.p.value.family=M.145.p.value.family, p.value.pop=p.value.pop)
 c(M.sum.p.value=M.sum.p.value, M.diff.p.value=M.diff.p.value, M.joint.p.value=M.joint.p.value, M.joint.wald.p.value=M.joint.wald.p.value, p.value.pop=p.value.pop)
})
# sim.result
alpha=0.05
# apply(sim.result[c(1:4),], 1, function(x) mean(x)) #mean of betas
# apply(sim.result[c(1:4),], 1, function(x) sd(x)) #se of betas
(power.result <- apply(sim.result, 1, function(x) mean(x<alpha, na.rm=T)))
power.omnibus.test <- mean(sim.result[1,] < alpha/2 | sim.result[2,] < alpha/2) 
power.result <- c(power.result[1:2], power.omnibus.test, power.result[-c(1:2)])
xx <- barplot(power.result, names=c("sib_sum", "sib_diff", "Omnibus","joint-LRT", "joint-Wald","Population"), ylim=c(0,1), ylab="power")
text(x = xx, y = power.result, label=round(power.result,3), pos = 3, col = "red")
