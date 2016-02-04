source("function.r")

#parameter
beta0=110; beta1=1; var=5; f=0.01; n_family=1000; n_pop=2*n_family
n_rep=1000

#generate familyi
generated_family <- gene.data.family(beta0=beta0, beta1=beta1, var=var, f=f, n_family=n_family, structure="2g2c")
#number of S=0,1,2 siblings
table(generated_family$S)
prop.table(table(generated_family$S))

#simulation to see the distribution of beta1 and beta3(interaction) and chisquare test
n_replicate <- 1
sim.result <- replicate(n_rep, {
  if(n_replicate%%10==0) print(n_replicate)
  n_replicate <<- n_replicate + 1
  generated_family <- gene.data.family(beta0=beta0, beta1=beta1, var=var, f=f, n_family=n_family, structure="2g2c")
  generated_family_informative <- generated_family[which(apply(generated_family[,c("H1", "H2", "H1_s", "H2_s")], 1, sum)>0),]
  test.data <-apply(generated_family_informative, 1, function(x) {
  #association test
#   test.data <-apply(generated_family, 1, function(x) {
    n_chrm <- x["n_IBD_chrm"] + x["n_nonIBD_chrm"]
    carrier <- c(rep(c(1,0), c(x["n_IBD_carrier"], x["n_IBD_chrm"]-x["n_IBD_carrier"])), rep(c(1,0), c(x["n_nonIBD_carrier"], x["n_nonIBD_chrm"]-x["n_nonIBD_carrier"])))
    IBD <- c(rep(1, x["n_IBD_chrm"]), rep(0, x["n_nonIBD_chrm"]))
    sib_pheno_diff <- x["sib_pheno_diff"]
    sib_pheno_sum <- x["sib_pheno_sum"]
    cbind(carrier, IBD, sib_pheno_diff=x["sib_pheno_diff"], sib_pheno_sum=x["sib_pheno_sum"])
  })
  test.data <- as.data.frame(do.call(rbind, test.data))
  max.sum <- max(test.data$sib_pheno_sum)
  # sum_weight <- mean(test.data$sib_pheno_sum)
  # weight <- test.data$sib_pheno_sum/sum_weight
  #test
  M.full <- glm(IBD~carrier+sib_pheno_diff+I(sib_pheno_sum-max.sum)+carrier*sib_pheno_diff+carrier*I(sib_pheno_sum-max.sum), family=binomial(link="logit"), data=test.data)
  (M.full.coef <- coef(summary(M.full)))
  M.full.beta <- M.full.coef[c(2,5,6),1]
  M.full.beta.p.value <- M.full.coef[c(2,5,6),4]

  M.14 <- glm(IBD~sib_pheno_diff+I(sib_pheno_sum-max.sum)+carrier*I(sib_pheno_sum-max.sum), family=binomial(link="logit"), data=test.data)
  M.14.result.family <- anova(M.14, M.full, test = "LRT")
  (M.14.p.value.family <- M.14.result.family$Pr[2])

  M.15 <- glm(IBD~sib_pheno_diff+I(sib_pheno_sum-max.sum)+carrier*sib_pheno_diff, family=binomial(link="logit"), data=test.data)
  M.15.result.family <- anova(M.15, M.full, test = "LRT")
  (M.15.p.value.family <- M.15.result.family$Pr[2])
  
  M.45 <- glm(IBD~carrier+sib_pheno_diff+I(sib_pheno_sum-max.sum), family=binomial(link="logit"), data=test.data)
  M.45.result.family <- anova(M.45, M.full, test = "LRT")
  (M.45.p.value.family <- M.45.result.family$Pr[2])
  
  M.145 <- glm(IBD~sib_pheno_diff+I(sib_pheno_sum-max.sum), family=binomial(link="logit"), data=test.data)
  M.145.result.family <- anova(M.145, M.full, test = "LRT")
  (M.145.p.value.family <- M.145.result.family$Pr[2])
  
  
  # M2 <- glm(IBD~carrier+sib_pheno_diff+sib_pheno_sum, family=binomial(link="logit"), data=test.data)
  # summary(M2)
  
  #simulation in population to see power
  generated_pop <- gene.data.pop(beta0=beta0, beta1=beta1, var=var, f=f, n=n_pop)
  result.pop <- lm(pheno~I(H1+H2), data=generated_pop)
  coef(summary(result.pop))
  (p.value.pop <- coef(summary(result.pop))[2,4])
  
#   c(M.full.beta=M.full.beta, n_info=nrow(generated_family_informative), M.full.beta.p.value=M.full.beta.p.value, M.14.p.value.family=M.14.p.value.family, M.15.p.value.family=M.15.p.value.family, M.45.p.value.family=M.45.p.value.family, M.145.p.value.family=M.145.p.value.family, p.value.pop=p.value.pop)
c(M.full.beta=M.full.beta, n_info=nrow(generated_family_informative), M.full.beta.p.value=M.full.beta.p.value, M.45.p.value.family=M.45.p.value.family, M.145.p.value.family=M.145.p.value.family, p.value.pop=p.value.pop)
})
sim.result
alpha=0.05
apply(sim.result[c(1:4),], 1, function(x) mean(x)) #mean of betas
apply(sim.result[c(1:4),], 1, function(x) sd(x)) #se of betas
(power.result <- apply(sim.result[-c(1:4),], 1, function(x) mean(x<alpha, na.rm=T)))
hist(power.result, )
xx <- barplot(power.result, names=c("beta1","beta4","beta5","beta4+5","beta1+4+5","Population"), ylim=c(0,1), ylab="power")
text(x = xx, y = power.result, label=round(power.result,3), pos = 3, col = "red")

