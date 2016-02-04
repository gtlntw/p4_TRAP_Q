source("function.r")

#parameter
beta0=110; beta1=5; var=5; f=0.01; n_family=50000; n_pop=2*n_family

#generate familyi
generated_family <- gene.data.family(beta0=beta0, beta1=beta1, var=var, f=f, n_family=n_family, structure="2g2c")
#number of S=0,1,2 siblings
table(generated_family$S)
prop.table(table(generated_family$S))

#difference in trait between two siblings
sub1 <- subset(generated_family, generated_family$sib_pheno_diff<=5)
sub2 <- subset(generated_family, 5<generated_family$sib_pheno_diff | generated_family$sib_pheno_diff<10)
sub3 <- subset(generated_family, 10<generated_family$sib_pheno_diff)

result.analysis <- function(data) {
  attach(data)
  result <- list((table(S)), (prop.table(table(S))),c(sum(n_IBD_carrier), sum(n_nonIBD_carrier)), c(sum(n_IBD_chrm)-sum(n_IBD_carrier), sum(n_nonIBD_chrm)-sum(n_nonIBD_carrier)), sum(n_IBD_carrier)*(sum(n_nonIBD_chrm)-sum(n_nonIBD_carrier))/sum(n_nonIBD_carrier)/(sum(n_IBD_chrm)-sum(n_IBD_carrier)),prop.test(c(sum(n_IBD_carrier), sum(n_nonIBD_carrier)), c(sum(n_IBD_chrm), sum(n_nonIBD_chrm))))
  detach()
  result
}

result.analysis(sub1)
result.analysis(sub2)
result.analysis(sub3)

#simulation to see the distribution of beta1 and beta3(interaction) and chisquare test
generated_family <- gene.data.family(beta0=beta0, beta1=beta1, var=var, f=f, n_family=n_family, structure="2g2c")
#association test
test.data <-apply(generated_family, 1, function(x) {
  n_chrm <- x["n_IBD_chrm"] + x["n_nonIBD_chrm"]
  carrier <- c(rep(c(1,0), c(x["n_IBD_carrier"], x["n_IBD_chrm"]-x["n_IBD_carrier"])), rep(c(1,0), c(x["n_nonIBD_carrier"], x["n_nonIBD_chrm"]-x["n_nonIBD_carrier"])))
  IBD <- c(rep(1, x["n_IBD_chrm"]), rep(0, x["n_nonIBD_chrm"]))
  sib_pheno_diff <- x["sib_pheno_diff"]
  sib_pheno_sum <- x["sib_pheno_sum"]
  cbind(carrier, IBD, sib_pheno_diff=x["sib_pheno_diff"], sib_pheno_sum=x["sib_pheno_sum"])
})
test.data <- as.data.frame(do.call(rbind, test.data))
# sum_weight <- mean(test.data$sib_pheno_sum)
# weight <- test.data$sib_pheno_sum/sum_weight

M0 <- glm(IBD~sib_pheno_diff, family=binomial(link="logit"), data=test.data)
summary(M0)
M1 <- glm(IBD~carrier+sib_pheno_diff+carrier*sib_pheno_diff, family=binomial(link="logit"), data=test.data)
summary(M1)
anova(M0, M1, test = "Chisq")
# M2 <- glm(IBD~carrier+sib_pheno_diff+sib_pheno_sum, family=binomial(link="logit"), data=test.data)
# summary(M2)

#simulation in population to see power
generated_pop <- gene.data.pop(beta0=beta0, beta1=beta1, var=var, f=f, n=n_pop)
result <- lm(pheno~I(H1+H2), data=generated_pop)
summary(result)
