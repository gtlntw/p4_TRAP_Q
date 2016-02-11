beta0=110; beta1=0; var=5; f=0.20; n_family=5000; n_pop=2*n_family;
generated_family <- gene.data.family(beta0=beta0, beta1=beta1, var=var, f=f, n_family=n_family, structure="2g2c")

idx=1
sim.result <- replicate(1000, {
  print(idx)
  idx<<-idx+1
  generated_family <- gene.data.family(beta0=beta0, beta1=beta1, var=var, f=f, n_family=n_family, structure="2g2c")
  generated_family <- subset(generated_family, S==2)
  test <- data.frame(
  pheno=c(generated_family$sib1_pheno, generated_family$sib2_pheno),
  allele=c(generated_family$H1+generated_family$H2, generated_family$H1_s+generated_family$H2_s)
  )
  with(test, sum(allele)/(nrow(generated_family)*4))
  #joint tset
  #sib_sum test
  M.sum <- summary(lm(pheno~allele, data=test))
  (M.sum.coef <- coef(M.sum))
  M.sum.p.value <- M.sum.coef[2,4]
  
  
  #population
  generated_pop <- gene.data.pop(beta0=beta0, beta1=beta1, var=var, f=f, n=n_pop)
  result.pop <- lm(pheno~I(H1+H2), data=generated_pop)
  coef(summary(result.pop))
  (p.value.pop <- coef(summary(result.pop))[2,4])
  
  c(M.sum.p.value, p.value.pop)
}
)
apply(sim.result, 1, function(x) mean(x< .05))
