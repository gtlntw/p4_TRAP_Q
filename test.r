
result <- replicate(20, {
  family_generated <<- gene_family_pe(family_strct=family_strct_ped, n_family=n_family, trait_mean=trait_mean, Beta=Beta, dis_cutoff = dis_cutoff, exact_affected = exact_affected) 
  #add pseudo null variant since FB-SKAT can only work with 2+ in a gene
  family_generated_diploid <- hap2dip(data=family_generated, risk.variant.id=risk.variant.id, save.file=F)
  trap <- family.test(data=family_generated, f=risk.variant.id, summary.stat="5")$p.value
  trap_center <- family.score.test(data=family_generated, f=risk.variant.id, summary.stat="5", lm.test=F)$p.value
  trap_within_center <- family.score.test(data=family_generated, f=risk.variant.id, summary.stat="5", lm.test=F, within_center = T)$p.value
  trap_center_comb <- family.score.test(data=family_generated, f=risk.variant.id, summary.stat="5", lm.test=T)$p.value.comb
  trap_within_center_comb <- family.score.test(data=family_generated, f=risk.variant.id, summary.stat="5", lm.test=T, within_center = T)$p.value.comb
  pedgene <- pedgene(ped=family_generated_diploid$ped, geno=family_generated_diploid$geno)$pgdf$pval.burden
  print(c(trap, trap_center, trap_within_center, trap_center_comb, trap_within_center_comb, pedgene))
  c(trap, trap_center, trap_within_center, trap_center_comb, trap_within_center_comb, pedgene)
})
result
