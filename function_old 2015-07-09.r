gene.data.family <- function(beta0=110, beta1=5, var=5, f=0.01, n_family=10, structure="2g2c") {
  if(!(structure %in% c("2g2c", "3g", "2g3c"))) stop("need to specify the family structure")
  
  ##generate sibpairs
  if(structure=="2g2c") {
    variable_names <- c("H1_f", "H2_f", "f_pheno",
      "H1_m", "H2_m", "m_pheno",
      "H1", "H2", "sib1_pheno",
      "H1_s", "H2_s", "sib2_pheno",
      "S","S0", "S1", "S2", "S11", "S12",
      "n_IBD_chrm", "n_nonIBD_chrm", "n_IBD_carrier", "n_nonIBD_carrier","n_indep_carrier","n_dep_carrier","n_carrier_diff",
      "sib_pheno_diff","sib_pheno_sum")
    data_family <- array(NA, c(n_family, length(variable_names)), dimnames=list(NULL,variable_names))
    n_family_count <- 1
    while(n_family_count <= n_family) {
      #father
      H1_f <- rbinom(1,1,p=f) #if the first haplotype carries risk raviant
      H2_f <- rbinom(1,1,p=f) #if the second haplotype carries risk raviant
      f_pheno <- beta0+beta1*(H1_f+H2_f)+rnorm(1,0,sd=sqrt(var))
      #mother
      H1_m <- rbinom(1,1,p=f) #if the first haplotype carries risk raviant
      H2_m <- rbinom(1,1,p=f) #if the second haplotype carries risk raviant
      m_pheno <- beta0+beta1*(H1_m+H2_m)+rnorm(1,0,sd=sqrt(var))
      #first sib
      H1_inh <- sample(c(1,2), 1, replace=TRUE) 
      H1 <- ifelse(H1_inh==1, H1_f, H2_f)
      H2_inh <- sample(c(1,2), 1, replace=TRUE) 
      H2 <- ifelse(H2_inh==1, H1_m, H2_m)
      sib1_pheno <- beta0+beta1*(H1+H2)+rnorm(1,0,sd=sqrt(var))
      #second sib
      H1_inh_s <- sample(c(1,2),1, replace=TRUE) 
      H1_s <- ifelse(H1_inh_s==1, H1_f, H2_f)
      H2_inh_s <- sample(c(1,2),1, replace=TRUE) 
      H2_s <- ifelse(H2_inh_s==1, H1_m, H2_m)
      sib2_pheno <- beta0+beta1*(H1_s+H2_s)+rnorm(1,0,sd=sqrt(var))
      
      S0 <- (H1_inh_s != H1_inh) & (H2_inh_s != H2_inh)
      S2 <- (H1_inh_s == H1_inh) & (H2_inh_s == H2_inh)
      S1 <- ((H1_inh_s == H1_inh) & (H2_inh_s != H2_inh)) | ((H1_inh_s != H1_inh) & (H2_inh_s == H2_inh))
      S11 <- ((H1_inh_s == H1_inh) & (H2_inh_s != H2_inh))
      S12 <- ((H1_inh_s != H1_inh) & (H2_inh_s == H2_inh))
      
      S <- (S0==TRUE)*0 + (S1==TRUE)*1 + (S2==TRUE)*2
      
      #modified to unaffected father, affected mother and both affected children
      #print(n_family_count)
      #count the T statistics
      n_IBD_chrm <- S1 + 2*S2
      n_nonIBD_chrm <- 4*S0 + 2*S1
      n_IBD_carrier <- S2*(H1+H2) + ifelse(S11==T, H1, 0) + ifelse(S12==T, H2, 0)
      n_nonIBD_carrier <- S0*(H1+H2+H1_s+H2_s) + ifelse(S11==T, H2+H2_s, 0) + ifelse(S12==T, H1+H1_s, 0)
      n_indep_carrier <- n_IBD_carrier + n_nonIBD_carrier
      n_dep_carrier <- 2*n_IBD_carrier + n_nonIBD_carrier
      n_carrier_diff <- ((H1+H2)-(H1_s+H2_s))
      n_nonIBD_carrier <- S0*(H1+H2+H1_s+H2_s) + ifelse(S11==T, H2+H2_s, 0) + ifelse(S12==T, H1+H1_s, 0)
      sib_pheno_diff = (sib1_pheno - sib2_pheno)
      sib_pheno_sum = (sib1_pheno+sib2_pheno)
      
      data_family[n_family_count,] <- c(H1_f, H2_f, f_pheno,
                                        H1_m, H2_m, m_pheno,
                                        H1, H2, sib1_pheno,
                                        H1_s, H2_s, sib2_pheno,
                                        S,S0, S1, S2, S11, S12,
                                        n_IBD_chrm, n_nonIBD_chrm, n_IBD_carrier, n_nonIBD_carrier,n_indep_carrier,n_dep_carrier,n_carrier_diff,
                                        sib_pheno_diff,sib_pheno_sum)
      n_family_count <- n_family_count + 1
    }
    as.data.frame(data_family)
  }
}
generated_family <- gene.data.family(beta0=110, beta1=5, var=5, f=0.01, n_family=10, structure="2g2c")
generated_family

gene.data.family.3sib <- function(beta0=110, beta1=5, var=5, f=0.01, n_family=10, structure="2g3c") {
  if(!(structure %in% c("2g2c", "3g", "2g3c"))) stop("need to specify the family structure")
  
  ##generate sibpairs
  if(structure=="2g3c") {
    variable_names <- c("H1_f", "H2_f", "f_pheno",
                        "H1_m", "H2_m", "m_pheno",
                        "H1", "H2", "sib1_pheno",
                        "H1_s", "H2_s", "sib2_pheno",
                        "H1_s2", "H2_s2", "sib3_pheno",
                        "n_dep_carrier","n_carrier_diff",
                        "sib_pheno_diff","sib_pheno_sum")
    data_family <- array(NA, c(n_family, length(variable_names)), dimnames=list(NULL,variable_names))
    n_family_count <- 1
    while(n_family_count <= n_family) {
      #father
      H1_f <- rbinom(1,1,p=f) #if the first haplotype carries risk raviant
      H2_f <- rbinom(1,1,p=f) #if the second haplotype carries risk raviant
      f_pheno <- beta0+beta1*(H1_f+H2_f)+rnorm(1,0,sd=sqrt(var))
      #mother
      H1_m <- rbinom(1,1,p=f) #if the first haplotype carries risk raviant
      H2_m <- rbinom(1,1,p=f) #if the second haplotype carries risk raviant
      m_pheno <- beta0+beta1*(H1_m+H2_m)+rnorm(1,0,sd=sqrt(var))
      #first sib
      H1_inh <- sample(c(1,2), 1, replace=TRUE) 
      H1 <- ifelse(H1_inh==1, H1_f, H2_f)
      H2_inh <- sample(c(1,2), 1, replace=TRUE) 
      H2 <- ifelse(H2_inh==1, H1_m, H2_m)
      sib1_pheno <- beta0+beta1*(H1+H2)+rnorm(1,0,sd=sqrt(var))
      #second sib
      H1_inh_s <- sample(c(1,2),1, replace=TRUE) 
      H1_s <- ifelse(H1_inh_s==1, H1_f, H2_f)
      H2_inh_s <- sample(c(1,2),1, replace=TRUE) 
      H2_s <- ifelse(H2_inh_s==1, H1_m, H2_m)
      sib2_pheno <- beta0+beta1*(H1_s+H2_s)+rnorm(1,0,sd=sqrt(var))
      #third sib
      H1_inh_s2 <- sample(c(1,2),1, replace=TRUE) 
      H1_s2 <- ifelse(H1_inh_s2==1, H1_f, H2_f)
      H2_inh_s2 <- sample(c(1,2),1, replace=TRUE) 
      H2_s2 <- ifelse(H2_inh_s2==1, H1_m, H2_m)
      sib3_pheno <- beta0+beta1*(H1_s2+H2_s2)+rnorm(1,0,sd=sqrt(var))
      
      #modified to unaffected father, affected mother and both affected children
      #print(n_family_count)
      #count the T statistics
      n_dep_carrier <- (H1+H2) + (H1_s+H2_s) + (H1_s2+H2_s2)
      n_carrier_diff <- ((H1+H2) - (H1_s+H2_s)) + ((H1+H2) - (H1_s2+H2_s2))
      sib_pheno_diff = (sib1_pheno - sib2_pheno) + (sib1_pheno - sib3_pheno)
      sib_pheno_sum = (sib1_pheno+sib2_pheno+sib3_pheno)
      
      data_family[n_family_count,] <- c(H1_f, H2_f, f_pheno,
                                        H1_m, H2_m, m_pheno,
                                        H1, H2, sib1_pheno,
                                        H1_s, H2_s, sib2_pheno,
                                        H1_s2, H2_s2, sib3_pheno,
                                        n_dep_carrier,n_carrier_diff,
                                        sib_pheno_diff,sib_pheno_sum)
      n_family_count <- n_family_count + 1
    }
    as.data.frame(data_family)
  }
}
generated_family.3sib <- gene.data.family.3sib(beta0=110, beta1=5, var=5, f=0.01, n_family=10, structure="2g3c")
generated_family.3sib



gene.data.pop <- function(beta0=110, beta1=5, var=5, f=0.01, n=10) {
  variable_names <- c("H1", "H2", "pheno")
  data_pop <- array(NA, c(n, length(variable_names)), dimnames=list(NULL,variable_names))
  n_count <- 1
  while(n_count <= n) {
    #father
    H1 <- rbinom(1,1,p=f) #if the first haplotype carries risk raviant
    H2 <- rbinom(1,1,p=f) #if the second haplotype carries risk raviant
    pheno <- beta0+beta1*(H1+H2)+rnorm(1,0,sd=sqrt(var))
    data_pop[n_count,] <- c(H1, H2, pheno)
    n_count <- n_count + 1
  }
  as.data.frame(data_pop)
}

generated_pop <- gene.data.pop(beta0=110, beta1=5, var=5, f=0.01, n=1000)
generated_pop

