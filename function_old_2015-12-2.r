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


##gene drop simulation given the number of childern in every generation
# family_strct.2g2c <- data.frame(family=c(1,1,1,1), person=c(1,2,3,4), father=c(0,0,1,1), 
#                                    mother=c(0,0,2,2), sex=c(1,2,1,1), affect=c(1,2,2,2)) #1=male, 2=female, 1=unaffected, 2=affected
#use the file format as in Merlin
#simulate the tranmission vector by generation and determine the affected status
#this generate the whole family first and then keep only those matched input
gene_family <- function(family_strct=family_strct.2g2c, n_family=1000, beta0=110, var=5, haplotype.risk=haplotype.risk) {
  n_family_member <- length(family_strct$person)
  data_family <- matrix(NA, nrow=n_family*n_family_member, ncol=(6+n_snp*2))
  tran_vec <- matrix(NA, nrow=n_family*n_family_member, ncol=3)
  #the strategy is generate each individual one by one and check affected status if fail then restart from the first individual
  data_family.idx <- 1 #index of the current family member being generated
  n_family.idx <- 1 #index of how many families have been generated
  #   affect_idx <- which(family_strct$affect!=0) #who is affected status is not missing
  while(n_family.idx <= n_family) {
    disease_vec <- matrix(NA, nrow=n_family_member, ncol=1) #store generated affected status
    family.haplo <- matrix(NA, nrow=n_family_member, ncol=2) #store haplotype info.
    ind.idx <- 1 # which individual is generating
    while(ind.idx <= n_family_member) { #until the generated matches the input
      #generate the current individual
      if(family_strct$father[ind.idx]==0 & family_strct$mother[ind.idx]==0) { #if founder
        haplo.id <- sample(1:n_haplo, 2, replace=T)
        disease_prob <- sum(haplotype.risk[haplo.id])
        disease_vec[ind.idx] <- beta0+disease_prob+rnorm(1,0,sd=sqrt(var))
      }
      else{ #if not founder
        haplo.id <- c(sample(family.haplo[family_strct$father[ind.idx],], 1), sample(family.haplo[family_strct$mother[ind.idx],], 1))
        disease_prob <- sum(haplotype.risk[haplo.id])
        disease_vec[ind.idx] <- beta0+disease_prob+rnorm(1,0,sd=sqrt(var))
      }
      family.haplo[ind.idx,] <- haplo.id
      ind.idx <- ind.idx + 1 #move to the next
    }
    #save the haplotype file
    letter.idx <- 1 #indicator used in the transmission vector
    for(i in 1:n_family_member) {
      #store transmission vector
      data_family[data_family.idx, ] <- unlist(c(n_family.idx, family_strct[i,2:5], disease_vec[i], haplotype[family.haplo[i,1],-c(1:2)], haplotype[family.haplo[i,2],-c(1:2)]))
      if(family_strct$father[i]==0 & family_strct$mother[i]==0) { #if founder
        tran_vec[data_family.idx,] <- c(n_family.idx, LETTERS[letter.idx], LETTERS[letter.idx+1])
        letter.idx <- letter.idx + 2
      }
      else{ #if not founder then compare with his/her parents
        current_row <- (n_family.idx-1)*n_family_member
        #print(current_row)
        tran_vec[data_family.idx,] <- c(n_family.idx, ifelse(family.haplo[i,1] == family.haplo[family_strct$father[i],1], tran_vec[family_strct$father[i]+current_row, 2], tran_vec[family_strct$father[i]+current_row, 3]) 
                                        ,ifelse(family.haplo[i,2] == family.haplo[family_strct$mother[i],1], tran_vec[family_strct$mother[i]+current_row, 2], tran_vec[family_strct$mother[i]+current_row, 3])) 
      }
      data_family.idx <- data_family.idx + 1
    }
    #     print(n_family.idx)
    n_family.idx <- n_family.idx + 1
  }
  colnames(data_family) <- c("family","person","father","mother","sex","affect",rep(paste("SNP", 1:n_snp, sep=""),2))
  colnames(tran_vec) <- c("family","h1","h2")
  return(list(data_family=data.frame(data_family, stringsAsFactors=F), tran_vec=data.frame(tran_vec, stringsAsFactors=F)))
}
# family_generated <- gene_family()    

##regular family test trap
trapq <- function(data=family_generated_2g2c, f=risk.variant.id, weight=F) {
  ##hypothesis testing count the number of carrier haplotypes transmitted to the offsprings
  n_family <- max(data$data_family$family)
  n_family_member <- table(data$data_family$family)
  #check if founder's haplotype carries any variant's with f < 0.1
  if(length(f)==1 & f[1] <1) { #support allele frequency or list of snps
    snp2look.idx <-  which(snp$FREQ1 < f) # snp to look for
  } else(snp2look.idx <-  f)
  
  
  #start looking at each family
  test.stat <- sapply(1:n_family, function(x) {
    family.idx=x 
    #     print(family.idx)
    current_row=sum(n_family_member[1:family.idx]) - n_family_member[family.idx]
    h1 <- data$data_family[(current_row+1):(current_row+n_family_member[family.idx]),7:(6+n_snp)]
    h2 <- data$data_family[(current_row+1):(current_row+n_family_member[family.idx]),-c(1:(6+n_snp))]
    #adjust here for more founders
    affect <- data$data_family[(current_row+1):(current_row+n_family_member[family.idx]),"affect"]
    tran_vec <- data$tran_vec[(current_row+1):(current_row+n_family_member[family.idx]),]
    family_strct <- data$data_family[(current_row+1):(current_row+n_family_member[family.idx]),1:6]
    person <- c(0,family_strct[,2]) #adding 0 to check if the chromosome if from missing founder
    
    #define who is a founder
    founder <- rep(0, length=n_family_member[family.idx])
    for(i in 1:n_family_member[family.idx]) {
      founder[i] <- ifelse((family_strct$father[i]==0 & family_strct$mother[i]==0), 1,  #full founder
                           ifelse(!(family_strct$father[i] %in% person), 0.5, #half founder from father
                                  ifelse(!(family_strct$mother[i] %in% person), -0.5, 0))) #half founder from mother
    }
    founder_idx <- which(founder!=0)
    n_founder <- sum(abs(founder))
    carrier <- as.vector(unlist(sapply(founder_idx, function(y){
      if(founder[y]==1) {
        c(ifelse(sum(h1[y, snp2look.idx]==1)>0, tran_vec[y, "h1"], NA),
          ifelse(sum(h2[y, snp2look.idx]==1)>0, tran_vec[y, "h2"], NA))
      }else{
        if(founder[y]==0.5) { #from father
          ifelse(sum(h1[y, snp2look.idx]==1)>0, tran_vec[y, "h1"], NA)
        }else{
          if(founder[y]==-0.5) { #from mother
            ifelse(sum(h2[y, snp2look.idx]==1)>0, tran_vec[y, "h2"], NA)
          }
        }
      }  
    })))
    #     print(carrier)
    
    #     carrier <- c( #indicator of which haplotypes is carrier
    #       ifelse(sum(h1[1, snp2look.idx]==1)>0, tran_vec[1, "h1"], NA) #check first founder's h1
    #       ,ifelse(sum(h2[1, snp2look.idx]==1)>0, tran_vec[1, "h2"], NA) #check first founder's h2
    #       ,ifelse(sum(h1[2, snp2look.idx]==1)>0, tran_vec[2, "h1"], NA) #check second founder's h1
    #       ,ifelse(sum(h2[2, snp2look.idx]==1)>0, tran_vec[2, "h2"], NA) #check second founder's h2
    #     )
    n_carrier <- sum(!is.na(carrier)) #no. of carrier haplotypes
    
    if(!(n_carrier==(2*n_founder) | n_carrier==0)) { #skip families with 0 or 4 carrier haplotypes in founders
      summary_stat = 0
      n_carrier_each <- apply(tran_vec[ , c("h1","h2")], 1, function(x) sum(x %in% carrier))
      for(i in 1:n_family_member[family.idx]) {
        for(j in i:n_family_member[family.idx]) {
          summary_stat <- summary_stat + abs(n_carrier_each[i] - n_carrier_each[j])*abs(affect[i]-affect[j])
        }
      }
      #calculate expectation and variance
      founder <- t(combn(LETTERS[1:(2*n_founder)], n_carrier))
      S <- apply(founder, 1, function(x) {
        carrier <- x #founder's haplotype
        #       print(carrier)
        summary_stat_observed = 0
        n_carrier_each <- apply(tran_vec[ , c("h1","h2")], 1, function(x) sum(x %in% carrier))
        for(i in 1:n_family_member[family.idx]) {
          for(j in i:n_family_member[family.idx]) {
            summary_stat_observed <- summary_stat_observed + abs(n_carrier_each[i] - n_carrier_each[j])*abs(affect[i]-affect[j])
          }
        }
        summary_stat_observed
      })
      
      c(observed=summary_stat, mean=mean(S), var=sum((S-mean(S))^2)/nrow(founder), n_carrier=n_carrier, family.idx=family.idx)
    }
  }
  )
  test.stat <- data.frame(do.call(rbind, test.stat))
  
  if(weight==F){ #regular TRAP
    v <- test.stat$var
    se <- sqrt(sum(v))
    #e_avg <- mean(e) #overall expectation
    final.test.stat <- sum(test.stat$observed - test.stat$mean)/se
    #c(t, n_test.data)#T and number of informative families
    p.value <- 2*pnorm(abs(final.test.stat), lower=F) #p-value
    #   print(c(final.test.stat=final.test.stat, sum_e=sum(test.stat$mean), se=se, mean_observed=mean(test.stat$observed), mean_mean=mean(test.stat$mean), mean_var=mean(test.stat$var), p.value=p.value, n_info_family=sum(test.stat$var!=0)))
    list(final.test.stat=final.test.stat, sum_e=sum(test.stat$mean), se=se, mean_observed=mean(test.stat$observed), mean_mean=mean(test.stat$mean), mean_var=mean(test.stat$var), p.value=p.value, n_info_family=sum(test.stat$var!=0))
  }else{
    #weighted test still applying Lyapunoc's CLT
    idx.informative <- 
    v <- test.stat$var
    se <- sqrt(nrow(test.stat))
    #e_avg <- mean(e) #overall expectation
    final.test.stat <- sum((test.stat$observed - test.stat$mean)/sqrt(v))/se
    #c(t, n_test.data)#T and number of informative families
    p.value <- 2*pnorm(abs(final.test.stat), lower=F) #p-value
    #   print(c(final.test.stat=final.test.stat, sum_e=sum(test.stat$mean), se=se, mean_observed=mean(test.stat$observed), mean_mean=mean(test.stat$mean), mean_var=mean(test.stat$var), p.value=p.value, n_info_family=sum(test.stat$var!=0)))
    list(final.test.stat=final.test.stat, sum_e=sum(test.stat$mean), se=se, mean_observed=mean(test.stat$observed), mean_mean=mean(test.stat$mean), mean_var=mean(test.stat$var), p.value=p.value, n_info_family=sum(test.stat$var!=0))
  }
}
# trapq()
