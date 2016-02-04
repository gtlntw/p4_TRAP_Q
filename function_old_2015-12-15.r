##parse parameters from command
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
select_family_strct <- function(family_strct) {
	switch(family_strct, #1=male, 2=female, 1=unaffected, 2=affected
				 "2g.2f.2c"=data.frame(family=c(1,1,1,1), person=c(1,2,3,4), father=c(0,0,1,1), 
                             mother=c(0,0,2,2), sex=c(1,2,1,1), affect=c(0,0,0,0)),
				 "2g.2a.2u"=data.frame(family=c(1,1,1,1), person=c(1,2,3,4), father=c(0,0,1,1), 
                             mother=c(0,0,2,2), sex=c(1,2,1,1), affect=c(1,2,1,2)),
				 "2g.3a.1u"=data.frame(family=c(1,1,1,1), person=c(1,2,3,4), father=c(0,0,1,1), 
                             mother=c(0,0,2,2), sex=c(1,2,1,1), affect=c(1,2,2,2)),
				 "2g.3a.2u"=data.frame(family=c(1,1,1,1,1), person=c(1,2,3,4,5), father=c(0,0,1,1,1), 
                             mother=c(0,0,2,2,2), sex=c(1,2,1,1,1), affect=c(1,1,2,2,2)),
				 "2g.4a.1u"=data.frame(family=c(1,1,1,1,1), person=c(1,2,3,4,5), father=c(0,0,1,1,1), 
                                mother=c(0,0,2,2,2), sex=c(1,2,1,1,1), affect=c(1,2,2,2,2)),
				 "3g.3a.4u"=data.frame(family=c(1,1,1,1,1,1,1), person=c(1,2,3,4,5,6,7), father=c(0,0,1,1,0,4,4), 
                             mother=c(0,0,2,2,0,5,5), sex=c(1,2,1,1,2,1,1), affect=c(1,2,1,2,1,2,1)),
				 "3g.2a.5u"=data.frame(family=c(1,1,1,1,1,1,1), person=c(1,2,3,4,5,6,7), father=c(0,0,1,1,0,4,4), 
                                 mother=c(0,0,2,2,0,5,5), sex=c(1,2,1,1,2,1,1), affect=c(1,2,1,1,1,2,1)),
				 "3g.3a.5u"=data.frame(family=c(1,1,1,1,1,1,1,1), person=c(1,2,3,4,5,6,7,8), father=c(0,0,1,1,0,4,4,4), 
                                mother=c(0,0,2,2,0,5,5,5), sex=c(1,2,1,1,2,1,1,1), affect=c(1,2,1,2,1,2,1,1))
				 )
}

##gene drop simulation given the number of childern in every generation
# family_strct.2g2c <- data.frame(family=c(1,1,1,1), person=c(1,2,3,4), father=c(0,0,1,1), 
#                                    mother=c(0,0,2,2), sex=c(1,2,1,1), affect=c(1,2,2,2)) #1=male, 2=female, 1=unaffected, 2=affected
#use the file format as in Merlin
#simulate the tranmission vector by generation and determine the affected status
#this generate the whole family first and then keep only those matched input
##consider polygenic effect when generating family data
##3. simulate data using the new model
gene_family_pe <- function(family_strct=family_strct_ped, n_family= 1000, trait_mean=0, Beta=Beta, dis_cutoff=NA) {
	library(kinship2) #load kinship function to calculate kinship matrix
	#get genotype matrix
	n_haplo <- nrow(haplotype)
	n_snp <- nrow(snp)
  n_family_member <- length(family_strct$person)
  data_family <- matrix(NA, nrow=n_family*n_family_member, ncol=(6+n_snp*2))
  tran_vec <- matrix(NA, nrow=n_family*n_family_member, ncol=3)
  affect_spec <- family_strct$affect-1 # afffect statuts so that 1 is affected 0 is unaffected
  haplotype.effect <- as.matrix(2- haplotype[, -(1:2)]) %*% Beta #precalculate to speed up
  #the strategy is generate the whole family and check affected status if fail then restart from the first individual
  data_family.idx <- 1 #index of the current family member being generated
  n_family.idx <- 1 #index of how many families have been generated
  
  ##parameters to generate phenotype
  pe_var <- 0.5 #p*(1-p)/(exp(2*p)/(exp(p)+1)^4) to achieve 50% heritability
	e_var <- 0.5 
  beta0 <- optimize(function(y) (trait_mean-integrate(function(x) x*dnorm(x, y, pe_var + e_var), lower = -Inf, upper = Inf)$value)^2, c(-10,10))$minimum #baseline prevalence
  K <- kinship(with(family_strct, pedigree(person, father, mother, sex)))
	PE <- 2*K*pe_var #polygenic effect
	E <- diag(nrow(K))*e_var
	PE_E <- PE + E

  while(n_family.idx <= n_family) {
    disease_vec <- matrix(NA, nrow=n_family_member, ncol=1) #store generated affected status
    family.haplo <- matrix(NA, nrow=n_family_member, ncol=2) #store haplotype info.
    ind.idx <- 1 # which individual is generating
    while(ind.idx <= n_family_member) { #until the generated matches the input
      #generate the current individual
      if(family_strct$father[ind.idx]==0 & family_strct$mother[ind.idx]==0) { #if founder
        family.haplo[ind.idx,] <- sample(1:n_haplo, 2, replace=T)
      } else{ #if not founder
        family.haplo[ind.idx,] <- c(sample(family.haplo[family_strct$father[ind.idx],], 1), sample(family.haplo[family_strct$mother[ind.idx],], 1))
      }
    	ind.idx <- ind.idx + 1
    }
		#save haplotype effect
    family.haplo.effect <- apply(family.haplo, 1, function(x) haplotype.effect[x[1]] + haplotype.effect[x[2]])
    
    #generate phenotype
    mu <- beta0 + family.haplo.effect
		trait <- MASS::mvrnorm(1, mu=mu, Sigma=PE_E)
		dis_status <- (trait > dis_cutoff) +0
		
		if(all(dis_status==affect_spec) | all(affect_spec==-1)) {
	    #save the haplotype file
	    letter.idx <- 1 #indicator used in the transmission vector
	    for(i in 1:n_family_member) {
	      #store genope and transmission vector
	      data_family[data_family.idx, ] <- unlist(c(n_family.idx, family_strct[i,2:5], trait[i], haplotype[family.haplo[i,1],-c(1:2)], haplotype[family.haplo[i,2],-c(1:2)]))
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
  }
  colnames(data_family) <- c("family","person","father","mother","sex","affect",rep(paste("SNP", 1:n_snp, sep=""),2))
  colnames(tran_vec) <- c("family","h1","h2")
  return(list(data_family=data.frame(data_family, stringsAsFactors=F), tran_vec=data.frame(tran_vec, stringsAsFactors=F)))
}
# Beta <- rep(0, length.out=50)
# Beta[7] <- log(1) #effect size of OR=2
# gene_family_pe(family_strct=family_strct_ped, n_family= 100, trait_mean=0.1, Beta=Beta)


#convernt family data from haploid to diploid format
hap2dip <- function(data=family_generated_2g3c, risk.variant.id, save.file=F, psuedo.snp.include=F) {
  n_snp <- (ncol(data$data_family)-6)/2
  n_test_snp <- length(risk.variant.id)
  n_sample <- nrow(data$data_family)
  geno <- matrix(NA, nrow=n_sample, ncol=n_snp)
  psuedo_snp <- rep(0, length.out = n_sample)
  
  ped <- data$data_family[, 1:6]
  names(ped) <- c("ped", "person", "father", "mother", "sex", "trait")
  #change to 1:affect 0:unaffected
  ped$trait <- ped$trait - 1
  
  geno <- 4 - data$data_family[, 6 + risk.variant.id] - data$data_family[, (6+n_snp+risk.variant.id)] 
  if(psuedo.snp.include==T) {
    geno <- cbind(ped[, c("ped","person")], geno, psuedo_snp)
  } else {
    geno <- cbind(ped[, c("ped","person")], geno)
  }             
  
  if(save.file==T) {
    ##output only snps of interest + psuedo_snp for cases there are only one variant of interest
    # geno_full <- 4 - data$data_family[, 7:(n_snp+6)] - data$data_family[, (6+n_snp+1):(6+2*n_snp)] 
    # ped_geno <- cbind(ped,geno_full)
    
    ##look up each offspring's father and mother and save the index
    ped_temp <- cbind(as.numeric(rownames(ped)), ped) #add row number
    idx <- apply(ped_temp, 1, function(x) if(x[3]!=0 & x[4]!=0) c(x[1], which(ped$ped==x[2] & ped$person==x[4]), which(ped$ped==x[2] & ped$person==x[5]))) #when non-founder, output the rownumber of offspring, father, and mother
    idx <- unname(unlist(idx))
        
    ped_geno <- cbind(ped,geno[, -c(1:2)])
    ped_geno <- ped_geno[idx, ]
    write.table(ped_geno, file = paste("data_",seed,".ped", sep=""), quote = F, sep = " ", row.names = F, col.names = F)
  }
  
  #return result
  list(ped=ped, geno=geno)
}
#hap2dip()

##regular family test trap
family.test <- function(data=family_generated, f=risk.variant.id, nofounderphenotype=F, conditional.test=F, prop.test=F, lm.test=F, noIBD=F, summary.stat="1") {
  ##hypothesis testing count the number of carrier haplotypes transmitted to the offsprings
  n_family <- max(data$data_family$family)
  n_family_member <- table(data$data_family$family)
  #check if founder's haplotype carries any variant's with f < 0.1
  if(length(f)==1 & f[1] <1) { #support allele frequency or list of snps
    snp2look.idx <-  which(snp$FREQ1 < f) # snp to look for
  } else(snp2look.idx <-  f)
  
  n_carrier_family <- vector("integer", n_family) #number of carrier founder haplotype in each family

  summary.stat <- switch (summary.stat,
  	"1" = function(tran_vec, carrier, affect) {sum(tran_vec[, c("h1","h2")] %in% carrier)*(affect-1)}, #only affected
  	"2" = function(tran_vec, carrier, affect) {sum((tran_vec[, c("h1","h2")] %in% carrier)*2-1)*(affect)}, #both affected and unaffected by halotype
  	"3" = function(tran_vec, carrier, affect) {(any(tran_vec[, c("h1","h2")] %in% carrier)*2-1)*(affect)} #both affected and unaffected by individual
  )
  
  
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
    n_carrier_family[family.idx] <<- n_carrier <- sum(!is.na(carrier)) #no. of carrier haplotypes
    
    if(noIBD==F) {
      if(conditional.test==F) {
        start.idx <- ifelse(nofounderphenotype==F, 1, 3)
        criteria <- !(n_carrier==(2*n_founder) | n_carrier==0)
      }else{
        n_carrier_offspring <- 0
        start.idx <- 3 #force conditonal test to ignore founder's phenotype
        for(i in start.idx:n_family_member[family.idx]) {
          #print(c(sum(tran_vec[current_row+i, c("h1","h2")] %in% carrier),(affect[current_row+i]-1)))
          n_carrier_offspring <- n_carrier_offspring + summary.stat(tran_vec[i, ], carrier, affect[i])
        }
        criteria <- n_carrier_offspring>0
      }
  
      if(criteria) { #skip families with 0 or 4 carrier haplotypes in founders or no carrier in children for conditonal test
        IBD_haplotype_observed = 0
        for(i in start.idx:n_family_member[family.idx]) {
          #print(c(sum(tran_vec[current_row+i, c("h1","h2")] %in% carrier),(affect[current_row+i]-1)))
            IBD_haplotype_observed <- IBD_haplotype_observed + summary.stat(tran_vec[i, ], carrier, affect[i]) 
        }
        observed <- IBD_haplotype_observed
        
        #calculate expectation and variance
        founder <- t(combn(LETTERS[1:(2*n_founder)], n_carrier))
        S <- apply(founder, 1, function(x) {
          carrier <- x #founder's haplotype
          #       print(carrier)
          IBD_haplotype_observed = 0
          for(i in start.idx:n_family_member[family.idx]) {
            #print(c(sum(tran_vec[current_row+i, c("h1","h2")] %in% carrier),(affect[current_row+i]-1)))
            IBD_haplotype_observed <- IBD_haplotype_observed + summary.stat(tran_vec[i, ], carrier, affect[i])
          }
          IBD_haplotype_observed
        }
        )
        if(conditional.test==T) {
          S <- S[which(S!=0)] #only count the configuration when there is at least one variant observed in the sibpair
        }
        
        mean_S <- mean(S)
        var_S <- sum((S-mean(S))^2)/nrow(founder)
        c(observed=observed, mean=mean_S, var=var_S, n_carrier=n_carrier, family.idx=family.idx)
      }
    }else{
      if(noIBD==T) { #assuming only genotype information is known          
        if(n_carrier == 1) { #only consider the family with one carrier haplotype in founder
          IBD_haplotype_observed = 0
          for(i in 3:n_family_member[family.idx]) { #check offspring's genotype
            #print(c(sum(tran_vec[current_row+i, c("h1","h2")] %in% carrier),(affect[current_row+i]-1)))
            IBD_haplotype_observed <- IBD_haplotype_observed + summary.stat(tran_vec[i, ], carrier, affect[i])
          }
          observed <- IBD_haplotype_observed
          
          mean_S <- 1
          if(observed==2){
            var_S = 0.75
          }
          if(observed==1){
            var_S = 0.25
          }
          if(observed==0){
            var_S = 0.75
          }
          c(observed=observed, mean=mean_S, var=var_S, n_carrier=n_carrier, family.idx=family.idx)
        }
      }
    }
  })  
  
  test.stat <- data.frame(do.call(rbind, test.stat))
  
  v <- test.stat$var
  se <- sqrt(sum(v))
  #e_avg <- mean(e) #overall expectation
  final.test.stat <- sum(test.stat$observed - test.stat$mean)/se
  
  #c(t, n_test.data)#T and number of informative families
  p.value <- 2*pnorm(abs(final.test.stat), lower=F) #p-value


  if(prop.test==T) { #fisher's method to combine two p-values
    #association test using founder treadted as all affected
  	founder.idx <- which(data$data_family$father==0 & data$data_family$mother==0)
  	h1 <- data$data_family[founder.idx,7:(6+n_snp)]
    h2 <- data$data_family[founder.idx,-c(1:(6+n_snp))]
   	carrier <- apply(cbind(h1[, snp2look.idx], h2[, snp2look.idx]), 1, function(x) sum(x==1)>0) 
    p.value.association <- prop.test(sum(carrier), length(carrier), p=1-dbinom(0, 2,risk.haplo.f))$p.value #assume the know pop f
    
    #fisher's method
    p <- c(p.value, p.value.association)
    Xsq <- -2*sum(log(p))
    p.val <- pchisq(Xsq, df = 2*length(p), lower.tail = FALSE)
    return(list(Xsq = Xsq, p.value = p.val))
  }
  
  if(lm.test==T) { #fisher's method to combine two p-values
    #association test using founder treadted as all affected 
		founder.idx <- which(data$data_family$father==0 & data$data_family$mother==0)
  	
		data.cc.diploid <- hap2dip(data=list(data_family=data$data_family[founder.idx, ]), risk.variant.id=risk.variant.id, save.file=F)
	  obj <- SKAT_Null_Model(data.cc.diploid$ped$trait ~ 1, out_type="C")
  	p.value.between <- SKAT(as.matrix(data.cc.diploid$geno[, -c(1:2)]), obj, weights.beta = c(1,1), r.corr = 1)$p.value

    #fisher's method
    p <- c(p.value, p.value.between)
    Xsq <- -2*sum(log(p))
    p.val <- pchisq(Xsq, df = 2*length(p), lower.tail = FALSE)
    return(list(Xsq = Xsq, p.value = p.val))
  }

  list(final.test.stat=final.test.stat, sum_e=sum(test.stat$mean), se=se, mean_observed=mean(test.stat$observed), mean_mean=mean(test.stat$mean), mean_var=mean(test.stat$var), p.value=p.value, n_info_family=length(test.stat$n_carrier))
}
# family.test()

##family test in TRAFIC spirit
family.test.trafic.ext <- function(data=family_generated, f=risk.variant.id, nofounderpheno=F, test_sq=F, prop.test=F) {
  ##hypothesis testing count the number of carrier haplotypes transmitted to the offsprings
  n_family <- max(data$data_family$family)
  n_family_member <- table(data$data_family$family)
  #check if founder's haplotype carries any variant's with f < 0.1
  if(length(f)==1 & f[1] <1) { #support allele frequency or list of snps
    snp2look.idx <-  which(snp$FREQ1 < f) # snp to look for
  } else(snp2look.idx <-  f)
  
  
  #start looking at each family
  data.family <- lapply(1:n_family, function(x) {
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
    
    #number of unique haplotype among affected family members
    if(nofounderpheno==F) {
      affect_id <- which(affect==2)
    }else {
      affect_id <- which(affect[-c(1:2)]==2)+2
    }
    unique_haplotype <-sort(unique(as.vector(as.matrix(tran_vec[affect_id,c("h1", "h2")]))))
    
    
    #calculate for each unique haplotype the number of affeted member and report if it is a carrier haplotype
    stat <- lapply(unique_haplotype, function(x) {
      idx <- which(tran_vec[,c("h1", "h2")]==x, arr.ind=T)[1,]
      haplotype <- switch(idx[2], "1"=h1[idx[1],], "2"=h2[idx[1],])
      carrier <- sum(haplotype[, snp2look.idx]==1)>0
      haplotype_on_affect <- sum(tran_vec[affect_id,c("h1", "h2")] == x)
      data.frame(family.idx=family.idx, x=x, haplotype_on_affect=haplotype_on_affect, carrier=carrier)
    })
    
    
    return(do.call(rbind,stat))
  })
  data.family <- do.call(rbind,data.family) #convert a list of dataframes to a dataframe  
  data.family <- transform(data.family, haplotype_on_affect_sq=haplotype_on_affect^2)
#   table(data.family$carrier, data.family$haplotype_on_affect)

  #fit a logistic regression
  if(test_sq==F) {
    glm.result <- summary(glm(carrier ~ haplotype_on_affect, family=binomial(link = "logit"), data=data.family))  
  }else{
    glm.result <- summary(glm(carrier ~ haplotype_on_affect_sq, family=binomial(link = "logit"), data=data.family))
  }
  p.value <- glm.result$coefficients[2, "Pr(>|z|)"]   

  if(prop.test==T) { #fisher's method to combine two p-values
    n_carrier_family <- vector("integer", n_family) #number of carrier founder haplotype in each family
    
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
      n_carrier_family[family.idx] <<- n_carrier <- sum(!is.na(carrier)) #no. of carrier haplotypes
    })
      
    #association test using founder treadted as all affected 
    p.value.association <- prop.test(sum(n_carrier_family), 4*n_family, p=risk.haplo.f)$p.value #assume the know pop f
    
    #fisher's method
    p <- c(p.value, p.value.association)
    Xsq <- -2*sum(log(p))
    p.val <- pchisq(Xsq, df = 2*length(p), lower.tail = FALSE)
    return(list(Xsq = Xsq, p.value = p.val))
  }

  
  list(p.value=glm.result$coefficients[2, "Pr(>|z|)"])
}
# family.test.trafic.ext()


##generate population sample
gene_case_control <- function(n_case_control_pair=1000, n_case=NA, n_control=NA) {
  if(is.na(n_case)) {
    n_case <- n_case_control_pair
    n_control <- n_case_control_pair
  }
  data_case_control <- matrix(NA, nrow=n_case+n_control, ncol=(2+n_snp*2))
  #generate cases first
  current_row <- 1 #which row is generating
  n_case.idx <- 1 #index of how many families have been generated
  while(n_case.idx <= n_case) {
    disease <<- 0
    while(disease!=2) {
      haplo.id <- sample(1:n_haplo, 2, replace=T)
      disease_prob <- prod(haplotype.risk[haplo.id])
      disease_prob <- ifelse(disease_prob <= 1, disease_prob, 1)
      disease <<- rbinom(1,1, prob=disease_prob) + 1
    }
    #save the haplotype file
    data_case_control[current_row, ] <- unlist(c(n_case.idx, disease, haplotype[haplo.id[1],-c(1:2)], haplotype[haplo.id[2],-c(1:2)]))
    current_row <- current_row + 1
    n_case.idx <- n_case.idx + 1
  }
  #generate controls first
  n_control.idx <- 1 #index of how many families have been generated
  while(n_control.idx <= n_control) {
    disease <<- 0
    while(disease!=1) {
      haplo.id <- sample(1:n_haplo, 2, replace=T)
      disease_prob <- prod(haplotype.risk[haplo.id])
      disease_prob <- ifelse(disease_prob <= 1, disease_prob, 1)
      disease <<- rbinom(1,1, prob=disease_prob) + 1
    }
    #save the haplotype file
    data_case_control[current_row, ] <- unlist(c(n_control.idx, disease, haplotype[haplo.id[1],-c(1:2)], haplotype[haplo.id[2],-c(1:2)]))
    current_row <- current_row + 1
    n_control.idx <- n_control.idx + 1
  }
  
  colnames(data_case_control) <- c("id","affect", rep(paste("SNP", 1:n_snp, sep=""),2))
  return(data_case_control=data.frame(data_case_control, stringsAsFactors=F))
}
# generated_case_control <- gene_case_control()

gene_case_control_pe <- function(n_case_control_pair=NA, n_case=NA, n_control=NA, trait_mean=0, Beta=Beta, dis_cutoff=NA) {
  #initialization
	n_haplo <- nrow(haplotype)
	n_snp <- nrow(snp)
	haplotype.effect <- as.matrix(2- haplotype[, -(1:2)]) %*% Beta #precalculate to speed up
	
	if(is.na(n_case)) { #assume equal case and conrtol
    n_case <- n_case_control_pair
    n_control <- n_case_control_pair
  }
	
	if(is.na(dis_cutoff) | dis_cutoff=="NA") { #do not specify cutoff
		dis_cutoff <- -9999999
		n_case <- n_case + n_control
		n_control <- 0
	}
	
	##parameters to generate phenotype
	pe_var <- 0.5 #p*(1-p)/(exp(2*p)/(exp(p)+1)^4) to achieve 50% heritability
	e_var <- 0.5 
	beta0 <- optimize(function(y) (trait_mean-integrate(function(x) x*dnorm(x, y, pe_var + e_var), lower = -Inf, upper = Inf)$value)^2, c(-10,10))$minimum #baseline prevalence
	PE <- pe_var #polygenic effect
	E <- e_var
	PE_E <- PE + E
	
  data_case_control <- matrix(NA, nrow=n_case+n_control, ncol=(6+n_snp*2))
  #generate cases first
  current_row <- 1 #which row is generating
  n_case.idx <- 1 #index of how many families have been generated
  while(n_case.idx <= n_case) {
    trait <<- -Inf
    while(trait < dis_cutoff) {
	    haplo.id <- sample(1:n_haplo, 2, replace=T)
	    sample.haplo.effect <- sum(haplotype.effect[haplo.id])
	    mu <- beta0 + sample.haplo.effect
	    trait <<- rnorm(1, mean=mu, sd=sqrt(PE_E))
    }
    #save the haplotype file
    data_case_control[current_row, ] <- unlist(c(current_row,1,0,0,1,trait, haplotype[haplo.id[1],-c(1:2)], haplotype[haplo.id[2],-c(1:2)]))
    current_row <- current_row + 1
    n_case.idx <- n_case.idx + 1
  }
  #generate controls 
  n_control.idx <- 1 #index of how many families have been generated
  while(n_control.idx <= n_control) {
    trait <<- Inf
    while(trait > dis_cutoff) {
	    #generate phenotype
      haplo.id <- sample(1:n_haplo, 2, replace=T)
      sample.haplo.effect <- sum(haplotype.effect[haplo.id])
      mu <- beta0 + sample.haplo.effect
      trait <<- rnorm(1, mean=mu, sd=sqrt(PE_E))
    }
    #save the haplotype file
    data_case_control[current_row, ] <- unlist(c(current_row,1,0,0,1, trait, haplotype[haplo.id[1],-c(1:2)], haplotype[haplo.id[2],-c(1:2)]))
    current_row <- current_row + 1
    n_control.idx <- n_control.idx + 1
  }
  
  colnames(data_case_control) <- c("family","person","father","mother","sex","affect", rep(paste("SNP", 1:n_snp, sep=""),2))
  return(data_case_control=data.frame(data_case_control, stringsAsFactors=F))
}
# generated_case_control_pe <- gene_case_control_pe()

case_control.test <- function(data=generated_case_control, f=0.01) {
  ##hypothesis testing count the number of carrier haplotypes transmitted to the offsprings
  n_case <- sum(data$affect==2)
  n_control <- nrow(data) - n_case
  #check if founder's haplotype carries any variant's with f < 0.1
  if(length(f)==1 & f[1] <1) { #support allele frequency or list of snps
    snp2look.idx <-  which(snp$FREQ1 < f) # snp to look for
  } else(snp2look.idx <-  f)
  
  #check if carrier haplotype
  carrier <- apply(data, 1, function(x) {
    h1 <- x[3:(2+n_snp)]
    h2 <- x[-(1:(2+n_snp))]
    #     sum(h1[snp2look.idx]==1)>0
    #     sum(h2[snp2look.idx]==1)>0
    c(h1.carrier=sum(h1[snp2look.idx]==1)>0, h2.carrier=sum(h2[snp2look.idx]==1)>0)
  })
  
  carrier.case <- sum(carrier[,1:n_case])
  carrier.control <- sum(carrier[,(n_case+1):(nrow(data))])
  test.result <- fisher.test(matrix(c(carrier.case, carrier.control, 2*n_case - carrier.case, 2*n_control - carrier.control), nrow=2)) #turn off correct to avoid a conservative test 
  c(test.result$estimate[1],test.result$estimate[2], p.value=test.result$p.value)
}
# case_control.test()

