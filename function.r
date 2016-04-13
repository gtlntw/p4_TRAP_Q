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
gene_family_pe <- function(family_strct=family_strct_ped, n_family= 1000, trait_mean=0, Beta=Beta, dis_cutoff=NA, exact_affected=F) {
	library(kinship2) #load kinship function to calculate kinship matrix
	#get genotype matrix
	n_haplo <- nrow(haplotype)
	n_snp <- nrow(snp)
  n_family_member <- length(family_strct$person)
  data_family <- matrix(NA, nrow=n_family*n_family_member, ncol=(6+n_snp*2))
  tran_vec <- matrix(NA, nrow=n_family*n_family_member, ncol=3)
  affect_spec <- family_strct$affect-1 # afffect statuts so that 1 is affected 0 is unaffected -1 is unspecified
  n_affect_spec <- sum(affect_spec) # number of affected
  haplotype.effect <- as.matrix(2- haplotype[, -(1:2)]) %*% Beta #precalculate to speed up
  #the strategy is generate the whole family and check affected status if fail then restart from the first individual
  data_family.idx <- 1 #index of the current family member being generated
  n_family.idx <- 1 #index of how many families have been generated
  if(dis_cutoff == "NA") {
  	dis_cutoff <- NA
  	affect_spec <- rep(-1,n_family_member)
  	exact_affected <- F
  }
  
  #generate exact affected members or just the total number of affected
  if(exact_affected == T) {
    criteria <- expression(all(dis_status==affect_spec))
  } else {
    criteria <- expression(n_affect == n_affect_spec)
  }
  
  
  ##parameters to generate phenotype
  pe_var <- 0.5 #p*(1-p)/(exp(2*p)/(exp(p)+1)^4) to achieve 50% heritability
	e_var <- 0.5 
  beta0 <- trait_mean
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
		dis_status <- (trait > dis_cutoff + trait_mean) + 0
		(n_affect <- sum(dis_status))
		
		if(eval(criteria) | all(affect_spec==-1)) {
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
if(FALSE) {
  f=0.01
  seed=1000
  n_family=1000
  family_strct = "2g.2f.2c" #"2g.2a.2u" "2g.2f.2c" "3g.3a.4u"
  trait_mean=0
  r=0.3
  dis_cutoff = "NA" #"NA" qnorm(0.01, lower=F)
  exact_affected=F
}

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
family.test <- function(data=family_generated, f=risk.variant.id, nofounderphenotype=F, summary.stat="3", lm.test=F) {
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
  	"3" = function(tran_vec, carrier, affect) {(any(tran_vec[, c("h1","h2")] %in% carrier)*2-1)*(affect)}, #both affected and unaffected by individual
  	"4" = function(tran_vec, carrier, affect) {(any(tran_vec[, c("h1","h2")] %in% carrier)*2)*(affect)} #both affected and unaffected by individual with 2,0 coding for carrier
  	
  )
  
  if(lm.test==T) {nofounderphenotype <- T}
  
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
    
    start.idx <- ifelse(nofounderphenotype==F, 1, 3)
    criteria <- !(n_carrier==(2*n_founder) | n_carrier==0)

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
      mean_S <- mean(S)
      var_S <- sum((S-mean(S))^2)/nrow(founder)
      c(observed=observed, mean=mean_S, var=var_S, n_carrier=n_carrier, family.idx=family.idx)
    }
  })  
  
  test.stat <- data.frame(do.call(rbind, test.stat))
  carrier.family.idx <- test.stat$family.idx
  
  v <- test.stat$var
  se <- sqrt(sum(v))
  #e_avg <- mean(e) #overall expectation
  final.test.stat <- sum(test.stat$observed - test.stat$mean)/se
  
  #c(t, n_test.data)#T and number of informative families
  p.value <- 2*pnorm(abs(final.test.stat), lower=F) #p-value

  ##reactivated 04/07/2016!! this part was obsolete
  if(TRUE) {
    prop.test=F
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
      #lm test on all founders and child in the family with non-carrier founder
  		lm.idx <- which((data$data_family$family %in% carrier.family.idx & data$data_family$father==0 & data$data_family$mother==0) |
  												 	(!data$data_family$family %in% carrier.family.idx))
    	data.cc.diploid <- hap2dip(data=list(data_family=data$data_family[lm.idx, ]), risk.variant.id=risk.variant.id, save.file=F)
  	  result.pedgene <- pedgene(ped=data.cc.diploid$ped, geno=data.cc.diploid$geno)
    	p.value.between <- result.pedgene$pgdf$pval.burden
    	
    	#SKAT with all founders only
    	lm.idx <- which(data$data_family$father==0 & data$data_family$mother==0)
    	data.cc.diploid <- hap2dip(data=list(data_family=data$data_family[lm.idx, ]), risk.variant.id=risk.variant.id, save.file=F)
  	  obj <- SKAT_Null_Model(data.cc.diploid$ped$trait ~ 1, out_type="C")
      p.value.only.founder.SKAT <- SKAT(as.matrix(data.cc.diploid$geno[, -c(1:2)]), obj, weights.beta = c(1,1), r.corr = 1)$p.value
  
    	#fisher's method lm test on all founders and child in the family with non-carrier founder
      p <- c(p.value, p.value.between)
      Xsq <- -2*sum(log(p))
      p.val <- pchisq(Xsq, df = 2*length(p), lower.tail = FALSE)
      
      #fisher's method SKAT with all founders only
      p <- c(p.value, p.value.only.founder.SKAT)
      Xsq <- -2*sum(log(p))
      p.value.only.founder <- pchisq(Xsq, df = 2*length(p), lower.tail = FALSE)
      
      
      return(list(Xsq = Xsq, p.value = p.val, p.value.trap = p.value, p.value.lm = p.value.between,
      						p.value.only.founder = p.value.only.founder, p.value.only.founder.SKAT = p.value.only.founder.SKAT))
    }
  }

  list(final.test.stat=final.test.stat, sum_e=sum(test.stat$mean), se=se, mean_observed=mean(test.stat$observed), mean_mean=mean(test.stat$mean), mean_var=mean(test.stat$var), p.value=p.value, n_info_family=length(test.stat$n_carrier), info_family=test.stat$family.idx)
}
# family.test()

##rank_based family test trap
family.rank.test <- function(data=family_generated, f=risk.variant.id, nofounderphenotype=F) {
  ##hypothesis testing count the number of carrier haplotypes transmitted to the offsprings
  n_family <- max(data$data_family$family)
  n_family_member <- table(data$data_family$family)
  #check if founder's haplotype carries any variant's with f < 0.1
  if(length(f)==1 & f[1] <1) { #support allele frequency or list of snps
    snp2look.idx <-  which(snp$FREQ1 < f) # snp to look for
  } else(snp2look.idx <-  f)
  
  n_carrier_family <- vector("integer", n_family) #number of carrier founder haplotype in each family
  
  rank_sum <- function(x,y){
		n_g1 <- length(x)
		n_g2 <- length(y)
		U_sum <- n_g1*n_g2
		U1 <- sum(sapply(x, function(z) sum(z>y, 0.5*(z==y))))
		U2 <- U_sum -U1
		U <- U1
		U_std <- U1/U_sum
		list(U=U,U_std=U_std)
  }
  
  summary.stat <- function(tran_vec, carrier, affect) {
  	IBD_haplotype_observed = list()
  	for(i in start.idx:n_family_member[family.idx]) {
  		#print(c(sum(tran_vec[current_row+i, c("h1","h2")] %in% carrier),(affect[current_row+i]-1)))
  		IBD_haplotype_observed[[i]] <- data.frame(carrier=(tran_vec[i, c("h1","h2")] %in% carrier), affect=affect[i]) 
  	}
  	observed <- do.call(rbind, IBD_haplotype_observed)
  	
  	carrier.idx <- which(observed$carrier==T)
  	carrier_trait <- observed$affect[carrier.idx]
  	non_carrier.idx <- which(observed$carrier==F)
  	non_carrier_trait <- observed$affect[non_carrier.idx]
  	
  	if(any(length(carrier.idx)==0, length(non_carrier.idx)==0)) return(NA)
  	rank_sum(carrier_trait, non_carrier_trait)$U_std
  }

  #start looking at each family
  T_list <- list()
  test_stat <- list()
  for(x in 1:n_family) {
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
    n_carrier_family[family.idx] <- n_carrier <- sum(!is.na(carrier)) #no. of carrier haplotypes
    
    start.idx <- ifelse(nofounderphenotype==F, 1, 3)
    criteria <- !(n_carrier==(2*n_founder) | n_carrier==0)

    if(criteria) { #skip families with 0 or 4 carrier haplotypes in founders or no carrier in children for conditonal test
      observed <- summary.stat(tran_vec, carrier, affect)
      
      #calculate expectation and variance
      founder <- t(combn(LETTERS[1:(2*n_founder)], n_carrier))
      S <- apply(founder, 1, function(x) {
				summary.stat(tran_vec, x, affect)
      })
      rank_S <- rank(S)
      rank_observed <- rank_S[match(observed, S)]
      T_list[[family.idx]] <- S 
      test_stat[[family.idx]] <- data.frame(observed=observed, rank_observed=rank_observed, n_carrier=n_carrier, family.idx=family.idx, min_S=min(S), max_S=max(S), min_rank=min(rank_S), max_rank=max(rank_S))
    }
  } 
  
  T.comb <- T_list[!sapply(T_list, is.null)]
  T.comb <- lapply(T.comb, function(x) x[!is.na(x)])
  test.stat.table <- data.frame(do.call(rbind, test_stat))
  n_info_family <- length(test.stat.table$family.idx)
  
  final.test.stat <- sum(test.stat.table$observed)
  sum_min_S <- sum(test.stat.table$min_S)
  sum_max_S <- sum(test.stat.table$max_S)
  null.dist <- replicate(100000, {
  	sum(sapply(T.comb, function(x) sample(x,1)))
  })
  p.value <- mean(min(final.test.stat, sum_max_S-final.test.stat+sum_min_S) >= null.dist) +
  	mean(max(final.test.stat, sum_max_S-final.test.stat+sum_min_S) <= null.dist)
  
  final.test.stat.rank <- sum(test.stat.table$rank_observed)
  # max_rank_table <- table(test.stat.table$max_rank)
  sum_min_rank <- sum(test.stat.table$min_rank)
  sum_max_rank <- sum(test.stat.table$max_rank)
#   nrow_comb <- nrow(max_rank_table)
#   max_rank_count <- names(max_rank_table)
#   null.dist.rank <- x <- replicate(100000, {
#   	temp <- list()
#   	for(i in 1:nrow_comb) {
#   		temp[[i]] <- sum(sample(1:max_rank_count[i], max_rank_table[i], replace = T))
#   	}
#   	unlist(temp)
#   })
  
  T.comb.rank <- lapply(T.comb, function(x) sort(rank(x)))
  null.dist.rank <- replicate(100000, {
  	sum(sapply(T.comb.rank, function(x) sample(x,1)))
  })
  p.value.rank <- mean(min(final.test.stat.rank, sum_max_rank-final.test.stat.rank+sum_min_rank) >= null.dist.rank) +
  	mean(max(final.test.stat.rank, sum_max_rank-final.test.stat.rank+sum_min_rank) <= null.dist.rank)
  
  list(final.test.stat=final.test.stat, p.value=p.value, final.test.stat.rank=final.test.stat.rank, p.value.rank=p.value.rank, n_info_family=n_info_family, info_family=test.stat.table$family.idx)
}
# family.rank.test()


##rank_based family test trap using all samples
family.rank.test.all <- function(data=family_generated, f=risk.variant.id, nofounderphenotype=F) {
  ##hypothesis testing count the number of carrier haplotypes transmitted to the offsprings
  n_family <- max(data$data_family$family)
  n_family_member <- table(data$data_family$family)
  #check if founder's haplotype carries any variant's with f < 0.1
  if(length(f)==1 & f[1] <1) { #support allele frequency or list of snps
    snp2look.idx <-  which(snp$FREQ1 < f) # snp to look for
  } else(snp2look.idx <-  f)
  
  #calculate the mean rank sum of carrier
  rank_sum <- function(stats){
    carrier.idx <- which(stats$carrier==T)
  	carrier_trait <- stats$affect[carrier.idx]
  	non_carrier.idx <- which(stats$carrier==F)
  	non_carrier_trait <- stats$affect[non_carrier.idx]
  	
		n_g1 <- length(carrier_trait)
		n_g2 <- length(non_carrier_trait)
		U_sum <- n_g1*n_g2
		U1 <- sum(sapply(carrier_trait, function(z) sum(z>non_carrier_trait, 0.5*(z==non_carrier_trait))))
		U2 <- U_sum -U1
		U <- U1
		U_std <- U1/U_sum
		list(U=U,U_std=U_std)
  }
  
  #for within-family test
  summary.stat <- function(tran_vec, carrier, affect) {
  	IBD_haplotype_observed = list()
  	for(i in start.idx:n_family_member[family.idx]) {
  		#print(c(sum(tran_vec[current_row+i, c("h1","h2")] %in% carrier),(affect[current_row+i]-1)))
  		IBD_haplotype_observed[[i]] <- data.frame(carrier=(tran_vec[i, c("h1","h2")] %in% carrier), affect=affect[i]) 
  	}
  	observed <- do.call(rbind, IBD_haplotype_observed)
  	
  	carrier.idx <- which(observed$carrier==T)
  	non_carrier.idx <- which(observed$carrier==F)
  	if(any(length(carrier.idx)==0, length(non_carrier.idx)==0)) return(NA)
  	rank_sum(observed)$U_std
  }
  
  #calculate between-family statistic given carrier_list
  stat_between <- function(carrier_list) {
    start.idx <- ifelse(nofounderphenotype==F, 1, 3)
    # criteria <- !(n_carrier==(2*n_founder) | n_carrier==0)

    stat_list <- list()
    for(x in 1:n_family) {
      if(length(carrier_list[[x]])!=0) {
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

        n_founder <- length(carrier_list[[family.idx]])
        carrier <- LETTERS[carrier_list[[family.idx]]]
 
        #datafame of carrier status and trait
        IBD_haplotype_observed = list()
      	for(i in start.idx:n_family_member[family.idx]) {
      		#print(c(sum(tran_vec[current_row+i, c("h1","h2")] %in% carrier),(affect[current_row+i]-1)))
      		IBD_haplotype_observed[[i]] <- data.frame(carrier=(tran_vec[i, c("h1","h2")] %in% carrier), affect=affect[i]) 
      	}
      	observed <- do.call(rbind, IBD_haplotype_observed)
      	
        stat_list[[family.idx]] <- observed
      }
    }
    stats <- do.call(rbind, stat_list)
    rank_sum(stats)$U_std
  }
  
  #randomize carrier list
  random_carrier <- function() {
    n_carrier <- length(unlist(carrier_list))
    carrier_status <- rep(0, length=2*sum((n_founder_list))) #vector of founder chromosome
    carrier_idx <- sample(1:length(carrier_status), n_carrier) #sample which founder chromosome is carrier
    carrier_status[carrier_idx] <- 1 #update the carrier_status vector
    new_carrier_list <- list() #new carrier list
    founder_idx <- c(0,cumsum(2*n_founder_list)) #indexing for processing
    for(x in 1:n_family) {
      idx <- ((founder_idx[x]+1):founder_idx[x+1]) #process each family
      new_carrier_list[[x]] <- which(carrier_status[idx]==1) #return which founder chromosome is carrier in each family
    }
    return(new_carrier_list)
  }
  
  
  #create carrier list and n_carrier list
  #within
  T_list <- list()
  stat_list_within <- list()
  #between
  carrier_list <- list()
  stat_list_between <- list()
  n_founder_list <- vector("integer", n_family)
  for(x in 1:n_family) {
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
        c(ifelse(sum(h1[y, snp2look.idx]==1)>0, 2*y-1, NA),
          ifelse(sum(h2[y, snp2look.idx]==1)>0, 2*y, NA))
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
    
    carrier_list[[family.idx]] <- carrier[!is.na(carrier)]
    n_founder_list[family.idx] <- n_founder #no. of carrier haplotypes
    start.idx <- ifelse(nofounderphenotype==F, 1, 3)
    carrier <- LETTERS[carrier_list[[family.idx]]] #convert number to letter 

    ##betwee-family
    #datafame of carrier status and trait
    IBD_haplotype_observed = list()
  	for(i in start.idx:n_family_member[family.idx]) {
  		#print(c(sum(tran_vec[current_row+i, c("h1","h2")] %in% carrier),(affect[current_row+i]-1)))
  		IBD_haplotype_observed[[i]] <- data.frame(carrier=(tran_vec[i, c("h1","h2")] %in% carrier), affect=affect[i]) 
  	}
  	stat_list_between[[family.idx]] <- do.call(rbind, IBD_haplotype_observed)
    
    ##within-family
    n_carrier <- length(carrier) #no. of carrier haplotypes
    criteria <- !(n_carrier==(2*n_founder) | n_carrier==0)

    if(criteria) { #skip families with 0 or 4 carrier haplotypes in founders or no carrier in children for conditonal test
      observed <- summary.stat(tran_vec, carrier, affect)
      
      #calculate expectation and variance
      founder <- t(combn(LETTERS[1:(2*n_founder)], n_carrier))
      S <- apply(founder, 1, function(x) {
				summary.stat(tran_vec, x, affect)
      })
      rank_S <- rank(S)
      rank_observed <- rank_S[match(observed, S)]
      T_list[[family.idx]] <- S 
      stat_list_within[[family.idx]] <- data.frame(observed=observed, rank_observed=rank_observed, n_carrier=n_carrier, family.idx=family.idx, min_S=min(S), max_S=max(S), min_rank=min(rank_S), max_rank=max(rank_S))
    }
  } 
  
  ##between
  stat_temp <- do.call(rbind, stat_list_between)
  test_stat_between <- rank_sum(stat_temp)$U_std
  
  #permutation test
  # null.dist_between <- replicate(100, {
  # 	new_carrier_list <- random_carrier()
  #   stat_between(new_carrier_list)
  # })
  # p.value.between <- mean(min(test.stat, 1-test.stat) >= null.dist_between) +
  # 	mean(max(test.stat, 1-test.stat) <= null.dist_between)
  # 
  # list(test_stat_between=test_stat_between, p.value.between=p.value.between)
  # 
  ##within
  T.comb <- T_list[!sapply(T_list, is.null)]
  T.comb <- lapply(T.comb, function(x) x[!is.na(x)])
  test.stat.table <- data.frame(do.call(rbind, stat_list_within))
  n_info_family <- length(test.stat.table$family.idx)
  
  test_stat_within <- sum(test.stat.table$observed)
  sum_min_S <- sum(test.stat.table$min_S)
  sum_max_S <- sum(test.stat.table$max_S)
  # null.dist_within <- replicate(100000, {
  # 	sum(sapply(T.comb, function(x) sample(x,1)))
  # })
  # p.value_within <- mean(min(test_stat_within, sum_max_S-test_stat_within+sum_min_S) >= null.dist_within) +
  # 	mean(max(test_stat_within, sum_max_S-test_stat_within+sum_min_S) <= null.dist_within)

  test_stat_within.rank <- sum(test.stat.table$rank_observed)
  sum_min_rank <- sum(test.stat.table$min_rank)
  sum_max_rank <- sum(test.stat.table$max_rank)
  T.comb.rank <- lapply(T.comb, function(x) sort(rank(x)))
  # null.dist_within.rank <- replicate(100000, {
  # 	sum(sapply(T.comb.rank, function(x) sample(x,1)))
  # })
  # p.value_within.rank <- mean(min(test_stat_within.rank, sum_max_rank-test_stat_within.rank+sum_min_rank) >= null.dist_within.rank) +
  # 	mean(max(test_stat_within.rank, sum_max_rank-test_stat_within.rank+sum_min_rank) <= null.dist_within.rank)
  # 
  # list(test_stat_within=test_stat_within, p.value=p.value, test_stat_within.rank=test_stat_within.rank, p.value.rank=p.value.rank, n_info_family=n_info_family, info_family=test.stat.table$family.idx)
  
  list(test_stat_between, test_stat_within/n_info_family, test_stat_within.rank/n_info_family)
}
# family.rank.test.all()


#exact p-value calculation
exact.permutation <- function(t.obs=0, t.current=0, level=0, assignment=0, prune=T) {
  # print(paste("level=",level))
  #add the function call
  count.call <<- count.call + 1
  #initialize the calculation
 	if(level==0 & assignment==0) {
		t.current <- 0
	} else{
		t.current <- t.current + t.info.family[[level]][assignment]
	}
  #check if at the tip and as extreme or more extreme than t.obs
  if(level==n_info_family) {
    if(t.current > t.obs) {
		  return(1/n.permutation)
    } else if(t.current == t.obs) {
        return(0.5/n.permutation)
    } else if(t.current < t.obs) {
	      return(0)
      }
  } 
	t.max <- t.current + t.info.family.cummax[level+1] #max in this branch sum the max in the remaining levels
	t.min <- t.current + t.info.family.cummin[level+1]
	# print(paste(level, assignment, t.current, t.max, t.min))
	#prune the branch to speed up
	if(prune==T) {
	  #if min is greater than t.obs all the descendant counted as hits
  	if(t.min > t.obs) {
  		return(4^(n_info_family-level)/n.permutation)
  	}
  	#if max is smaller than t.obs, discard this branch
  	if(t.max < t.obs) {
  		return(0)
  	} 
	}
	#descend
	p.list <- list()
	for(i in 1:length(t.info.family[[level+1]])) {
	  p.list[[i]] <- exact.permutation(t.obs=t.obs, t.current=t.current, level=level+1, assignment=i, prune=prune)
	}
	p1 <- exact.permutation(t.obs=t.obs, t.current=t.current, level=level+1, assignment=1, prune=prune)
	p2 <- exact.permutation(t.obs=t.obs, t.current=t.current, level=level+1, assignment=2, prune=prune)
	p3 <- exact.permutation(t.obs=t.obs, t.current=t.current, level=level+1, assignment=3, prune=prune)
	p4 <- exact.permutation(t.obs=t.obs, t.current=t.current, level=level+1, assignment=4, prune=prune)
	p.value <- (p1+p2+p3+p4)
	return(p.value)
}

exact.permutation.2tails <- function(t.obs=0, t.current=0, level=0, assignment=0, prune=T) {
  # print(paste("level=",level))
  #add the function call
  count.call <<- count.call + 1
  #initialize the calculation
 	if(level==0 & assignment==0) {
		t.current <- 0
	} else{
		t.current <- t.current + t.info.family[[level]][assignment]
	}
  #check if at the tip and as extreme or more extreme than t.obs
  if(level==n_info_family) {
    if(t.current > right.cutoff | t.current < left.cutoff) {
		  return(1/n.permutation)
    } else if(t.current == right.cutoff | t.current == left.cutoff) {
        return(0.5/n.permutation)
    } else if(t.current < t.obs) {
	      return(0)
      }
  } 
	t.max <- t.current + t.info.family.cummax[level+1] #max in this branch sum the max in the remaining levels
	t.min <- t.current + t.info.family.cummin[level+1]
	# print(paste(level, assignment, t.current, t.max, t.min, t.obs))
	#prune the branch to speed up
	if(prune==T) {
	  #if min is greater than t.obs all the descendant counted as hits
  	if(t.min >= right.cutoff | t.max <= left.cutoff) {
  		return(4^(n_info_family-level)/n.permutation)
  	}
  	#if max is smaller than t.obs, discard this branch
  	if((t.max < right.cutoff) & (t.min > left.cutoff)) {
  		return(0)
  	} 
	}
	#descend
	p1 <- exact.permutation.2tails(t.obs=t.obs, t.current=t.current, level=level+1, assignment=1, prune=prune)
	p2 <- exact.permutation.2tails(t.obs=t.obs, t.current=t.current, level=level+1, assignment=2, prune=prune)
	p3 <- exact.permutation.2tails(t.obs=t.obs, t.current=t.current, level=level+1, assignment=3, prune=prune)
	p4 <- exact.permutation.2tails(t.obs=t.obs, t.current=t.current, level=level+1, assignment=4, prune=prune)
	p.value <- (p1+p2+p3+p4)
	return(p.value)
}

exact.pvalue <- function(t.stat, prune=T) {
  n.permutation <<- 4^length(t.info.family)
  n_info_family <<- length(t.info.family)
  
  count.call <<- 0
  result <- exact.permutation(t.obs=t.stat, prune=prune)
  number.function.calls <- count.call
  
  list(p.value=result, number.function.calls=number.function.calls)
}

compare_result <- function(n_info_family, t.stat) {
  n_info_family <- n_info_family
  t.stat=3.5*n_info_family
  t.info.family <<- rep(list(1:4), n_info_family) #place to save t at each level
  t.mean <- sum(sapply(t.info.family, mean))
  t.info.family.cummax <- rev(cumsum(rev(sapply(t.info.family, max))))
  t.info.family.cummin <- rev(cumsum(rev(sapply(t.info.family, min))))
  n.permutation <<- 4^length(t.info.family)
  n_info_family <<- length(t.info.family)
  left.cutoff <- min(t.stat, 2*t.mean - t.stat)
  right.cutoff <- max(t.stat, 2*t.mean - t.stat)
  exact.pvalue(t.stat=t.stat, prune=T)
  exact.permutation.2tails(t.obs=t.stat, prune=T)
  
  result <- list()
  time <- list() 
  time[["result"]] <- system.time(result[["result"]] <- unlist(exact.pvalue(t.stat=t.stat, prune=F)))
  time[["result.prune"]] <- system.time(result[["result.prune"]] <- unlist(exact.pvalue(t.stat=t.stat, prune=T)))
  
  time[["result.c"]] <- system.time(result[["result.c"]] <- unlist(exact_pvalue_c(t=t.info.family,t_stat=t.stat,prune=F)))
  time[["result.prune.c"]] <- system.time(result[["result.prune.c"]] <- unlist(exact_pvalue_c(t=t.info.family,t_stat=t.stat,prune=T)))
  list(result=result, time=time)
}
# compare_result(n_info_family=5, t.stat=2.5*n_info_family)

exact.permutation.collapse <- function(t.obs=0, t.current=0, level=0, assignment=0, n_assigned=0, prune=T) {
  #add the function call
  count.call <<- count.call + 1
  #initialize the calculation
 	if(level==0 & assignment==0) {
		t.current <- 0
	} else{
		t.current <- t.current + level*assignment
		n_assigned <- n_assigned + assignment
	}
  # if(level==4) print(paste("level=",level, " t.current=", t.current))
  #check if at the tip and as extreme or more extreme than t.obs
  if(level==4 & t.current >= t.obs) {
		return(1/n.permutation)
  } else if(level==4 & t.current < t.obs) {
	  return(0)
	}
	t.max <- t.current + (n_info_family-n_assigned)*4 #max in this branch sum the max in the remaining levels
	t.min <- t.current + (n_info_family-n_assigned)*(level+1)
	# print(paste(level, assignment, n_assigned, t.current, t.max, t.min))
	#prune the branch to speed up
	if(prune==T) {
	  #if min is greater than t.obs all the descendant counted as hits
  	# if(t.min >= t.obs) {
  	# 	return(4^(n_info_family-level)/n.permutation) #multiply by combination factor
  	# }
  	#if max is smaller than t.obs, discard this branch
  	if(t.max < t.obs) {
  		return(0)
  	} 
	}
	#descend
	p.value <- 0
	if(level==3) { #only consider the possible assignment in the 4th level if at level 3
	  p.value <- p.value + exact.permutation.collapse(t.obs=t.obs, t.current=t.current, level=level+1, assignment=(n_info_family-n_assigned), n_assigned=n_assigned, prune=prune) 
	} else{
  	for(i in 0:(n_info_family-n_assigned)) {
  	  p.value <- p.value + exact.permutation.collapse(t.obs=t.obs, t.current=t.current, level=level+1, assignment=i, n_assigned=n_assigned, prune=prune)*choose(n_info_family-n_assigned, i)  #multiply by combination factor
  	}
	}  
	return(p.value)
}
# count.call <<-0
# n_info_family <- 50
# t.stat=2.5*n_info_family
# t.info.family <<- rep(list(1:4), n_info_family) #place to save t at each level
# n.permutation <<- 4^length(t.info.family)
# exact.permutation.collapse(t.obs=t.stat)


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
	beta0 <- trait_mean
	PE <- pe_var #polygenic effect
	E <- e_var
	PE_E <- PE + E
	
  data_case_control <- matrix(NA, nrow=n_case+n_control, ncol=(6+n_snp*2))
  #generate cases first
  current_row <- 1 #which row is generating
  n_case.idx <- 1 #index of how many families have been generated
  while(n_case.idx <= n_case) {
    trait <<- -Inf
    while(trait < dis_cutoff + trait_mean) {
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
    while(trait > dis_cutoff  + trait_mean) {
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
  test.result <- fisher.test(matrix(c(carrier.case, carrier.control, n_case - carrier.case, n_control - carrier.control), nrow=2)) #turn off correct to avoid a conservative test 
  c(test.result$estimate[1],test.result$estimate[2], p.value=test.result$p.value)
}
# case_control.test()

