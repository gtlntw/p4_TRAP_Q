#one affected founder and three affected 
#2g.3a.1u exact_phenotype=True, no-tie
rank_sum <- function(x,y){
	n_g1 <- length(x)
	n_g2 <- length(y)
	U_sum <- n_g1*n_g2
	U1 <- sum(sapply(x, function(z) sum(z>y)))
	U2 <- U_sum -U1
	U <- U1
	U_std <- U1/U_sum
	list(U=U,U_std=U_std)
}
bias_stat <- function(y=c(99,101.1, 101.3, 101.2), x=c(1,-1), center_shift=0){
	y <- y - center_shift
	a <- y[1]*x[1] + y[2]*x[2] + y[3]*x[2] + y[4]*x[1]
	b <- y[1]*x[1] + y[2]*x[2] + y[3]*x[1] + y[4]*x[2]
	c <- y[1]*x[2] + y[2]*x[1] + y[3]*x[1] + y[4]*x[1]
	d <- y[1]*x[2] + y[2]*x[1] + y[3]*x[2] + y[4]*x[2]
	obs <- c(a,b,c,d)
	mu <- mean(obs)
	var <- sum((obs-mu)^2)/4
	stat <- 2*(c-mu)/sqrt(2*var)
	
	a <- rank_sum(x=c(y[1],y[4]), y=c(y[2],y[3]))
	b <- rank_sum(x=c(y[1],y[3]), y=c(y[2],y[4]))
	c <- rank_sum(x=c(y[2],y[3],y[4]), y=c(y[1]))
	d <- rank_sum(x=c(y[2]), y=c(y[1],y[3],y[4]))
	obs_rank <- c(a$U,b$U,c$U,d$U)
	obs_rank_std <- c(a$U_std,b$U_std,c$U_std,d$U_std)
	list(obs=obs, mu=mu, var=var, stat=stat, obs_rank=obs_rank, obs_rank_std=obs_rank_std)
}

bias_stat(x=c(1,-1))
bias_stat(x=c(1,-1), center_shift = 100)



#one affected founder and three affected 
#2g.1a.3u exact_phenotype=True, no-tie
rank_sum <- function(x,y){
	n_g1 <- length(x)
	n_g2 <- length(y)
	U_sum <- n_g1*n_g2
	U1 <- sum(sapply(x, function(z) sum(z>y)))
	U2 <- U_sum -U1
	U <- U1
	U_std <- U1/U_sum
	list(U=U,U_std=U_std)
}
bias_stat <- function(y=c(101,99.1, 99.3, 99.2), x=c(1,-1), center_shift=0){
	y <- y - center_shift
	a <- y[1]*x[1] + y[2]*x[2] + y[3]*x[2] + y[4]*x[1]
	b <- y[1]*x[1] + y[2]*x[2] + y[3]*x[1] + y[4]*x[2]
	c <- y[1]*x[2] + y[2]*x[1] + y[3]*x[1] + y[4]*x[1]
	d <- y[1]*x[2] + y[2]*x[1] + y[3]*x[2] + y[4]*x[2]
	obs <- c(a,b,c,d)
	mu <- mean(obs)
	var <- sum((obs-mu)^2)/4
	stat <- 2*(c-mu)/sqrt(2*var)
	
	a <- rank_sum(x=c(y[1],y[4]), y=c(y[2],y[3]))
	b <- rank_sum(x=c(y[1],y[3]), y=c(y[2],y[4]))
	c <- rank_sum(x=c(y[2],y[3],y[4]), y=c(y[1]))
	d <- rank_sum(x=c(y[2]), y=c(y[1],y[3],y[4]))
	obs_rank <- c(a$U,b$U,c$U,d$U)
	obs_rank_std <- c(a$U_std,b$U_std,c$U_std,d$U_std)
	list(obs=obs, mu=mu, var=var, stat=stat, obs_rank=obs_rank, obs_rank_std=obs_rank_std)
}

bias_stat(x=c(1,-1))
bias_stat(y=c(101,99.1, 99.3, 99.2), x=c(1,-1), center_shift = 100)


