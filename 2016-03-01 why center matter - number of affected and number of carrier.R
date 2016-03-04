#one affected founder and one affected child

bias_stat <- function(y=c(1,-1), x=c(1,-1)){
	a <- y[2]*x[1] + y[1]*x[2] + y[1]*x[1] + y[2]*x[2]
	b <- y[2]*x[1] + y[1]*x[2] + y[1]*x[2] + y[2]*x[1]
	c <- y[2]*x[2] + y[1]*x[1] + y[1]*x[1] + y[2]*x[1]
	d <- y[2]*x[2] + y[1]*x[1] + y[1]*x[2] + y[2]*x[2]
	obs <- c(a,b,c,d)
	mu <- mean(obs)
	var <- sum((obs-mu)^2)/4
	stat <- 2*(c-mu)/sqrt(2*var)
	list(obs=obs, mu=mu, var=var, stat=stat)
}

bias_stat()


##variance plot as a function a
fn <- function(a, b) {
8*a^2 +8*a*b + 8*b^2
}
fn <- Vectorize(fn)

curve(sapply(x, function(x) fn(x, x-2)), -10,10, xlab="a", ylab="variance", main="a-b=2")

#find minimum
a <- seq(-10,10,by=1)
b <- a -2
fn(a,b)

##signal/noise ratio plot
t_fn <- function(a, b) {
  (2*a)/(8*a^2 +8*a*b + 8*b^2)
}
t_fn <- Vectorize(t_fn)

curve(sapply(x, function(x) t_fn(x, x-2)), -10,10, xlab="a", ylab="signal/noise ratio", main="a-b=2")
abline(v=1)
x_tick <- seq(-10, 10, by=1)
axis(1, at=x_tick)


##for families with two founder carrier haplotype
bias_stat <- function(y=c(1,-1), x=c(1,-1)){
	a <- y[2]*x[1] + y[1]*x[2] + y[1]*x[1] + y[1]*x[1]
	b <- y[2]*x[1] + y[1]*x[1] + y[1]*x[1] + y[1]*x[1]
	c <- y[2]*x[1] + y[1]*x[1] + y[1]*x[2] + y[1]*x[1]
	d <- y[2]*x[1] + y[1]*x[1] + y[1]*x[1] + y[1]*x[1]
	e <- y[2]*x[1] + y[1]*x[1] + y[1]*x[1] + y[1]*x[2]
	f <- y[2]*x[2] + y[1]*x[1] + y[1]*x[1] + y[1]*x[1]
	obs <- c(a,b,c,d,e,f)
	mu <- mean(obs)
	var <- sum((obs-mu)^2)/6
	stat <- 2*(c-mu)/sqrt(2*var)
	list(obs=obs, mu=mu, var=var, stat=stat)
}

bias_stat()


##variance plot as a function a
fn <- function(a,b, x=c(1,-1)) {
	y <- c(a,b)
	a <- y[2]*x[1] + y[1]*x[2] + y[1]*x[1] + y[1]*x[1]
	b <- y[2]*x[1] + y[1]*x[1] + y[1]*x[1] + y[1]*x[1]
	c <- y[2]*x[1] + y[1]*x[1] + y[1]*x[2] + y[1]*x[1]
	d <- y[2]*x[1] + y[1]*x[1] + y[1]*x[1] + y[1]*x[1]
	e <- y[2]*x[1] + y[1]*x[1] + y[1]*x[1] + y[1]*x[2]
	f <- y[2]*x[2] + y[1]*x[1] + y[1]*x[1] + y[1]*x[1]
	obs <- c(a,b,c,d,e,f)
	mu <- mean(obs)
	var <- sum((obs-mu)^2)/6
	var
}
fn <- Vectorize(fn)
curve(sapply(x, function(x) fn(x, x-2)), -10,10, xlab="a", ylab="variance", main="a-b=2")


#find minimum
a <- seq(-10,10,by=1)
b <- a -2
a[which(fn(a,b)==min(fn(a,b)))]

##signal/noise ratio plot
t_fn <- function(a,b, x=c(1,-1)) {
	y <- c(a,b)
	a <- y[2]*x[1] + y[1]*x[2] + y[1]*x[1] + y[1]*x[1]
	b <- y[2]*x[1] + y[1]*x[1] + y[1]*x[1] + y[1]*x[1]
	c <- y[2]*x[1] + y[1]*x[1] + y[1]*x[2] + y[1]*x[1]
	d <- y[2]*x[1] + y[1]*x[1] + y[1]*x[1] + y[1]*x[1]
	e <- y[2]*x[1] + y[1]*x[1] + y[1]*x[1] + y[1]*x[2]
	f <- y[2]*x[2] + y[1]*x[1] + y[1]*x[1] + y[1]*x[1]
	obs <- c(a,b,c,d,e,f)
	mu <- mean(obs)
	var <- sum((obs-mu)^2)/6
	b/var
}
t_fn <- Vectorize(t_fn)

curve(sapply(x, function(x) t_fn(x, x-2)), -10,10, xlab="a", ylab="signal/noise ratio", main="a-b=2")
abline(v=1)
x_tick <- seq(-10, 10, by=1)
axis(1, at=x_tick)

