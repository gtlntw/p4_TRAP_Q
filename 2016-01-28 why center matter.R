bias_stat <- function(y=c(1,-1), x=c(1,-1)){
	a <- y[2]*x[1] + y[1]*x[2] + y[1]*x[1] + y[1]*x[2]
	b <- y[2]*x[1] + y[1]*x[2] + y[1]*x[2] + y[1]*x[1]
	c <- y[2]*x[2] + y[1]*x[1] + y[1]*x[1] + y[1]*x[1]
	d <- y[2]*x[2] + y[1]*x[1] + y[1]*x[2] + y[1]*x[2]
	obs <- c(a,b,c,d)
	mu <- mean(obs)
	var <- sum((obs-mu)^2)/4
	stat <- 2*(c-mu)/sqrt(2*var)
	list(obs=obs, mu=mu, var=var, stat=stat)
}

bias_stat()


##variance plot as a function a
fn <- function(a, b) {
3*a^2 -2*a*b + b^2
}
fn <- Vectorize(fn)

curve(sapply(x, function(x) fn(x, x-2)), -10,10, xlab="a", ylab="variance", main="a-b=2")

#find minimum
a <- seq(-10,10,by=1)
b <- a -2
fn(a,b)

##signal/noise ratio plot
t_fn <- function(a, b) {
  (3*a-b)/(3*a^2 -2*a*b + b^2)
}
t_fn <- Vectorize(t_fn)

curve(sapply(x, function(x) t_fn(x, x-2)), -10,10, xlab="a", ylab="signal/noise ratio", main="a-b=2")
abline(v=1)
x_tick <- seq(-10, 10, by=1)
axis(1, at=x_tick)



#3dplot
library(plotly)
x = seq(-10,10, length.out = 100)
y = seq(-10,10, length.out = 100)
z = outer(x,y,fn)
colnames(z) <- y
rownames(z) <- x
plot_ly(z=z, type="surface", )



##for sibpair with IBD=2
bias_stat <- function(y=c(1,-1), x=c(1,-1)){
	a <- y[2]*x[1] + y[1]*x[2] + y[1]*x[1] + y[1]*x[1]
	b <- y[2]*x[1] + y[1]*x[2] + y[1]*x[2] + y[1]*x[2]
	c <- y[2]*x[2] + y[1]*x[1] + y[1]*x[1] + y[1]*x[1]
	d <- y[2]*x[2] + y[1]*x[1] + y[1]*x[2] + y[1]*x[2]
	obs <- c(a,b,c,d)
	mu <- mean(obs)
	var <- sum((obs-mu)^2)/4
	stat <- 2*(c-mu)/sqrt(2*var)
	list(obs=obs, mu=mu, var=var, stat=stat)
}

bias_stat()


##variance plot as a function a
fn <- function(a, b) {
4*b^2 -8*a*b + 20*a^2
}
fn <- Vectorize(fn)

curve(sapply(x, function(x) fn(x, x-2)), -10,10, xlab="a", ylab="variance", main="a-b=2")

#find minimum
a <- seq(-10,10,by=1)
b <- a -2
[which(fn(a,b)==min(fn(a,b)))]

##signal/noise ratio plot
t_fn <- function(a, b) {
  (-b + 3*a)/(4*b^2 -8*a*b + 20*a^2)
}
t_fn <- Vectorize(t_fn)

curve(sapply(x, function(x) t_fn(x, x-2)), -10,10, xlab="a", ylab="signal/noise ratio", main="a-b=2")
abline(v=1)
x_tick <- seq(-10, 10, by=1)
axis(1, at=x_tick)

