library(rethinking)
sppnames <- c( "afarensis","africanus","habilis","boisei",
               "rudolfensis","ergaster","sapiens")
brainvolcc <- c( 438 , 452 , 612, 521, 752, 871, 1350 )
masskg <- c( 37.0 , 35.5 , 34.5 , 41.5 , 55.5 , 61.0 , 53.5 )
d <- data.frame( species=sppnames , brain=brainvolcc , mass=masskg )
str(d)

plot(brain ~ mass, data =d)

d$mass_std <- (d$mass - mean(d$mass))/sd(d$mass)
d$brain_std <- d$brain/max(d$brain)

m7.1 <- quap(
  alist(
    brain_std ~ dnorm(mu, exp(log_sigma)),
    mu <- a + bM*mass_std,
    a ~ dnorm(0.5,1),
    bM ~ dnorm(0,1),
    log_sigma ~ dnorm(0,1)
  ), data = d
)
precis(m7.1)

set.seed(12)
s <- sim(m7.1)
r <- apply(s, 2, mean) - d$brain_std
resid_var <- var2(r)
outcome_var <- var2(d$brain_std)
1 - (resid_var/outcome_var)

R2_is_bad <- function(quap_fit){
  s <- sim(quap_fit, refresh = 0)
  r <- apply(s, 2, mean) - d$brain_std
  1 - var2(r)/var2(d$brain_std)
}

m7.2 <- quap(
  alist(
    brain_std ~ dnorm(mu, exp(log_sigma)),
    mu <- a + b[1]*mass_std + b[2]*mass_std^2,
    a ~ dnorm(0.5,1),
    b ~ dnorm(0,1),
    log_sigma ~ dnorm(0,1)
  ), data = d, start = list(b=rep(0,2))
)
m7.3 <- quap(
  alist(
    brain_std ~ dnorm(mu, exp(log_sigma)),
    mu <- a + b[1]*mass_std + b[2]*mass_std^2 + b[3]*mass_std^3,
    a ~ dnorm(0.5,1),
    b ~ dnorm(0,1),
    log_sigma ~ dnorm(0,1)
  ), data = d, start = list(b=rep(0,3))
)
m7.4 <- quap(
  alist(
    brain_std ~ dnorm( mu , exp(log_sigma) ),
    mu <- a + b[1]*mass_std + b[2]*mass_std^2 +
      b[3]*mass_std^3 + b[4]*mass_std^4,
    a ~ dnorm( 0.5 , 1 ),
    b ~ dnorm( 0 , 10 ),
    log_sigma ~ dnorm( 0 , 1 )
  ), data=d , start=list(b=rep(0,4)))
m7.5 <- quap(
  alist(
    brain_std ~ dnorm( mu , exp(log_sigma) ),
    mu <- a + b[1]*mass_std + b[2]*mass_std^2 +
      b[3]*mass_std^3 + b[4]*mass_std^4 +
      b[5]*mass_std^5,
    a ~ dnorm( 0.5 , 1 ),
    b ~ dnorm( 0 , 10 ),
    log_sigma ~ dnorm( 0 , 1 )
  ), data=d , start=list(b=rep(0,5)) )

m7.6 <- quap(
  alist(
    brain_std ~ dnorm(mu, 0.001),
    mu <- a + b[1]*mass_std + b[2]*mass_std^2 +
      b[3]*mass_std^3 + b[4]*mass_std^4 +
      b[5]*mass_std^5 + b[6]*brain_std^6,
    a ~ dnorm( 0.5 , 1 ),
    b ~ dnorm( 0 , 10 )
  ), data = d,
  start = list(b=rep(0,6))
)
a <- list(m7.1,m7.2, m7.3,m7.4,m7.5,m7.6)

for(i in 1:6){
  print(R2_is_bad(a[[i]]))
}

par(mfrow = c(2,3))
for(i in 1:5){
  post <- extract.samples(a[[i]])
  mass_seq <- seq( from=min(d$mass_std) , to=max(d$mass_std) , length.out=100 )
  l <- link( a[[i]] , data=list( mass_std=mass_seq ) )
  mu <- apply( l , 2 , mean )
  ci <- apply( l , 2 , PI )
  plot( brain_std ~ mass_std , data=d )
  lines( mass_seq , mu )
  shade( ci , mass_seq )
  mtext(paste("mass_std7.",i,": R_squared = ",round(R2_is_bad(a[[i]]),2)))
}

set.seed(1)
lppd(m7.1, n=1e4)

##Function for lppd####
set.seed(1)
logprob <- sim(m7.1, ll = T, n = 1e4)
n <- ncol(logprob)
ns <- nrow(logprob)
f <- function(i) log_sum_exp(logprob[,i]) - log(ns)
(lppd <- sapply(1:n, f))


###################
sapply(a,function(m) sum(lppd(m)))


N <- 20
kseq <- 1:5
dev <- sapply( kseq , function(k) {
  print(k);
  r <- replicate( 1e4 , sim_train_test( N=N, k=k ) );
  c( mean(r[1,]) , mean(r[2,]) , sd(r[1,]) , sd(r[2,]) )
} )

plot( 1:5 , dev[1,] , ylim=c( min(dev[1:2,])-5 , max(dev[1:2,])+10 ) ,
      xlim=c(1,5.1) , xlab="number of parameters" , ylab="deviance" ,
      pch=16 , col=rangi2 )
mtext( concat( "N = ",N ) )
points( (1:5)+0.1 , dev[2,] )
for ( i in kseq ) {
  pts_in <- dev[1,i] + c(-1,+1)*dev[3,i]
  pts_out <- dev[2,i] + c(-1,+1)*dev[4,i]
  lines( c(i,i) , pts_in , col=rangi2 )
  lines( c(i,i)+0.1 , pts_out )
}

library(rethinking)
###WAIC Calculations
data(cars)
m <- quap(
  alist(
    dist ~ dnorm(mu, sigma),
    mu <- a + b*speed,
    a ~ dnorm(0,100),
    b ~ dnorm(0, 10),
    sigma ~ dexp(1)
  ), data = cars
)
set.seed(94)
post <- extract.samples(m,n=1000)
n_samples <- 1000

logprop <- sapply(1:n_samples,
                  function(s){
                    mu <- post$a[s] + post$b[s]*cars$speed
                    dnorm(cars$dist,mu,post$sigma[s],log=T)
                  })

n_cases <- nrow(cars)
lppd <- sapply(1:n_cases,
               function(i){
                 log_sum_exp(logprop[i,])-log(n_samples)
               })


##Penalty term
pWAIC <- sapply(1:n_cases, function(i) var(logprop[i,]))
sum(pWAIC)

-2*(sum(lppd)-sum(pWAIC))

waic_vec <- -2*( lppd - pWAIC )
sqrt( n_cases*var(waic_vec) )
