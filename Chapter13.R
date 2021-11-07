rm(list = ls())
library(rethinking)
data("reedfrogs")
d <- reedfrogs
str(d)
d

# A multilevel model, in which we simultaneously estimate both an intercept for each tank
# and the variation among tanks, is what we want. This will be a varying intercepts model.
# Varying intercepts are the simplest kind of varying effects.

##Heres an example without using a multilevel model, just a regular one

d$tank <- 1:nrow(d)

dat <- list(
  S = d$surv,
  N = d$density,
  tank = d$tank
)

m13.1 <- ulam(
  alist(
    S ~ dbinom(N, p),
    logit(p) <- a[tank],
    a[tank] ~ dnorm(0,1.5)
  ), data = dat, chains = 4, log_lik = T
)

precis(m13.1,2)

##now we do a multilevel model for the same problem
# All that us required to enable adaptive pooling is to make the prior for the 'a' parameter
# a function of some new paramters

m13.2 <- ulam(
  alist(
    S ~ dbinom(N, p),
    logit(p) <- a[tank],
    a[tank] ~ dnorm(a_bar, sigma),
    a_bar ~ normal(0, 1.5),
    sigma ~ dexp(1)
  ), data = dat, chains = 4, log_lik = T
)

precis(m13.2)
compare(m13.1, m13.2)

par(mfrow = c(2,1))

##visualise the two posteriors to compare them

#extract stan samples
post <- extract.samples(m13.2)

##compute mean estimates for each tank
#also transform to probability with logistic

d$propsurv.est <- logistic(apply(post$a,2,mean))

##display raw proportions surviving in each tank
plot(d$propsurv, ylim = c(0,1), pch = 16, xaxt = 'n',
     xlab = 'tank', ylab = 'proportion survival', col = rangi2)
axis(1, at = c(1,16,32,48), labels = c(1,16,32, 48))

##overlaw posterior means
points(d$propsurv.est, col = col.alpha('black', 0.5), pch = 16 )

#mark posterior mean probability across tanks
abline(h=mean(inv_logit(post$a_bar)),lty=2)

#draw vertical dividers between tank densities
abline(v=16.5, lwd=0.5)
abline(v=32.5, lwd=0.5)
text(8,0,'small tanks')
text(16+8,0,'medium tanks')
text(32+8,0,'large tanks')
for(i in 1:length(d$propsurv.est)){
  lines(x=c(i,i), y =c(d$propsurv[i], d$propsurv.est[i]))
}

post <- extract.samples(m13.1)

##compute mean estimates for each tank
#also transform to probability with logistic

d$propsurv.est <- logistic(apply(post$a,2,mean))

##display raw proportions surviving in each tank
plot(d$propsurv, ylim = c(0,1), pch = 16, xaxt = 'n',
     xlab = 'tank', ylab = 'proportion survival', col = rangi2)
axis(1, at = c(1,16,32,48), labels = c(1,16,32, 48))

##overlaw posterior means
points(d$propsurv.est, col = col.alpha('black', 0.9), pch = 16 )

#mark posterior mean probability across tanks
abline(h=0.789,lty=2)

#draw vertical dividers between tank densities
abline(v=16.5, lwd=0.5)
abline(v=32.5, lwd=0.5)
text(8,0,'small tanks')
text(16+8,0,'medium tanks')
text(32+8,0,'large tanks')


# We can visualize
# it by sampling from the posterior distribution, as usual. First we'll plot 100 Gaussian distributions, 
# one for each of the first 100 samples from the posterior distribution of both ??
# and ??. Then we'll sample 8000 new log-odds of survival for individual tanks. The result will
# be a posterior distribution of variation in survival in the population of tanks. Before we do
# the sampling though, remember that "sampling" from a posterior distribution is not a simulation
# of empirical sampling. It's just a convenient way to characterize and work with the
# uncertainty in the distribution. 
par(mfrow = c(1,1))
##show first 100 populations in the posterior
post <- extract.samples(m13.2)
plot(NULL, xlim = c(-3,4), ylim = c(0,0.35),
     xlab = 'log-odds survive', ylab = 'Density')
for ( i in 1:100 ){
  curve( dnorm(x,post$a_bar[i],post$sigma[i]) , add=TRUE ,
         col=col.alpha("black",0.2) )
}
  

##sample 8000 imaginary tanks from the posterior distribution
sim_tanks <- rnorm(8000, post$a_bar, post$sigma)

##transform to probabilit and visualise
dens(inv_logit(sim_tanks), lwd = 2, adj = 0.1,
     xlab = 'probability survive')



##Here we are going to simulate the data and fit our model on the simulated data to
# make sure that our model is running correctly

a_bar <- 1.5
sigma <- 1.5
nponds <- 60
Ni <- as.integer(rep(c(5,10,25,35), each = 15))
set.seed(5005)
a_pond <- rnorm(nponds, mean = a_bar, sd = sigma)
dsim <- data.frame(pond = 1:nponds, Ni=Ni, true_a=a_pond)
dsim

##simulating survivors
dsim$Si <- rbinom(nponds, prob = logistic(dsim$true_a), size = dsim$Ni)


##There are 3 types of models we are trying to compare here
#The no pooling model, where an intercept is calculated for each pond seperately
# The one pool model, where the entire population is used for calculating an estimate of population
# The partial pooling model, where an intercept is calculated for each pond, but the estimate for each pond
# becomes a prior for the next pond that needs to be estimated. This becomes a multilevel model

# Computing by the no pooling model
dsim$p_nopool <- dsim$Si/dsim$Ni

# Compute the partial pooling estimates
dat <- list(Si=dsim$Si, Ni=dsim$Ni, pond=dsim$pond)
dat
m13.3 <- ulam(
  alist(
    Si ~ dbinom(Ni, p),
    logit(p) <- a_pond[pond],
    a_pond[pond] ~ dnorm(a_bar, sigma),
    a_bar ~ dnorm(0, 1.5),
    sigma ~ dexp(1)
  ), data = dat, chains = 4
)

# m13.3 <- stan( fit=m13.3@stanfit , data=dat , chains=4 )

precis(m13.3)

# Now let's compute the predicted survival proportions and add those proportions to our
# growing simulation data frame. To indicate that it contains the partial pooling estimates, I'll
# call the column p_partpool.
post <- extract.samples(m13.3)
dsim$p_partpool <- apply(inv_logit(post$a_pond),2,mean)
dsim$p_true <- inv_logit(dsim$true_a)
nopool_error <- abs(dsim$p_nopool - dsim$p_true)
partpool_error <- abs(dsim$p_partpool - dsim$p_true)

par(mfrow = c(1,1))
##Plot the graph now to get the point of the entire lesson
plot(1:60, nopool_error, xlab = 'pond', ylab = 'absolute error',
     col=rangi2, pch = 16, 
     ylim = c(0,max(max(nopool_error), max(partpool_error))))
points(1:60, partpool_error, pch = 16)
for(i in 1:length(partpool_error)){
  lines(x=c(i,i), y= c(nopool_error[i], partpool_error[i]))
}

abline(v=15, lwd=0.5)
abline(v=30, lwd=0.5)
abline(v=45, lwd=0.5)

text(7.5,0,'tiny ponds(5)', cex = 0.7)
text(22.5,0,'small ponds(10)', cex = 0.7)
text(37.5,0,'medium ponds(25)', cex = 0.7)
text(52.5,0,'large ponds(35)', cex = 0.7)

nopool_avg <- aggregate(nopool_error,list(dsim$Ni),mean)
partpool_avg <- aggregate(partpool_error,list(dsim$Ni),mean)

lines(x=1:15, y = rep(nopool_avg[1,2],15), col = rangi2, lwd = 2)
lines(x=1:15, y = rep(partpool_avg[1,2],15), lty=2, lwd = 2)
lines(x=16:30, y = rep(nopool_avg[2,2],15), col = rangi2, lwd = 2)
lines(x=16:30, y = rep(partpool_avg[2,2],15), lty=2, lwd = 2)
lines(x=31:45, y = rep(nopool_avg[3,2],15), col = rangi2, lwd = 2)
lines(x=31:45, y = rep(partpool_avg[3,2],15), lty=2, lwd = 2)
lines(x=46:60, y = rep(nopool_avg[4,2],15), col = rangi2, lwd = 2)
lines(x=46:60, y = rep(partpool_avg[4,2],15), lty=2, lwd = 2)




##Without having to compile the model again and fit it on new data,
# you can use the already compiled data to fit new data
a <- 1.5
sigma <- 1.5
nponds <- 60
Ni <- as.integer( rep( c(5,10,25,35) , each=15 ) )
a_pond <- rnorm( nponds , mean=a , sd=sigma )
dsim <- data.frame( pond=1:nponds , Ni=Ni , true_a=a_pond )
dsim$Si <- rbinom( nponds,prob=inv_logit( dsim$true_a ),size=dsim$Ni )
dsim$p_nopool <- dsim$Si / dsim$Ni
newdat <- list(Si=dsim$Si,Ni=dsim$Ni,pond=1:nponds)
m13.3new <- stan( fit=m13.3@stanfit , data=newdat , chains=4 )

post <- extract.samples( m13.3new )
dsim$p_partpool <- apply( inv_logit(post$a_pond) , 2 , mean )
dsim$p_true <- inv_logit( dsim$true_a )
nopool_error <- abs( dsim$p_nopool - dsim$p_true )
partpool_error <- abs( dsim$p_partpool - dsim$p_true )
plot( 1:60 , nopool_error , xlab="pond" , ylab="absolute error" , col=rangi2 , pch=16 )
points( 1:60 , partpool_error )


abline(v=15, lwd=0.5)
abline(v=30, lwd=0.5)
abline(v=45, lwd=0.5)
text(7.5,0.38,'tiny ponds(5)')
text(22.5,0.38,'small ponds(10)')
text(37.5,0.38,'medium ponds(25)')
text(52.5,0.37,'large ponds(35)')
lines(x=1:15, y = rep(mean(nopool_error[1:15]),15), col = rangi2, lwd = 2)
lines(x=1:15, y = rep(mean(partpool_error[1:15]),15), lty=2, lwd = 2)
lines(x=15:30, y = rep(mean(nopool_error[16:30]),16), col = rangi2, lwd = 2)
lines(x=15:30, y = rep(mean(partpool_error[16:30]),16), lty=2, lwd = 2)
lines(x=30:45, y = rep(mean(nopool_error[31:45]),16), col = rangi2, lwd = 2)
lines(x=30:45, y = rep(mean(partpool_error[31:45]),16), lty=2, lwd = 2)
lines(x=45:60, y = rep(mean(nopool_error[46:60]),16), col = rangi2, lwd = 2)
lines(x=45:60, y = rep(mean(partpool_error[46:60]),16), lty=2, lwd = 2)




##MORE THAN ONE TYPE OF CLUSTER
library(rethinking)
data("chimpanzees")
d <- chimpanzees
d$treatment <- 1 + d$prosoc_left + 2*d$condition

dat_list <- list(
  pulled_left = d$pulled_left,
  actor = d$actor,
  block_id = d$block,
  treatment = as.integer(d$treatment)
)

set.seed(13)
m13.4 <- ulam(
  alist(
    pulled_left ~ dbinom( 1 , p ) ,
    logit(p) <- a[actor] + g[block_id] + b[treatment] ,
    b[treatment] ~ dnorm( 0 , 0.5 ),
    ## adaptive priors
    a[actor] ~ dnorm( a_bar , sigma_a ),
    g[block_id] ~ dnorm( 0 , sigma_g ),
    ## hyper-priors
    a_bar ~ dnorm( 0 , 1.5 ),
    sigma_a ~ dexp(1),
    sigma_g ~ dexp(1)
  ) , data=dat_list , chains=4 , cores=4 , log_lik=TRUE )

precis(m13.4, 2)
plot(NULL, xlim = c(0,5), ylim = c(0,4))
post <- extract.samples(m13.4)
dens(post$sigma_a, add = T, col = rangi2)
dens(post$sigma_g, add = T)

set.seed(14)
m13.5 <- ulam(
  alist(
    pulled_left ~ dbinom( 1 , p ) ,
    logit(p) <- a[actor] + b[treatment] ,
    b[treatment] ~ dnorm( 0 , 0.5 ),
    a[actor] ~ dnorm( a_bar , sigma_a ),
    a_bar ~ dnorm( 0 , 1.5 ),
    sigma_a ~ dexp(1)
  ) , data=dat_list , chains=4 , cores=4 , log_lik=TRUE )

compare(m13.4, m13.5)

set.seed(16)
m13.6 <- ulam(
  alist(
    pulled_left ~ dnorm(1, p),
    logit(p) <- a[actor] + g[block_id] + b[treatment],
    b[treatment] ~ dnorm(0, sigma_b),
    a[actor] ~ dnorm(a_bar, sigma_a),
    g[block_id] ~ dnorm(0, sigma_g),
    a_bar ~ dnorm(0, 1.5),
    sigma_a ~ dexp(1),
    sigma_g ~ dexp(1),
    sigma_b ~ dexp(1)
  ), data = dat_list, chains =4, cores = 4, log_lik = T
)

coeftab(m13.4, m13.6)


## the devils funnel
m13.7 <- ulam(
  alist(
    v ~ normal(0,3),
    x ~ normal(0, exp(v))
  ), data = list(N=1), chains = 4
)


precis(m13.7)
traceplot(m13.7)


m13.7nc <- ulam(
  alist(
    v ~ normal(0, 3),
    z ~ normal(0,1),
    gq> real[1]:x <<- z*exp(v)
  ), data= list(N=1), chains = 4
)

precis(m13.7nc)
trankplot(m13.7nc, lwd = 2)


##Non centered chimpanzees
#First thing to avoid divergent transitions is to increase the acceptance rate of the
# metropolis proposal. This can be done by changing the adapt_delta parameter in ulam. The 
# default is 0.95. High adapt delta results in smaller leapfrog steps, leading to higher acceptance rate
set.seed(13)
m13.4b <- ulam(m13.4, chains = 4, cores = 4, control = list(adapt_delta = 0.99))
divergent(m13.4b)
precis(m13.4b)


# model 13.4:
# pulled_left ~ dbinom(1, p)
# logit(p) <- a[actor] + g[block_id] + b[treatment]
# b[treatment] ~ dnorm(0, 0.5)
# a[actor] ~ dnorm(a_bar, sigma_a)
# g[block_id] ~ dnorm(0, sigma_g)
# a_bar ~ dnorm(0, 1.5)
# sigma_a ~ dexp(1)
# sigma_g ~ dexp(1)

set.seed(13)
m13.4nc <- ulam(
  alist(
    pulled_left ~ dbinom(1, p),
    logit(p) <- a_bar + z[actor]*sigma_a + 
      x[block_id]*sigma_g + b[treatment],
    b[treatment] ~ dnorm(0, 0.5),
    z[actor] ~ dnorm(0,1),
    x[block_id]~ dnorm(0,1),
    a_bar ~ dnorm(0, 1.5),
    sigma_a ~ dexp(1),
    sigma_g ~ dexp(1),
    gq> vector[actor]:a <<-a_bar + z*sigma_a,
    gq>vector[block_id]:g <<- x*sigma_g
  ), data=  dat_list, chains = 4, cores = 4
)

precis_c <- precis(m13.4, depth=2)
precis_nc <- precis(m13.4nc,2)
pars <- c(paste('a[',1:7,']', sep = ''),paste('g[',1:6,']', sep = ''),
          paste('b[',1:4,']', sep = ''), 'a_bar', 'sigma_a','sigma_g' )
neff_table <- cbind(precis_c[pars,"n_eff"], precis_nc[pars, 'n_eff'])
par(mfrow = c(1,1))
plot(neff_table, xlim=range(neff_table), ylim=range(neff_table),
     xlab = 'n_eff (centered)', ylab = 'n_eff (non-centered)', lwd = 2)
abline(a=0,b=1, lty=2)


##MULTILEVEL POSTERIOR PREDICTIONS

##Making predictions in and out of sample

chimp <- 2
d_pred <- list(
  actor = rep(chimp,4),
  treatment = 1:4,
  block_id = rep(1,4)
)
p <- link(m13.4, data = d_pred)
p_mu <- apply(p,2,mean)
p_ci <- apply(p,2,PI)


post <- extract.samples(m13.4)
str(post)
dens(post$a[,5])

p_link <- function(treatment, actor = 1, block_id =1){
  logodds <- with(post,
                  a[,actor]+g[,block_id]+b[,treatment])
  return(inv_logit(logodds))
}

p_raw <- sapply(1:4, function(i){
  p_link(i,actor=4,block_id = 1)
})

p_mu <- apply(p_raw,2,mean)
p_ci <- apply(p_raw,2,PI)
plot(NULL,xlim = c(1,4), ylim = c(0,1))
lines(1:4, p_mu)
shade(p_ci, 1:4)



# Posterior predictions out of sample now

p_link_abar <- function(treatment){
  logodds <- with(post,
                  a_bar + b[,treatment])
  return(inv_logit(logodds))
}

post <- extract.samples(m13.4)
p_raw <- sapply(1:4, function(i)p_link_abar(i))
p_mu <- apply(p_raw,2,mean)
p_ci <- apply(p_raw,2,PI)
plot(NULL, xlab = 'treatment',
     ylab = 'proportion pulled left',
     ylim = c(0,1), xlim = c(1,4), xaxt = 'n')
axis(1, at = 1:4, labels = c('R/N', 'L/N', 'R/P', 'L/P'))
lines(1:4, p_mu)
shade(p_ci, 1:4)

a_sim <- with(post, rnorm(length(post$a_bar), a_bar, sigma_a))
p_link_asim <- function(treatment){
  logodds <- with(post, a_sim + b[,treatment])
  return(inv_logit(logodds))
}
p_raw_asim <- sapply(1:4, function(i) p_link_asim(i))

p_mu <- apply(p_raw_asim,2,mean)
p_ci <- apply(p_raw_asim,2,PI)
plot(NULL, xlab = 'treatment',
     ylab = 'proportion pulled left',
     ylim = c(0,1), xlim = c(1,4), xaxt = 'n')
axis(1, at = 1:4, labels = c('R/N', 'L/N', 'R/P', 'L/P'))
lines(1:4, p_mu, lwd =2)
shade(p_ci,1:4, alpha = 0.8)


a_sim <- with(post, rnorm(length(post$a_bar), a_bar, sigma_a))
p_link <- function(treatment){
  logodds <- with(post, a_sim + b[,treatment])
  return(inv_logit(logodds))
}
p_raw_asim <- sapply(1:4,function(i) p_link_asim(i))
plot( NULL , xlab="treatment" , ylab="proportion pulled left" ,
      ylim=c(0,1) , xaxt="n" , xlim=c(1,4) )
axis( 1 , at=1:4 , labels=c("R/N","L/N","R/P","L/P") )
for ( i in sample(1:2000,100)) lines( 1:4 , p_raw_asim[i,] , col=col.alpha('black',0.1) , lwd=2 )

##Post stratification