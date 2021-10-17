
a <- data.frame(val = seq(0,1,0.001), 
                prior = dnorm(mean = 0.8, sd = 0.08, seq(0,1,0.001)),
                likelihood = dbinom(6,9,seq(0,1,0.001)))

a$posterior <- a$prior*a$likelihood
a$posterior <- a$posterior/sum(a$posterior)
for(i in 1:1001 ){
  if(i==1){
    a$probcum[i] <- a$posterior[i]
  }
  else {
    a$probcum[i] <- a$probcum[i-1]+a$posterior[i]
  }
}

par(mfrow= c(1,2))
plot(a$val, a$posterior, type = 'l', xlab = 'proportion of water', ylab = 'posterior probability')
plot(a$val, a$probcum, type = 'l', , xlab = 'proportion of water', ylab = 'Cumulative posterior probability')

mean(a$probcum)
a$diffmeancum <- abs(a$probcum- 0.5)
subset(a, a$diffmeancum == min(a$diffmeancum))

#-------------GRID APPROXIMATION-----------------#
#(1) Define the grid. This means you decide how many points to use in estimating the posterior, and then you make a list of the parameter values on the grid.
# (2) Compute the value of the prior at each parameter value on the grid.
# (3) Compute the likelihood at each parameter value.
# (4) Compute the unstandardized posterior at each parameter value, by multiplying the prior by the likelihood.
# (5) Finally, standardize the posterior, by dividing each value by the sum of all values.

#define grid
proportion_water <- seq(from = 0, to = 1, length.out = 100)

#define prior
prior <- rep(1, 100)
#compute likelihood
likelihood <- dbinom(6, size = 9, prob = proportion_water)
#compute product
unstd.posterior <- likelihood*prior

#standardize posterior
posterior <- unstd.posterior/sum(unstd.posterior)

plot(proportion_water, prior, type = 'l')
plot(proportion_water, posterior, type = 'l',
     xlab = 'probability of water',
     ylab = 'posterior probability')
mtext('20 points')


#------MAXIMUM A POSTERIORI-------
library(rethinking)
globe.qa <- quap(
  alist(
    w ~ dbinom(9, p), #binomial likelihood
    p ~ dunif(0,1) #uniform prior
  ),
  data=list(w=6) )

precis(globe.qa)
p <- extract.samples(globe.qa)
dens(p, norm.comp = T)

#analytical calculation
w <-6
n <- 9
par(mfrow = c(1,2))
curve(dbeta(x,w+1,n-w+1), from = 0, to = 1)
curve(dnorm(x, 0.667, 0.16), lty =2, add = T)


#________CHAPTER 3 
##Sampling from grid approximation
proportion_water <- seq(from = 0, to = 1, length.out = 1000)
prior <- rep(1, 1000)
likelihood <- dbinom(6, size = 9, prob = proportion_water)
unstd.posterior <- likelihood*prior
posterior <- unstd.posterior/sum(unstd.posterior)

samples <- sample(proportion_water, prob=posterior, size=1e4, replace = T)
par(mfrow = c(1,2))
plot(samples)
dens(samples)


#interval of defined boundary
# Find the posterior probability that the proportion of wayer is less than 0.5
#add up the posterior prob where p<0.5
sum(posterior[proportion_water < .5])
# or
sum(samples < 0.5)/1e4

# What is the probability that proportion of water lies between 0.5 and 0.75
sum(samples >0.5 & samples < 0.75)/1e4

#credibility interval
quantile(samples,0.8)
quantile(samples, c(0,0.5))

#PI gives you the confidence interval, where first argument is the sample
#and the second argument is the amoung of confidence you want to have
#in the result. This is done by dividing the non-confidence regeion equally in both tails
PI(samples, prob=0.90)

#in contrast, the HDPI function displays the 50% highest posterior
# density interval (HPDI).51 The HPDI is the narrowest interval containing the specified
# probability mass.
HPDI(samples, prob=0.9)

#LOSS FUNCTION
sum(posterior*abs(0.5-proportion_water))
loss <- sapply(proportion_water, function(d) sum(posterior*abs(d-proportion_water)))
plot(loss, type= 'l')
proportion_water[which.min(loss)]


##Simulating data using baysian
#likelihood
dbinom(0:2, size=2, prob=0.7)
##sampling
rbinom(12, size= 2, prob = 0.7)


#generate 100k dummy observations just t verify that each value (0,1,2)
#appears in the proportion to its likelihood

#--this simulates 2 globe tosses 100k times
dummy_w <- rbinom(1e5, size =2, prob = 0.7)
prop.table(table(dummy_w))

#now lets simulate 9 tosses 100k times
dummy_w <- rbinom(10, size =9, prob = 0.7)
simplehist(dummy_w, xlab = 'dummy water count')
prop.table(table(dummy_w))

# ## For each possible value of the parameter p, there is an implied
# distribution of outcomes. So if you were to compute the sampling distribution of outcomes at
# each value of p, then you could average all of these prediction distributions together, using the
# posterior probabilities of each value of p, to get a posterior predictive distribution.

#STEP 1: CREATE POSTERIOR BASED ON OBSERVATIONS
proportion_water <- seq(from = 0, to = 1, length.out = 1000)
prior <- rep(1, 1000)
likelihood <- dbinom(6, size = 9, prob = proportion_water)
unstd.posterior <- likelihood*prior
posterior <- unstd.posterior/sum(unstd.posterior)

#SAMPLE FROM THE POSTERIOR DIFFERENT VALUES OF
#PROB OF SUCCESS
samples <- sample(proportion_water, prob=posterior, size=1e4, replace = T)

#uSE THE SAMPLE VALUES OF P FROM THE POSTERIOR TO GENERATE MORE SAMPLES
w <- rbinom(1e4, size = 9, prob = samples)
simplehist(w)
plot(prop.table(table(w)), type= 'l')
plot(posterior, type = 'l')


#########    PRACTICE QUESTIONS     ######################

p_grid <- seq( from=0 , to=1 , length.out=1000 )
prior <- rep( 1 , 1000 )
likelihood <- dbinom( 6 , size=9 , prob=p_grid )
posterior <- likelihood * prior
posterior <- posterior / sum(posterior)
set.seed(100)
samples <- sample( p_grid , prob=posterior , size=1e4 , replace=TRUE )

#posterior probability lies below 0.2
sum(samples<0.2)/1e4
#posterior that lies above 0.8
sum(samples>0.8)/1e4
#b/2 0.2 and 0.8
sum(samples>0.2 & samples <0.8)/1e4

#20% of posterior probability lies below which value of p
quantile(samples, 0.2)
#20% of posterior lies above which value of p
quantile(samples, 0.8)

#which value of p contains the narrowest interval qual to 0.66 of posterior
HPDI(samples, 0.66)
# 
# Which values of p contain 66% of the posterior probability,
# assuming equal posterior probability both below and above the interval?
PI(samples, 0.66)


###Medium questions
prior <- rep(1, 1000)
p_grid <- seq(from = 0, to = 1, length.out = 1000)
likelihood <- dbinom(8, size= 15, prob = p_grid)
posterior <- likelihood*prior
posterior <- posterior/sum(posterior)
plot(posterior, type = 'l')

samples <- sample(p_grid, size = 1e4, replace = T, prob = posterior )
hist(samples)
HPDI(samples, 0.9)

w <- rbinom(1e4,size =15, prob = samples)
mean(w==8)
plot(posterior, type = 'l')

#Using the posterior distribution constructed from the new (8/15) data,
# now calculate the probability of observing 6 water in 9 tosses.
w <- rbinom(1e4, size= 9, prob = samples)
mean(w==6)


prior <- rep(1, 1000)
prior[1:499] <- 0
plot(prior, type = 'l')
p_grid <- seq(from = 0, to = 1, length.out = 1000)
likelihood <- dbinom(8, size= 15, prob = p_grid)
posterior <- likelihood*prior
posterior <- posterior/sum(posterior)
plot(posterior, type = 'l')

samples <- sample(p_grid, size = 1e4, replace = T, prob = posterior )
HPDI(samples, 0.9)

w <- rbinom(1e4,size =15, prob = samples)
mean(w==8)

#Using the posterior distribution constructed from the new (8/15) data,
# now calculate the probability of observing 6 water in 9 tosses.
w <- rbinom(1e4, size= 9, prob = samples)
mean(w==6)

